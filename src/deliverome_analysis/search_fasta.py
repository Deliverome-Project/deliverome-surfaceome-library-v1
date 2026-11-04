"""Build, check and publish a versioned DIA search FASTA: UniProt human + contaminants + iRT.

A search is only reproducible if the database it ran against is pinned. Spectronaut records
the FASTA's file name, not where it came from, so the Astral benchmark ran against a 2018
file of unknown origin that never reached S3. This module builds the database from three
pinned sources, writes a manifest (release, URL, commit, SHA-256 and counts for every part)
and publishes it under ``external://search-fasta/<name>/`` so every search can name the
exact snapshot it used:

* **UniProt human** (UP000005640), reviewed (Swiss-Prot), with isoforms by default; the
  release number comes from the ``X-UniProt-Release`` response header.
* **Contaminants**: the Frankenfield et al., 2022 (J Proteome Res) libraries, pinned to a git
  commit: by default the universal library (which the authors recommend for every bottom-up
  analysis) plus their cell-culture library (medium and serum proteins), since our samples
  are cultured cells. Duplicates across libraries are dropped, as are entries whose sequence
  is identical to a human entry, so the human accession owns those peptides.
* **iRT standard**: the Biognosys iRT kit peptides, read from the kit's reference sheet (the
  sheet is archived on S3 too, so sequences are never typed by hand), joined into one
  synthetic fusion entry ``sp|Biognosys_iRT|IRT_KIT`` whose tryptic digest is exactly the
  kit's peptides.

Usage:

    uv run python -m deliverome_analysis.search_fasta build --irt-sheet irt-kit-reference-sheet.xls --out search_db/
    uv run python -m deliverome_analysis.search_fasta publish search_db/

UniProt's REST API always serves its current release, so re-running ``build`` later gives a
newer database by design; the published snapshot is what makes an old search reproducible.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
HUMAN_PROTEOME = "UP000005640"
CONTAMINANT_REPO = "HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics"
CONTAMINANT_COMMIT = "bddf6da9334f4300caed75f7a87c096dd96b6bca"  # main on 2026-10-07
CONTAMINANT_FILES = {
    "universal": "Universal protein contaminant FASTA/0602_Universal Contaminants.fasta",
    "cell_culture": "Sample-type specific contaminant FASTA/0602_Cell Culture Contaminants.fasta",
}
CONTAMINANT_CITATION = "Frankenfield et al., 2022, J Proteome Res 21:2104, PMID 35793413"
DEFAULT_CONTAMINANTS = ("universal", "cell_culture")
# The kit's reference sheet as Biognosys distributes it (last saved 2015-06-23).
IRT_KIT_KEY = "biognosys/irt-kit/2015-06-23/irt-kit-reference-sheet.xls"
IRT_ACCESSION = "Biognosys_iRT"
HUMAN_TAXON = "OX=9606"


# ----------------------------------------------------------------------------------------------
# FASTA helpers
# ----------------------------------------------------------------------------------------------


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """(header without '>', sequence) pairs."""
    records, header, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def write_fasta(records: list[tuple[str, str]], path: str | os.PathLike, width: int = 60) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def accession(header: str) -> str:
    """``sp|P02786|TFR1_HUMAN ...`` -> ``P02786``; other headers -> first token."""
    first = header.split()[0]
    parts = first.split("|")
    return parts[1] if len(parts) >= 3 else first


def tryptic_peptides(seq: str) -> list[str]:
    """Full trypsin digest: cleave after K/R unless followed by P."""
    return [p for p in re.split(r"(?<=[KR])(?!P)", seq) if p]


def shared_peptides(
    contaminants: list[tuple[str, str]],
    human: list[tuple[str, str]],
    *,
    min_len: int = 7,
    max_len: int = 30,
) -> pd.DataFrame:
    """Kept contaminants that share fully tryptic peptides with human entries (I = L).

    Such peptides cannot be assigned to the human protein or the contaminant from the MS data
    alone, so the human protein's grouping and quantification can change when the
    contaminant is searched. One row per contaminant: shared / total peptides and the human
    accessions (isoform suffix removed) it overlaps.
    """

    def peps(seq):
        return {p.replace("I", "L") for p in tryptic_peptides(seq) if min_len <= len(p) <= max_len}

    by_pep: dict[str, set[str]] = {}
    for h, sq in human:
        for p in peps(sq):
            by_pep.setdefault(p, set()).add(accession(h).split("-")[0])
    rows = []
    for h, sq in contaminants:
        own = peps(sq)
        hit = [p for p in own if p in by_pep]
        if hit:
            accs = sorted(set().union(*(by_pep[p] for p in hit)))
            rows.append((accession(h), h.split()[0], len(hit), len(own), ";".join(accs)))
    return pd.DataFrame(
        rows,
        columns=["contaminant", "entry", "shared_peptides", "total_peptides", "human_accessions"],
    )


def _rejoinable(pieces: list[str]) -> tuple[list[str], int]:
    """Order tryptic pieces so their concatenation digests back to exactly these pieces.

    Two constraints: a piece starting with P cannot be cleaved from a preceding K/R, so at most
    one such piece can be kept, at the front; a piece not ending in K/R (the C-terminus) must
    be last. Pieces that cannot be placed are dropped and counted.
    """
    p_start = [p for p in pieces if p.startswith("P")]
    rest = [p for p in pieces if not p.startswith("P")]
    tail = [p for p in rest if p[-1] not in "KR"]
    body = [p for p in rest if p[-1] in "KR"]
    head = [p for p in p_start if p[-1] in "KR"][:1] or (
        [p_start[0]] if p_start and not body and not tail else []
    )
    ordered = head + body + tail[-1:]
    if tryptic_peptides("".join(ordered)) != ordered:
        raise ValueError("re-joined contaminant does not digest back to its pieces")
    return ordered, len(pieces) - len(ordered)


def trim_shared(
    contaminants: list[tuple[str, str]],
    protected: list[tuple[str, str]],
    *,
    min_len: int = 7,
    max_len: int = 30,
) -> tuple[list[tuple[str, str]], list[dict]]:
    """Cut out of each contaminant the tryptic peptides it shares with ``protected`` entries.

    The remaining tryptic pieces are joined in their original order under the same accession,
    so trypsin regenerates them (checked); an entry left with no searchable peptide is
    dropped. Returns the new contaminant list and one record per changed entry.
    """

    def key(p):
        return p.replace("I", "L")

    prot_peps: dict[str, set[str]] = {}
    for h, sq in protected:
        for p in tryptic_peptides(sq):
            if min_len <= len(p) <= max_len:
                prot_peps.setdefault(key(p), set()).add(accession(h).split("-")[0])
    out, changes = [], []
    for h, sq in contaminants:
        pieces = tryptic_peptides(sq)
        shared = [p for p in pieces if min_len <= len(p) <= max_len and key(p) in prot_peps]
        if not shared:
            out.append((h, sq))
            continue
        drop = set(shared)
        keep, lost = _rejoinable([p for p in pieces if p not in drop])
        new = "".join(keep)
        searchable = [p for p in tryptic_peptides(new) if min_len <= len(p) <= max_len]
        record = {
            "accession": accession(h),
            "removed_peptides": sorted(drop),
            "human_accessions": sorted(set().union(*(prot_peps[key(p)] for p in drop))),
            "length_before": len(sq),
            "length_after": len(new),
            "searchable_peptides_left": len(searchable),
            "unshared_pieces_lost_to_rejoining": lost,
        }
        if searchable:
            note = " (shared peptides with protected human proteins removed)"
            out.append(
                (h.split(" OS=")[0] + note + " OS=" + h.split(" OS=", 1)[1], new)
                if " OS=" in h
                else (h + note, new)
            )
        else:
            record["dropped"] = True
        changes.append(record)
    return out, changes


def sha256(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------
# sources
# ----------------------------------------------------------------------------------------------


def read_irt_kit(path: str | os.PathLike) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Peptides (sequence, kit_id, irt, q1_mz) and transitions from the Biognosys kit sheet."""
    tr = pd.read_excel(path, sheet_name="iRT Kit transitions")
    tr = tr.rename(
        columns={
            "Q1 monoisotopic": "q1_mz",
            "Q3": "q3_mz",
            "relative intensity (approximate, TSQ-Vantage)": "relative_intensity",
            "precursor charge": "precursor_charge",
            "fragment type": "fragment_type",
            "fragment charge": "fragment_charge",
            "fragment number": "fragment_number",
            "nominal sequence": "sequence",
            "sequence id": "kit_id",
            "transition id": "transition_id",
        }
    )
    rt = pd.read_excel(path, sheet_name="RT->iRT->RT")
    rt = rt.rename(columns={"standard": "sequence", "Q1": "q1_mz", "iRT": "irt"})
    rt = rt.dropna(subset=["sequence", "irt"])
    rt = rt[rt["sequence"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")]
    ids = tr.drop_duplicates("sequence").set_index("sequence")["kit_id"]
    peps = rt[["sequence", "irt", "q1_mz"]].copy()
    peps.insert(1, "kit_id", peps["sequence"].map(ids))
    if peps["kit_id"].isna().any() or len(peps) != tr["sequence"].nunique():
        raise ValueError("iRT sheet: peptide list and transition list disagree")
    return peps.reset_index(drop=True), tr


def irt_fusion_record(peptides: list[str]) -> tuple[str, str]:
    """One synthetic entry whose tryptic digest is exactly ``peptides`` (in order)."""
    seq = "".join(peptides)
    if tryptic_peptides(seq) != list(peptides):
        raise ValueError("iRT peptides would not digest back to themselves when fused")
    header = (
        f"sp|{IRT_ACCESSION}|IRT_KIT Biognosys iRT kit peptides (synthetic fusion of "
        f"{len(peptides)} standards) OS=Synthetic construct OX=32630 GN=iRT PE=5 SV=1"
    )
    return header, seq


def _get(url: str) -> tuple[str, dict]:
    req = urllib.request.Request(url, headers={"User-Agent": "deliverome-analysis/search_fasta"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        return resp.read().decode(), {k.lower(): v for k, v in resp.headers.items()}


def fetch_uniprot(
    proteome: str = HUMAN_PROTEOME, *, reviewed: bool = True, isoforms: bool = True
) -> tuple[str, dict]:
    """FASTA text for a UniProt proteome plus provenance (release, date, URL, query)."""
    query = f"(proteome:{proteome})" + (" AND (reviewed:true)" if reviewed else "")
    params = {"format": "fasta", "query": query}
    if isoforms:
        params["includeIsoform"] = "true"
    url = f"{UNIPROT_STREAM}?{urllib.parse.urlencode(params)}"
    text, headers = _get(url)
    return text, {
        "source": "UniProtKB",
        "release": headers.get("x-uniprot-release"),
        "release_date": headers.get("x-uniprot-release-date"),
        "url": url,
        "query": query,
        "isoforms": isoforms,
    }


def fetch_contaminants(
    library: str = "universal", commit: str = CONTAMINANT_COMMIT
) -> tuple[str, dict]:
    path = CONTAMINANT_FILES[library]
    url = (
        f"https://raw.githubusercontent.com/{CONTAMINANT_REPO}/{commit}/{urllib.parse.quote(path)}"
    )
    text, _ = _get(url)
    return text, {
        "source": CONTAMINANT_CITATION,
        "repo": CONTAMINANT_REPO,
        "commit": commit,
        "file": path,
        "url": url,
    }


# ----------------------------------------------------------------------------------------------
# build / publish
# ----------------------------------------------------------------------------------------------


def database_name(release: str, isoforms: bool, contaminants: tuple[str, ...]) -> str:
    """e.g. ``human-sp-iso_cont-universal-and-cell-culture_irt/2026_03`` ({source}/{version}).

    Only [a-z0-9._-] are allowed in external keys, hence ``-and-`` rather than ``+``.
    """
    iso = "sp-iso" if isoforms else "sp"
    cont = "-and-".join(c.replace("_", "-") for c in contaminants)
    return f"human-{iso}_cont-{cont}_irt/{release}"


def build_search_fasta(
    out_dir: str | os.PathLike,
    irt_sheet: str | os.PathLike,
    *,
    isoforms: bool = True,
    contaminants: tuple[str, ...] = DEFAULT_CONTAMINANTS,
    contaminant_commit: str = CONTAMINANT_COMMIT,
    uniprot_text: str | None = None,
    uniprot_meta: dict | None = None,
    contaminant_texts: dict[str, str] | None = None,
    protect: dict[str, str] | None = None,
    protect_source: str | None = None,
) -> dict:
    """Write search.fasta, its three parts, the iRT tables and manifest.json; return the manifest.

    The ``*_text`` / ``*_meta`` arguments let tests (or an offline rebuild from an archived
    snapshot) supply sources instead of downloading them.

    ``protect`` maps human accession -> reason (e.g. canonical surface tier). A contaminant may
    not claim any tryptic peptide (7-30 aa, I = L) of a protected protein: such peptides are cut
    out of the contaminant entry (see :func:`trim_shared`), so the search gives them to the human
    protein, and the contaminant is identified by its own peptides only.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if uniprot_text is None:
        uniprot_text, uniprot_meta = fetch_uniprot(isoforms=isoforms)
    contaminants = tuple(contaminants)
    libs: dict[str, tuple[str, dict]] = {}
    for lib in contaminants:
        if contaminant_texts is not None:
            libs[lib] = (contaminant_texts[lib], {"file": CONTAMINANT_FILES[lib]})
        else:
            libs[lib] = fetch_contaminants(lib, contaminant_commit)

    # Archive every input exactly as used, so the snapshot can be rebuilt without the web.
    src = out / "sources"
    src.mkdir(exist_ok=True)
    (src / "uniprot_retrieval.json").write_text(json.dumps(uniprot_meta or {}, indent=2) + "\n")
    for lib, (text, _) in libs.items():
        (src / f"frankenfield2022_{lib}_{contaminant_commit[:7]}.fasta").write_text(text)
    (src / Path(irt_sheet).name).write_bytes(Path(irt_sheet).read_bytes())

    human = parse_fasta(uniprot_text)
    if not human:
        raise ValueError("UniProt returned no records")
    human_seqs = {s for _, s in human}
    human_by_seq: dict[str, str] = {}
    for h, sq in human:
        human_by_seq.setdefault(sq, accession(h))
    cont, kept_by_seq, filtered, per_lib = [], {}, [], {}
    dropped_human, dropped_dup = [], []
    for lib, (text, meta) in libs.items():
        recs = parse_fasta(text)
        kept = 0
        for h, sq in recs:
            row = {"library": lib, "accession": accession(h), "header": h, "length": len(sq)}
            if HUMAN_TAXON in h and sq not in human_by_seq:
                # A human protein is never labelled a contaminant: it could be endogenous or
                # dust, and the search cannot tell. Its peptides stay with the human proteome;
                # the entry goes on the watchlist for optional post-filtering.
                filtered.append(
                    {**row, "reason": "human protein (watchlist)", "kept_as": None, "sequence": sq}
                )
            elif sq in human_by_seq:
                dropped_human.append(row["accession"])
                filtered.append(
                    {**row, "reason": "identical to human entry", "kept_as": human_by_seq[sq]}
                )
            elif sq in kept_by_seq:
                dropped_dup.append(row["accession"])
                filtered.append(
                    {**row, "reason": "duplicate across libraries", "kept_as": kept_by_seq[sq]}
                )
            else:
                cont.append((h, sq))
                kept_by_seq[sq] = row["accession"]
                kept += 1
        per_lib[lib] = {**meta, "entries_in_library": len(recs), "entries_kept": kept}

    trimmed = []
    if protect:
        cont, trimmed = trim_shared(
            cont, [(h, sq) for h, sq in human if accession(h).split("-")[0] in protect]
        )
        for t in trimmed:
            t["protected"] = ";".join(sorted({protect.get(a, "") for a in t["human_accessions"]}))

    peps, transitions = read_irt_kit(irt_sheet)
    irt = [irt_fusion_record(peps["sequence"].tolist())]
    clash = [p for p in peps["sequence"] if any(p in s for s in human_seqs)]

    parts = {"human.fasta": human, "contaminants.fasta": cont, "irt.fasta": irt}
    for name, recs in parts.items():
        write_fasta(recs, out / name)
    write_fasta(human + cont + irt, out / "search.fasta")
    peps.to_csv(out / "irt_kit_peptides.tsv", sep="\t", index=False)
    # Human-origin entries that differ from every human entry (older or fragment versions):
    # map each to the human accession sharing the most tryptic peptides.
    # (dict: the same entry often appears in both libraries)
    watch_seqs = list({f["accession"]: f["sequence"] for f in filtered if "sequence" in f}.items())
    if watch_seqs:
        # count shared peptides per human accession and keep the best match (ties: lowest
        # accession), so a keratin maps to itself, not to the first paralog alphabetically
        def _digest(seq):
            return {x.replace("I", "L") for x in tryptic_peptides(seq) if 7 <= len(x) <= 30}

        by_pep: dict[str, set[str]] = {}
        for h, sq in human:
            for x in _digest(sq):
                by_pep.setdefault(x, set()).add(accession(h))
        best = {}
        for acc, sq in watch_seqs:
            own = _digest(sq)
            counts: dict[str, int] = {}
            for x in own:
                for a in by_pep.get(x, ()):
                    counts[a] = counts.get(a, 0) + 1
            if counts:
                top = min(counts, key=lambda a: (-counts[a], a))
                best[acc] = (top, f"{counts[top]}/{len(own)}")
        for f in filtered:
            if f.pop("sequence", None) is not None and f["accession"] in best:
                f["kept_as"], f["shared_peptides"] = best[f["accession"]]
    # Every contaminant entry left out, why, and which entry now carries its peptides.
    filt = pd.DataFrame(
        filtered,
        columns=[
            "library",
            "accession",
            "header",
            "length",
            "reason",
            "kept_as",
            "shared_peptides",
        ],
    )
    # Watchlist: human proteins the contaminant libraries list (keratins, KRTAPs, filaggrin,
    # ...), keyed on the current human accession. For optional post-filtering only; several
    # (e.g. KRT7/8/17/18/19) are expressed by epithelial cells such as HeLa.
    gene_of = {}
    for h, _ in human:
        m = re.search(r"GN=(\S+)", h)
        gene_of.setdefault(accession(h), m.group(1) if m else "")
    is_human = filt["header"].str.contains(HUMAN_TAXON, regex=False)
    w = filt[is_human & filt["reason"].ne("duplicate across libraries")].copy()
    w["match"] = w.apply(
        lambda r: (
            "identical sequence"
            if r["reason"] == "identical to human entry"
            else f"shared peptides {r['shared_peptides']}"
        ),
        axis=1,
    )
    watch = (
        w.groupby("kept_as", dropna=False)
        .agg(
            contaminant_entries=("accession", lambda x: ";".join(sorted(set(x)))),
            libraries=("library", lambda x: ";".join(sorted(set(x)))),
            match=("match", lambda x: "; ".join(sorted(set(x)))),
        )
        .reset_index()
        .rename(columns={"kept_as": "human_accession"})
    )
    watch.insert(1, "human_gene", watch["human_accession"].map(gene_of))
    watch = watch.sort_values(["human_gene", "human_accession"])
    watch.to_csv(out / "human_contaminant_watchlist.tsv", sep="\t", index=False)
    filt.to_csv(out / "contaminants_filtered.tsv", sep="\t", index=False)
    shared = shared_peptides(cont, human)
    shared.to_csv(out / "contaminants_shared_with_human.tsv", sep="\t", index=False)
    transitions.to_csv(out / "irt_kit_transitions.tsv", sep="\t", index=False)

    release = (uniprot_meta or {}).get("release") or "unknown"
    manifest = {
        "name": database_name(release, isoforms, contaminants),
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "search_fasta": {
            "file": "search.fasta",
            "sha256": sha256(out / "search.fasta"),
            "entries": len(human) + len(cont) + len(irt),
        },
        "human": {
            **(uniprot_meta or {}),
            "file": "human.fasta",
            "entries": len(human),
            "sha256": sha256(out / "human.fasta"),
        },
        "contaminants": {
            "source": CONTAMINANT_CITATION,
            "repo": CONTAMINANT_REPO,
            "commit": contaminant_commit,
            "libraries": per_lib,
            "file": "contaminants.fasta",
            "entries": len(cont),
            "dropped_identical_to_human": dropped_human,
            "dropped_duplicate_across_libraries": dropped_dup,
            "filtered_table": "contaminants_filtered.tsv",
            "shared_with_human_table": "contaminants_shared_with_human.tsv",
            "human_watchlist_table": "human_contaminant_watchlist.tsv",
            "trimmed_for_protected_human": trimmed,
            "protected_human_accessions": len(protect or {}),
            "protected_source": protect_source,
            "human_watchlist_proteins": len(watch),
            "policy": "no human-origin entry is labelled a contaminant; human proteins the "
            "libraries list are kept as human and recorded on the watchlist",
            "entries_sharing_peptides_with_human": len(shared),
            "sha256": sha256(out / "contaminants.fasta"),
        },
        "irt": {
            "source": "Biognosys iRT kit reference sheet",
            "sheet": f"external://{IRT_KIT_KEY}",
            "sheet_sha256": sha256(irt_sheet),
            "accession": IRT_ACCESSION,
            "peptides": peps["sequence"].tolist(),
            "peptides_also_in_human": clash,
            "file": "irt.fasta",
            "sha256": sha256(out / "irt.fasta"),
        },
        "sources": {f"sources/{f.name}": sha256(f) for f in sorted((out / "sources").iterdir())},
        "spectronaut_parsing_rule": "UniProt FASTA (contaminants are sp|Cont_<acc>|..., "
        f"iRT is sp|{IRT_ACCESSION}|...)",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def publish(out_dir: str | os.PathLike, *, irt_sheet: str | os.PathLike | None = None) -> list[str]:
    """Upload a built database to ``external://search-fasta/<name>/`` (append-only bucket).

    Refuses if that prefix already holds files: a published snapshot is never overwritten. The
    iRT sheet is archived at ``IRT_KIT_KEY`` the first time.
    """
    from deliverome_analysis.s3 import external

    out = Path(out_dir)
    manifest = json.loads((out / "manifest.json").read_text())
    prefix = f"search-fasta/{manifest['name']}/"
    if external.list_keys(prefix=prefix):
        raise FileExistsError(f"external://{prefix} already exists; snapshots are immutable")
    keys = []
    if irt_sheet is not None and not external.exists(IRT_KIT_KEY):
        if sha256(irt_sheet) != manifest["irt"]["sheet_sha256"]:
            raise ValueError("iRT sheet differs from the one this database was built from")
        external.upload(irt_sheet, IRT_KIT_KEY)
        keys.append(IRT_KIT_KEY)
    for f in sorted(out.rglob("*")):
        if f.is_file():
            rel = f.relative_to(out).as_posix()
            external.upload(f, prefix + rel)
            keys.append(prefix + rel)
    return keys


def download(name: str, dest: str | os.PathLike | None = None) -> Path:
    """Fetch a published search.fasta, checking it against its manifest's SHA-256."""
    from deliverome_analysis.s3 import external

    prefix = f"search-fasta/{name}/"
    manifest = json.loads(Path(external.download(prefix + "manifest.json")).read_text())
    path = Path(external.download(prefix + "search.fasta"))
    if sha256(path) != manifest["search_fasta"]["sha256"]:
        raise ValueError(f"search.fasta for {name} does not match its manifest")
    if dest is not None:
        Path(dest).write_bytes(path.read_bytes())
        return Path(dest)
    return path


def write_review_xlsx(
    review: pd.DataFrame, path: str | os.PathLike, db_dir: str | os.PathLike, dataset: str = ""
) -> Path:
    """The review table as a formatted workbook: README, Review, Duplicates, Summary."""
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    d = Path(db_dir)
    m = json.loads((d / "manifest.json").read_text())
    filt = pd.read_csv(d / "contaminants_filtered.tsv", sep="\t")
    dups = filt[filt["reason"] == "duplicate across libraries"][
        ["library", "accession", "header", "kept_as"]
    ]
    summary = (
        review.groupby(["concern", "overlap_degree"], observed=True)
        .size()
        .rename("entries")
        .reset_index()
    )
    readme = [
        (
            "What this is",
            "Every protein-contaminant entry considered for the DIA search FASTA, "
            "ranked by how much it could interfere with endogenous human signal.",
        ),
        ("Search database", m["name"]),
        (
            "UniProt",
            f"{m['human'].get('release')} ({m['human'].get('release_date')}), "
            f"{m['human']['entries']:,} human entries, isoforms={m['human'].get('isoforms')}",
        ),
        (
            "Contaminant libraries",
            f"{m['contaminants']['source']}; commit "
            f"{m['contaminants']['commit']}; libraries: {', '.join(m['contaminants']['libraries'])}",
        ),
        (
            "In the FASTA",
            f"{m['contaminants']['entries']} contaminant entries are in search.fasta "
            "(in_search_fasta = Yes; also listed on their own sheet). Entries marked No were "
            "removed because their sequence is identical to a human entry, which now carries "
            "those peptides; cross-library duplicates are on 'Duplicates removed'. The FASTA "
            "used for the Astral benchmark search contained no contaminants.",
        ),
        ("Dataset used for 'detected'", dataset or "none"),
        (
            "Overlap",
            "Fully tryptic peptides 7-30 aa, I and L treated as equal. "
            "overlap_fraction = shared / total peptides of the contaminant; 'identical' = whole "
            "sequence equals a human entry, which was removed so the human accession owns it.",
        ),
        (
            "Concern: High",
            "Overlaps a canonical-tier surface protein (accessible-surfaceome catalog).",
        ),
        (
            "Concern: Medium",
            "Overlaps a likely-tier surface protein, or a human protein detected in the dataset.",
        ),
        ("Concern: Low", "Overlaps only human proteins not detected in the dataset."),
        ("Concern: No overlap", "No shared tryptic peptide with any human entry."),
        (
            "Human proteins",
            "No human protein is labelled a contaminant: a human protein in the sample could be "
            "endogenous or dust, and the search cannot tell. Human entries in the contaminant "
            "libraries are left out of the contaminant FASTA (their peptides stay with the human "
            "proteome) and listed on 'Human watchlist' for optional post-filtering. Several, e.g. "
            "KRT7/8/17/18/19, are expressed by epithelial cells such as HeLa.",
        ),
        (
            "Regenerate",
            "python -m deliverome_analysis.search_fasta review <db_dir> "
            "--matrix-dir <trimmed export> --out <file.xlsx>",
        ),
    ]
    fills = {"High": "F4C7C3", "Medium": "FCE8B2", "Low": "E2EFDA", "No overlap": "F2F2F2"}
    out = Path(path)
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        pd.DataFrame(readme, columns=["Item", "Detail"]).to_excel(
            xw, sheet_name="README", index=False
        )
        review.to_excel(xw, sheet_name="Review", index=False)
        review[review["in_search_fasta"] == "Yes"].to_excel(
            xw, sheet_name="In search FASTA", index=False
        )
        summary.to_excel(xw, sheet_name="Summary", index=False)
        if (d / "human_contaminant_watchlist.tsv").exists():
            pd.read_csv(d / "human_contaminant_watchlist.tsv", sep="\t").to_excel(
                xw, sheet_name="Human watchlist", index=False
            )
        dups.to_excel(xw, sheet_name="Duplicates removed", index=False)
        for ws in xw.book.worksheets:
            ws.freeze_panes = "A2"
            for c in ws[1]:
                c.font = Font(bold=True)
            ws.auto_filter.ref = ws.dimensions
            for i, col in enumerate(ws.iter_cols(min_row=1, max_row=min(ws.max_row, 200)), 1):
                width = max(len(str(c.value or "")) for c in col)
                ws.column_dimensions[get_column_letter(i)].width = min(max(10, width + 2), 60)
        for name in ("Review", "In search FASTA"):
            rv = xw.book[name]
            pct = next(c.column_letter for c in rv[1] if c.value == "overlap_fraction")
            for row in rv.iter_rows(min_row=2):
                row[0].fill = PatternFill("solid", fgColor=fills.get(str(row[0].value), "FFFFFF"))
                if row[1].value == "No":
                    row[1].font = Font(color="9C0006", bold=True)
            for c in rv[pct][1:]:
                c.number_format = "0%"
        xw.book["README"].column_dimensions["B"].width = 110
        for row in xw.book["README"].iter_rows(min_row=2):
            row[1].alignment = Alignment(wrap_text=True, vertical="top")
    return out


def detected_in_trimmed(matrix_dir: str | os.PathLike) -> dict[str, float]:
    """Base accession -> median precursor log2 intensity in a trim_spectronaut_report folder."""
    d = Path(matrix_dir)
    qc = pd.read_parquet(d / "precursor_qc.parquet", columns=["precursor", "log2_intensity"])
    ann = pd.read_parquet(d / "precursors.parquet").set_index("precursor")
    qc["acc"] = qc["precursor"].map(ann["PG.ProteinAccessions"]).str.split(";")
    qc = qc.explode("acc")
    qc["acc"] = qc["acc"].str.split("-").str[0]
    return qc.groupby("acc")["log2_intensity"].median().to_dict()


# ----------------------------------------------------------------------------------------------
# contaminant review table
# ----------------------------------------------------------------------------------------------

CONCERN_ORDER = ["High", "Medium", "Low", "No overlap"]


def _overlap_degree(frac: float) -> str:
    if frac >= 1.0:
        return "identical"
    if frac >= 0.5:
        return "high (>=50% of peptides)"
    if frac >= 0.1:
        return "partial (10-50%)"
    if frac > 0:
        return "low (<10%)"
    return "none"


def contaminant_review(
    db_dir: str | os.PathLike,
    *,
    surface_tier: dict[str, str] | None = None,
    detected: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Every contaminant entry with its overlap with human proteins and a concern level.

    ``surface_tier`` maps human UniProt accession -> surfaceome tier (canonical / likely);
    ``detected`` maps accession -> median log2 intensity in a real dataset. Concern:

    * High: overlaps (shares a tryptic peptide with, or is identical to) a canonical-tier
      surface protein;
    * Medium: overlaps a likely-tier surface protein, or a human protein detected in the
      dataset (endogenous signal at stake);
    * Low: overlaps only human proteins not detected there;
    * No overlap: no shared peptide.

    Duplicates across contaminant libraries are left out (they change nothing).
    """
    d = Path(db_dir)
    human = parse_fasta((d / "human.fasta").read_text())
    gene = {}
    for h, _ in human:
        m = re.search(r"GN=(\S+)", h)
        gene.setdefault(accession(h).split("-")[0], m.group(1) if m else "")
    kept = parse_fasta((d / "contaminants.fasta").read_text())
    shared = pd.read_csv(d / "contaminants_shared_with_human.tsv", sep="\t").set_index(
        "contaminant"
    )
    filt = pd.read_csv(d / "contaminants_filtered.tsv", sep="\t")
    surface_tier, detected = surface_tier or {}, detected or {}

    def info(header: str) -> dict:
        name = re.search(r"^\S+\s+(.+?)\s+OS=", header)
        org = re.search(r"OS=(.+?)\s+OX=", header)
        gn = re.search(r"GN=(\S+)", header)
        return {
            "protein": name.group(1) if name else "",
            "organism": org.group(1) if org else "",
            "contaminant_gene": gn.group(1) if gn else "",
        }

    rows = []
    for h, _sq in kept:
        acc = accession(h)
        if acc in shared.index:
            r = shared.loc[acc]
            hum = str(r["human_accessions"]).split(";")
            frac = r["shared_peptides"] / max(r["total_peptides"], 1)
            sp, tp = int(r["shared_peptides"]), int(r["total_peptides"])
        else:
            hum, frac, sp, tp = [], 0.0, 0, None
        rows.append(
            {
                "accession": acc,
                **info(h),
                "status": "in search FASTA",
                "shared_peptides": sp,
                "total_peptides": tp,
                "overlap_fraction": frac,
                "human_accessions": hum,
            }
        )
    removed = filt[filt["reason"].isin(["identical to human entry", "human protein (watchlist)"])]
    for r in removed.drop_duplicates("accession").itertuples():
        human_origin = HUMAN_TAXON in r.header
        kept_as = r.kept_as if isinstance(r.kept_as, str) else None
        rows.append(
            {
                "accession": r.accession,
                **info(r.header),
                "status": (
                    f"removed: human protein, on watchlist as {kept_as}"
                    if human_origin
                    else f"removed: identical to human {kept_as}"
                ),
                "shared_peptides": None,
                "total_peptides": None,
                "overlap_fraction": 1.0 if r.reason == "identical to human entry" else None,
                "human_accessions": [kept_as.split("-")[0]] if kept_as else [],
            }
        )
    df = pd.DataFrame(rows)
    df["overlap_degree"] = df["overlap_fraction"].map(_overlap_degree)
    df["human_genes"] = df["human_accessions"].map(lambda a: ", ".join(gene.get(x, x) for x in a))
    df["surface_overlap"] = df["human_accessions"].map(
        lambda a: ", ".join(f"{gene.get(x, x)} ({surface_tier[x]})" for x in a if x in surface_tier)
    )
    df["detected_overlap"] = df["human_accessions"].map(
        lambda a: ", ".join(f"{gene.get(x, x)} ({detected[x]:.1f})" for x in a if x in detected)
    )

    def concern(r) -> str:
        if not r["human_accessions"]:
            return "No overlap"
        if "(canonical)" in r["surface_overlap"]:
            return "High"
        if r["surface_overlap"] or r["detected_overlap"]:
            return "Medium"
        return "Low"

    df["concern"] = df.apply(concern, axis=1)
    df["n_human_proteins"] = df["human_accessions"].map(len)
    df["human_accessions"] = df["human_accessions"].map(";".join)
    df["concern"] = pd.Categorical(df["concern"], CONCERN_ORDER, ordered=True)
    df = df.sort_values(
        ["concern", "overlap_fraction", "shared_peptides"],
        ascending=[True, False, False],
        na_position="first",
    )
    df["in_search_fasta"] = df["status"].eq("in search FASTA").map({True: "Yes", False: "No"})
    man = json.loads((d / "manifest.json").read_text())
    for t in man["contaminants"].get("trimmed_for_protected_human", []):
        hit = df["accession"] == t["accession"]
        df.loc[hit, "status"] = (
            f"in search FASTA, trimmed: {len(t['removed_peptides'])} peptides shared with "
            f"{t.get('protected', 'protected')} protein(s) removed"
        )
    cols = [
        "concern",
        "in_search_fasta",
        "accession",
        "protein",
        "organism",
        "contaminant_gene",
        "status",
        "overlap_degree",
        "overlap_fraction",
        "shared_peptides",
        "total_peptides",
        "n_human_proteins",
        "human_genes",
        "surface_overlap",
        "detected_overlap",
        "human_accessions",
    ]
    return df[cols].reset_index(drop=True)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--irt-sheet", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--canonical-only", action="store_true", help="no isoforms")
    b.add_argument(
        "--protect",
        default="surface-canonical",
        choices=["surface-canonical", "surface", "none"],
        help="human proteins whose peptides contaminants may not claim",
    )
    b.add_argument(
        "--contaminants",
        default=",".join(DEFAULT_CONTAMINANTS),
        help=f"comma-separated, from {sorted(CONTAMINANT_FILES)}",
    )
    p = sub.add_parser("publish")
    p.add_argument("out")
    p.add_argument("--irt-sheet", default=None, help="archive the kit sheet if not yet on S3")
    r = sub.add_parser("review", help="contaminant review workbook")
    r.add_argument("db_dir")
    r.add_argument("--matrix-dir", default=None, help="trimmed export for 'detected'")
    r.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "build":
        protect = protect_source = None
        if a.protect != "none":
            from deliverome_analysis import surfaceome

            tiers = ("canonical",) if a.protect == "surface-canonical" else ("canonical", "likely")
            cat = surfaceome.membership(tiers)
            protect = {
                acc: f"{t} surface"
                for acc, t in zip(cat["uniprot_accession"], cat["tier"], strict=True)
            }
            protect_source = f"external://{surfaceome.CATALOG_KEY} tiers={','.join(tiers)}"
        m = build_search_fasta(
            a.out,
            a.irt_sheet,
            isoforms=not a.canonical_only,
            contaminants=tuple(a.contaminants.split(",")),
            protect=protect,
            protect_source=protect_source,
        )
        print(json.dumps({k: m[k] for k in ("name", "search_fasta")}, indent=2))
        print(
            f"human {m['human']['entries']:,} | contaminants {m['contaminants']['entries']} "
            f"(dropped {len(m['contaminants']['dropped_identical_to_human'])} identical to human, "
            f"{len(m['contaminants']['dropped_duplicate_across_libraries'])} duplicates) "
            f"| iRT {len(m['irt']['peptides'])} peptides"
        )
    elif a.cmd == "review":
        from deliverome_analysis import surfaceome

        cat = surfaceome.membership(("canonical", "likely"))
        det = detected_in_trimmed(a.matrix_dir) if a.matrix_dir else None
        rv = contaminant_review(
            a.db_dir,
            surface_tier=dict(zip(cat["uniprot_accession"], cat["tier"], strict=True)),
            detected=det,
        )
        dataset = ""
        if a.matrix_dir and (Path(a.matrix_dir) / "manifest.json").exists():
            dataset = (
                Path(
                    json.loads((Path(a.matrix_dir) / "manifest.json").read_text())["source_report"]
                ).name.split(".tsv")[0]
                + ".tsv"
            )
        print(write_review_xlsx(rv, a.out, a.db_dir, dataset))
        print(rv["concern"].value_counts().reindex(CONCERN_ORDER).to_string())
    else:
        for k in publish(a.out, irt_sheet=a.irt_sheet):
            print(f"external://{k}")
