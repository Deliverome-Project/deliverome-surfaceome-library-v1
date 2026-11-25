"""Put the tag cassette behind the signal peptide of every ready L1 "N after SP" member.

This is the N-after-SP rung of issue #2 ("assemble the tag cassette into every L1 insert").
Every input is fixed in the record; nothing is chosen here:

* **which members, how and where**: ``reviewer-package/nterm/nterm_decisions.tsv``, one
  decided arm per member (``nterm_decide.py``). ``native_site`` puts the tandem after native
  residue ``n``. ``igk_swap`` replaces residues ``1..c`` with the Igκ leader plus Asp
  (``build_gpcr_rescue_arms.LEADER`` / ``PLUS1``, as ``build_sp_swap_arms.py`` does for
  VSIG10) and puts the tandem after the same mature residue. ``no_nterm_arm`` members are
  not built. Each row carries the decision's ``decided_by``, ``reason`` and rejected
  alternative: Becca's rule for sites before the fold (2026-11-11), and this session's
  delegated calls for the rest (Becca, 2026-11-12, "use your best judgment").
* **the ORF**: the member's domesticated ORF from ``candidate_constructs.json`` (stage 10),
  which must translate to the member's designed protein plus a stop.
* **the cassette**: the 51-aa tandem (Becca, 2026-11-10: "we want the tandem"), on this arm's own
  linker encoding from ``l1_arm_cassette_encodings.tsv`` (``diversify_arm_cassettes.py``):
  SpyTag003 and P-ALFA-P on p003's codons, the four GGGGS blocks recoded per arm. Until
  2026-11-21 every N-terminal arm carried p003's one encoding, on the reasoning that an
  upstream cassette never enters the tract to the barcode; Becca, 2026-11-21: "we still don't
  want shared linkers in the middle", since a template switch inside a shared cassette joins two
  arms' ORFs.
  The full cassette, leading ``GGGGS`` included, is kept everywhere (Becca, 2026-11-11: "keep
  full"; relayed by the l1-recode-dups session). With SpyTag003 directly after Met1, DeepTMHMM v1
  called a spurious signal peptide inside the tag on 12 of 62 type III constructs and on none with
  the full cassette (PR #282). Do not trim it.
* **the layout**: the C-terminal builder's (``build_l1_cterm_arms.py``, PR #284), in the Esp3I
  wrapper ``assemble_candidates.py`` puts every L1 insert in::

      CGTCTCA | Kozak | ORF[:3n] | cassette | ORF[3n:] (native stop) | CC
              | T7 | 5' handle | barcode (20 N) | Feldman 3' | 10x CS1 | AGA GAGACG

Per member it asserts:
- the cassette goes in after nucleotide ``3n`` of the ORF;
- the protein is exactly the tagged form SignalP scored: ``native[:n] + cassette +
  native[n:]``, or for an Igκ arm ``leader + D + native[c:n] + cassette + native[n:]``. It must
  hash to the decision's ``tagged_sha256``;
- the coding sequence has one stop, at its end;
- the insert has the two wrapper Esp3I sites and no other BbsI/Esp3I/BsaI site, alone and in
  the p0039 backbone.

A seam that completes a Type IIS site is reported and the member held, never silently recoded.

**Attribution** (convention of 2026-11-10). Decided by Becca: the cassette on 2026-11-10, the
site rule on 2026-11-11, and on 2026-11-11 that this session build the arms. Implemented by the
Claude Code session ``nterm-reviewer-package`` (PR #254). The C-terminal layout and checks are
the l1-cterm-build session's (PR #284), restated here because that PR is not merged.

    uv run --locked python build_l1_nterm_arms.py            # build and write
    uv run --locked python build_l1_nterm_arms.py --check    # rebuild, compare with the record
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import tempfile
from collections import Counter
from pathlib import Path

from build_constructs import CASSETTE_AA
from diversify_arm_cassettes import arm_cassette

from deliverome_analysis import library_design as ld
from deliverome_analysis import plasmid_qc as qc
from deliverome_analysis import seqtools

LEADER = qc.IGK_LEADER_NT
LEADER_AA = qc.IGK_LEADER_AA
PLUS1 = qc.IGK_PLUS1_NT
PLUS1_AA = qc.IGK_PLUS1_AA

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "library/l1/proteoform-selection"))

MANIFEST = HERE / "candidate_manifest_tranches.tsv"
DECISIONS = HERE / "reviewer-package/nterm/nterm_decisions.tsv"
OUT_TSV = HERE / "constructs_l1_nterm.tsv"
OUT_FASTA = HERE / "constructs_l1_nterm.fasta"
OUT_JSON = HERE / "constructs_l1_nterm_summary.json"
#: The domesticated-ORF record versioned with the designed CDS (library_design.CONSTRUCTS_KEY).
CONSTRUCTS_KEY = ld.CONSTRUCTS_KEY

BARCODE_SLOT = qc.BARCODE_SLOT
ENZYMES = list(qc.DOMESTICATION_ENZYMES)
BUILT_ARMS = ("native_site", "igk_swap")
NON_BUILT_ARMS = ("no_nterm_arm", "routed_internal")
VALID_ARMS = BUILT_ARMS + NON_BUILT_ARMS

COLUMNS = [
    "candidate_id",
    "design_id",
    "gene_symbol",
    "candidate_set",
    "arm",
    "native_aa",
    "protein_aa",
    "orf_nt",
    "coding_nt",
    "insert_nt",
    "placement",
    "igk_replaces_residues",
    "cassette_after_residue",
    "native_after_residue",
    "junction",
    "cassette_start_nt",
    "type_iis_sites",
    "inside_predicted_fold",
    "tagged_cleavage_probability",
    "refold_verdict",
    "decided_by",
    "implemented_by",
    "reason",
    "rejected_alternative",
    "seam_edit",
    "protein_sha256",
]


def read_table(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _index_unique(rows: list[dict], key: str, label: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for row in rows:
        k = (row.get(key) or "").strip()
        if not k:
            raise ValueError(f"{label}: missing/empty {key!r}")
        if k in out:
            raise ValueError(f"{label}: duplicate {key} {k!r}")
        out[k] = row
    return out


def _parse_int_strict(raw: object, field: str, ctx: str) -> int:
    text = ("" if raw is None else str(raw)).strip()
    if not text:
        raise ValueError(f"{ctx}: missing/empty {field}")
    try:
        val = float(text)
    except ValueError as exc:
        raise ValueError(f"{ctx}: non-numeric {field} {raw!r}") from exc
    if not math.isfinite(val) or val != int(val):
        raise ValueError(f"{ctx}: invalid {field} {raw!r}")
    return int(val)


def build_one(d: dict, orf: str, protein: str, tag_fixed: dict) -> dict | str:
    """The insert for one decided member, or a reason it cannot be built as decided.

    ``native_site``: ORF[:3n] + cassette + ORF[3n:]. ``igk_swap``: Igκ leader + Asp codon
    replace codons 1..c, then ORF[3c:3n] + cassette + ORF[3n:] (n is the native residue).
    """
    cid = (d.get("candidate_id") or "").strip()
    assert cid, "nterm decision missing candidate_id"
    assert d.get("arm") in BUILT_ARMS, f"{cid}: unexpected arm {d.get('arm')!r}"
    n = _parse_int_strict(d.get("native_after_residue"), "native_after_residue", cid)
    assert 0 < n < len(protein), (
        f"{cid}: native_after_residue {n} out of range 1..{len(protein) - 1}"
    )
    assert seqtools.translate(orf) == protein + "*", f"{cid}: ORF does not encode its protein"
    cas = arm_cassette(cid)
    assert all(cas[3 * i : 3 * i + 3] == c for i, c in tag_fixed.items()), "not p003's tags"
    if d["arm"] == "igk_swap":
        parts = (d.get("igk_replaces_residues") or "").split("-")
        assert len(parts) == 2, (
            f"{cid}: malformed igk_replaces_residues {d.get('igk_replaces_residues')!r}"
        )
        c = _parse_int_strict(parts[1], "igk_replaces_residues", cid)
        assert 0 < c <= n, f"{cid}: igk_replaces_residues end {c} outside 1..{n}"
        head = LEADER + PLUS1 + orf[3 * c : 3 * n]
        expect = LEADER_AA + PLUS1_AA + protein[c:n] + CASSETTE_AA + protein[n:]
        after = len(LEADER_AA) + len(PLUS1_AA) + (n - c)
        assert after == _parse_int_strict(
            d.get("insert_after_residue"), "insert_after_residue", cid
        ), f"{cid}: Igκ site does not map"
    else:
        head, expect, after = orf[: 3 * n], protein[:n] + CASSETTE_AA + protein[n:], n
    coding = head + cas + orf[3 * n :]
    translated = seqtools.translate(coding)
    assert translated == expect + "*", f"{cid}: protein is not the decided tagged form"
    assert translated.count("*") == 1, f"{cid}: internal stop"
    sha = hashlib.sha256(expect.encode()).hexdigest()
    assert sha == d["tagged_sha256"], f"{cid}: not the protein SignalP scored"
    insert = qc.wrap_insert(coding)
    seam_edit = ""
    if not qc.wrapper_sites_only(insert):
        # Objective fix: a seam created a BbsI/Esp3I/BsaI site. Try one synonymous change in the
        # native codons flanking the cassette (2 either side), protein unchanged; the cassette's
        # own codons (p003 tag codons + per-arm diversified linkers) are never touched.
        head, rest, seam_edit = seam_fix(head, orf[3 * n :], cas)
        if not seam_edit:
            return f"Type IIS site at a seam: {seqtools.find_sites(insert, ENZYMES)}"
        coding = head + cas + rest
        assert seqtools.translate(coding) == expect + "*", f"{cid}: seam fix changed the protein"
        insert = qc.wrap_insert(coding)
    if not qc.wrapper_sites_only(insert):
        return f"Type IIS site at a seam: {seqtools.find_sites(insert, ENZYMES)}"
    lo = qc.CODING_START + len(head)
    return {
        "candidate_id": cid,
        "design_id": d["design_id"],
        "gene_symbol": d["gene_symbol"],
        "candidate_set": d["candidate_set"],
        "arm": d["arm"],
        "native_aa": len(protein),
        "protein_aa": len(expect),
        "orf_nt": len(orf),
        "coding_nt": len(coding),
        "insert_nt": len(insert),
        "placement": "N after SP",
        "igk_replaces_residues": d["igk_replaces_residues"],
        "cassette_after_residue": after,
        "native_after_residue": n,
        "junction": f"{protein[n - 1]}{n}|{protein[n]}",
        "cassette_start_nt": lo + 1,
        "type_iis_sites": "Esp3I x2 (wrapper)",
        "inside_predicted_fold": d["inside_predicted_fold"],
        "tagged_cleavage_probability": d["tagged_cleavage_probability"],
        "refold_verdict": d["refold_verdict"],
        "decided_by": d["decided_by"],
        "implemented_by": d["implemented_by"],
        "reason": d["reason"]
        + (
            f"; seam edit {seam_edit} (objective fix: a seam Type IIS site, protein unchanged; "
            "implemented by nterm-reviewer-package, PR #254)"
            if seam_edit
            else ""
        ),
        "rejected_alternative": d["rejected_alternative"],
        "seam_edit": seam_edit,
        "protein_sha256": sha,
        "insert": insert,
    }


def seam_fix(head: str, rest: str, cas: str) -> tuple[str, str, str]:
    """One synonymous codon change beside the cassette that clears every seam site."""

    def clean(h: str, r: str) -> bool:
        return qc.wrapper_sites_only(qc.wrap_insert(h + cas + r))

    syn = seqtools.synonymous_codons()
    for side, k in (("before", 1), ("before", 2), ("after", 0), ("after", 1)):
        seq = head if side == "before" else rest
        i = len(seq) - 3 * k if side == "before" else 3 * k
        if side == "before" and i < 3:  # never touch the start codon
            continue
        if side == "after" and i + 3 > len(seq) - 3:  # never touch the stop codon
            continue
        old = seq[i : i + 3]
        for codon, _ in syn[seqtools.translate(old)]:
            if codon == old:
                continue
            new = seq[:i] + codon + seq[i + 3 :]
            h, r = (new, rest) if side == "before" else (head, new)
            if clean(h, r):
                return h, r, f"{old}->{codon}, codon {k} {side} the cassette"
    return head, rest, ""


def vector_check(rows: list[dict]) -> int:
    """Every arm in the p0039 backbone: backbone sites plus the two wrapper Esp3I, no more."""
    import candidate_inputs
    from Bio import SeqIO

    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    return qc.check_inserts_in_backbone([(r["candidate_id"], r["insert"]) for r in rows], bb)


def build() -> tuple[list[dict], dict]:
    from deliverome_analysis.s3 import processed

    decisions = read_table(DECISIONS)
    _index_unique(decisions, "candidate_id", DECISIONS.name)
    _index_unique(decisions, "design_id", DECISIONS.name)
    bad_arms = [
        (d["candidate_id"], d.get("arm")) for d in decisions if d.get("arm") not in VALID_ARMS
    ]
    assert not bad_arms, f"{DECISIONS.name}: unknown arm values: {bad_arms[:5]}"
    manifest = _index_unique(
        [r for r in read_table(MANIFEST) if ld.is_active(r)], "candidate_id", MANIFEST.name
    )
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    constructs = json.loads(Path(processed.download(CONSTRUCTS_KEY)).read_text())
    tag_fixed = qc.p003_tag_codons()
    rows, held = [], []
    expected_built = 0
    for d in decisions:
        if d["arm"] not in BUILT_ARMS:
            continue
        expected_built += 1
        cid = d["candidate_id"]
        assert cid in manifest, f"{cid}: built N-after-SP decision missing from active manifest"
        assert cid in proteins, f"{cid}: built N-after-SP decision missing designed protein"
        orf = constructs[manifest[cid]["key"]]["construct"]
        out = build_one(d, orf, proteins[cid]["sequence"], tag_fixed)
        if isinstance(out, str):
            held.append({"candidate_id": cid, "gene_symbol": d["gene_symbol"], "reason": out})
        else:
            rows.append(out)
    assert len(rows) + len(held) == expected_built, (
        f"silent row drop: built {len(rows)} + held {len(held)} != expected {expected_built}"
    )
    n_bb = vector_check(rows)
    summary = {
        "tranche": "l1",
        "rung": "N after SP",
        "decisions": dict(sorted(Counter(d["arm"] for d in decisions).items())),
        "built": len(rows),
        "built_by_arm": dict(sorted(Counter(r["arm"] for r in rows).items())),
        "built_inside_predicted_fold": sum(r["inside_predicted_fold"] == "1" for r in rows),
        "seam_failures": held,
        "no_nterm_arm": sorted(
            f"{d['candidate_id']} {d['gene_symbol']}"
            for d in decisions
            if d["arm"] == "no_nterm_arm"
        ),
        "by_set": dict(sorted(Counter(r["candidate_set"] for r in rows).items())),
        "type_iis": f"every insert: 2 Esp3I (wrapper), 0 BbsI, 0 BsaI; in p0039: backbone's "
        f"{n_bb} + 2",
        "insert_nt_range": [min(r["insert_nt"] for r in rows), max(r["insert_nt"] for r in rows)],
        "cassette": "full p003 tandem, this arm's own linker encoding (l1_arm_cassette_encodings.tsv)",
        "barcode": "20-N placeholder; barcodes not reserved (issue #1)",
        "swap": "linkers recoded per arm (Becca, 2026-11-21: no shared linkers in the middle)",
    }
    return rows, summary


def write(rows: list[dict], summary: dict, tsv: Path, fasta: Path, js: Path) -> None:
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)
    fasta.write_text(
        "".join(
            f">{r['candidate_id']}|{r['gene_symbol']}|{r['arm']}|after{r['cassette_after_residue']}"
            f"|{r['insert_nt']}nt\n"
            + "\n".join(r["insert"][i : i + 70] for i in range(0, len(r["insert"]), 70))
            + "\n"
            for r in rows
        )
    )
    js.write_text(json.dumps(summary, indent=2) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="rebuild and compare with the record")
    args = ap.parse_args(argv)
    rows, summary = build()
    if args.check:
        tmp = Path(tempfile.mkdtemp())
        write(rows, summary, tmp / "t.tsv", tmp / "f.fasta", tmp / "s.json")
        for new, rec in (
            (tmp / "t.tsv", OUT_TSV),
            (tmp / "f.fasta", OUT_FASTA),
            (tmp / "s.json", OUT_JSON),
        ):
            assert new.read_bytes() == rec.read_bytes(), f"{rec.name} does not reproduce"
        print(f"all {len(rows)} N-after-SP arms rebuild byte for byte")
        return 0
    write(rows, summary, OUT_TSV, OUT_FASTA, OUT_JSON)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k not in ("seam_failures", "no_nterm_arm")},
            indent=2,
        )
    )
    for h in summary["seam_failures"]:
        print(f"HELD {h['candidate_id']} {h['gene_symbol']}: {h['reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
