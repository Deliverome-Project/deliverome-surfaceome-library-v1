"""Which structure the viewer should show for each L1 member, and which it should not.

The viewer currently shows every member the AlphaFold DB model for its canonical
accession. For 296 of 2,725 that is not the protein we are ordering: 70 had no model at
all, 200 have a model of a different proteoform, and 26 have one that is our protein
embedded in a longer sequence and therefore numbered wrong without an offset.

This emits one row per (member, structure). The row the viewer should render carries
``is_chosen=1``; the AlphaFold DB model it replaces is still emitted with ``is_chosen=0``
and a note saying why, so an alternate stays visible and legible rather than disappearing.

    uv run --locked python build_library_structure_payload.py
"""

from __future__ import annotations

import csv
import io
import re
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld

HERE = Path(__file__).parent
OUT = HERE / "library_structure_payload.tsv"
RELEASE = "colabfold-2026-10-14"
R2_PREFIX = "library-colabfold"
AFDB = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb"

COLS = [
    "release",
    "candidate_id",
    "uniprot_acc",
    "gene_symbol",
    "species",
    "candidate_set",
    "source",
    "r2_key",
    "s3_key",
    "model_url",
    "proteoform_match",
    "residue_offset",
    "ordered_length_aa",
    "model_length_aa",
    "is_chosen",
    "mean_plddt",
    "note",
]


S3_FOLDS = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/"
#: ColabFold runs mirrored to the viewer's R2 bucket under ``library-colabfold/<run>/``. Later
#: runs are on S3 only until they are mirrored, so they get no ``r2_key``.
R2_RUNS = frozenset({"2026-10-14", "2026-10-14-mismatched"})


def fold_locations() -> dict[str, list[str]]:
    """Every published fold: ``{name: [run, ...]}``, newest run first.

    ``name`` is the accession (missing-model folds) or the design id (mismatched-model folds).
    Read from S3 so the payload names the run that actually holds the file, rather than
    assuming the first release's folder (Becca's attribution rule: l1-reviewer-package fixed
    this on 2026-11-12, PR #283, after the membrane-side gate found five keys that did not
    exist).
    """
    from deliverome_analysis.core.s3 import processed

    out: dict[str, list[str]] = {}
    for key in processed.list_keys(prefix=S3_FOLDS):
        if key.endswith(".pdb"):
            run, name = key[len(S3_FOLDS) :].split("/", 1)
            out.setdefault(name[: -len(".pdb")], []).append(run)
    for runs in out.values():
        runs.sort(reverse=True)
    return out


def fold_keys(locations: dict[str, list[str]], name: str, *, stale: bool = False) -> dict:
    """``r2_key`` and ``s3_key`` for a fold: the newest run, or for a stale fold the oldest."""
    runs = locations.get(name, [])
    if not runs:
        return {"r2_key": "", "s3_key": ""}
    run = runs[-1] if stale else runs[0]
    return {
        "r2_key": f"{R2_PREFIX}/{run}/{name}.pdb" if run in R2_RUNS else "",
        "s3_key": f"processed://{S3_FOLDS}{run}/{name}.pdb",
    }


def rows(name: str) -> list[dict]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def plddt_map() -> dict[str, float]:
    """Mean pLDDT per fold, keyed the way that fold was named (accession or design id).

    Read from the ``fold_summary.tsv`` published beside the models in S3 -- their durable
    home, and written by the same run that wrote the coordinates -- together with the
    committed summaries.

    The mismatched-run scores used to be read from
    ``~/.cache/deliverome-analysis/predicted_structures_mismatched/``, which made every
    score depend on whether the machine building this table happened to be the machine that
    ran the fold. It usually was not: on 2026-11-15 all 182 chosen models of the
    ``2026-10-14-mismatched`` run carried a score in S3 and an empty ``mean_plddt`` here,
    so the public gene page offered 182 genes a chosen structure with no confidence to
    judge it by. S3 is authoritative where a name appears in both.
    """
    from deliverome_analysis.core.s3 import processed

    out: dict[str, float] = {}

    def absorb(reader) -> None:
        for r in reader:
            if r.get("mean_plddt"):
                out[r["name"]] = float(r["mean_plddt"])

    for path in (
        HERE / "predicted_structures/fold_summary.tsv",
        HERE / "predicted_structures/fold_summary_mismatched.tsv",
    ):
        if path.exists():
            with path.open() as fh:
                absorb(csv.DictReader(fh, delimiter="\t"))
    for key in sorted(processed.list_keys(prefix=S3_FOLDS)):
        if key.endswith("/fold_summary.tsv"):
            with processed.open_raw(key) as fh:
                absorb(csv.DictReader(io.StringIO(fh.read().decode()), delimiter="\t"))
    return out


def folded_sequences() -> dict[str, str]:
    """The protein each mismatched-model fold was run on, keyed by design id."""
    out: dict[str, str] = {}
    name = None
    for line in (HERE / "fold_inputs_mismatched.fasta").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[0]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


def ordered_proteins() -> dict[str, str]:
    prot = ld.designed_proteins(HERE / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    return {cid: m["sequence"] for cid, m in prot.items()}


def ordered_lengths() -> dict[str, int]:
    """Residue count of the protein each member encodes, from the designed CDS itself.

    This used to be computed from the manifest as ``orf_nt // 3 - 1``, which assumes every
    designed ORF ends in a terminator codon. Thirteen of the 7,016 do not -- three of them
    in L1 (KIR2DL5A, KIR2DS2, KIR2DS5) -- so for those the formula returned a length one
    residue short of the protein. That is not a cosmetic error: this table is paired with
    ``library_topology_index.tsv``, whose length comes from the translated sequence, and a
    consumer that trusts both reads the per-residue topology one residue out of register
    against the model. Taking the length from the sequence removes the assumption and makes
    the two tables agree by construction rather than by coincidence.

    Members with no designed CDS get no length. They are the rescue arms whose sequence is
    still an open question, and a computed-looking number for them is worse than a blank.
    """
    prot = ld.designed_proteins(HERE / "candidate_manifest_tranches.tsv", library=ld.LIBRARY_L1)
    return {cid: len(m["sequence"]) for cid, m in prot.items()}


#: Untagged native folds of whole designed proteins, published by the tag-refold work.
S3_TAGGED = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/"
#: The payload release a dated native-fold batch is reported under (the labels its publisher
#: used); any other batch is reported as ``tagged-<batch>``.
BATCH_RELEASE = {
    "2026-11-21-coreceptor": "coreceptor-additions-2026-11-21",
    "2026-11-22-natives": "natives-2026-11-22",
}
#: Release for an AlphaFold DB model adopted by sequence rather than through the gap table.
ADOPTED_AFDB_RELEASE = "natives-2026-11-22"


def native_models() -> dict[str, dict]:
    """Every model of a whole native protein, keyed by the SHA-256 of the protein it models.

    Two committed indexes, both matched on the sequence rather than an accession:
    ``tag_refold_native_models.tsv`` rows with ``match == identical`` (AlphaFold DB v6
    models and ColabFold folds whose own residues were compared), then the untagged
    ``kind == native`` folds in ``tag_refold_structures.tsv``. An AlphaFold DB model of the
    exact protein wins over a fold of it. ``ours_by_construction`` index rows are skipped:
    they were inferred from an earlier payload by length, not compared residue by residue.
    """
    out: dict[str, dict] = {}
    for r in rows("tag_refold_native_models.tsv"):
        rank = 0 if r["source"] == "alphafold_db" else 1
        h = r["native_sha256"]
        if r["match"] == "identical" and (h not in out or rank < out[h]["rank"]):
            out[h] = {"rank": rank, "source": r["source"], "ref": r["ref"]}
    for r in rows("tag_refold_structures.tsv"):
        if r["kind"] == "native" and r["native_sha256"] not in out:
            out[r["native_sha256"]] = {"rank": 1, "source": "colabfold", "ref": r["pdb_key"]}
    return out


def _fold_plddt(ref: str) -> str:
    """Mean CA pLDDT of a published fold, read from its PDB (B-factor column)."""
    from deliverome_analysis.core.s3 import processed

    path = Path(processed.download(ref.removeprefix("processed://")))
    b = [
        float(ln[60:66])
        for ln in path.read_text().splitlines()
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA"
    ]
    return f"{sum(b) / len(b):.2f}"


def adopt_native_models(out: list[dict], gap: dict[str, dict], ordered: dict[str, str]) -> int:
    """Choose a model, by sequence, for every member the rules above left without one.

    The rules above know the 2026-10-14 AlphaFold DB comparison and the 2026-10/11 ColabFold
    runs. Members added or re-referenced since (the 2026-11-21 co-receptor additions and the
    2026-11-22 rebuild: gap genes, allele restorations, isoform refills, ortholog
    restorations, GTF MANE Select switches) were folded untagged by the tag-refold work
    instead, and are found here in :func:`native_models` by the SHA-256 of the designed
    protein. A member with no model of its exact protein anywhere keeps no chosen row: it is
    reported, never given a near miss. Rows are reordered so each member's chosen row leads.
    """
    import hashlib

    natives = native_models()
    chosen = {r["candidate_id"] for r in out if r["is_chosen"] == 1}
    order = {cid: i for i, cid in enumerate(gap)}
    adopted, unmodelled = 0, []
    for cid, g in gap.items():
        if cid in chosen or cid not in ordered:
            continue
        seq = ordered[cid]
        m = natives.get(hashlib.sha256(seq.encode()).hexdigest())
        if m is None:
            unmodelled.append(cid)
            continue
        row = {
            "candidate_id": cid,
            "uniprot_acc": g["accession"],
            "gene_symbol": g["symbol"],
            "species": g["species"],
            "candidate_set": g["set"],
            "r2_key": "",
            "proteoform_match": "identical",
            "residue_offset": 0,
            "ordered_length_aa": len(seq),
            "model_length_aa": len(seq),
            "is_chosen": 1,
        }
        if m["source"] == "alphafold_db":
            row |= {
                "release": ADOPTED_AFDB_RELEASE,
                "uniprot_acc": m["ref"].rsplit("AF-", 1)[1].split("-F1-")[0],
                "source": "alphafold_db",
                "s3_key": "",
                "model_url": m["ref"],
                "mean_plddt": "",
                "note": "AlphaFold DB v6 model identical to the designed protein (by sequence)",
            }
        else:
            batch = m["ref"].removeprefix(f"processed://{S3_TAGGED}").split("/")[0]
            row |= {
                "release": BATCH_RELEASE.get(batch, f"tagged-{batch}"),
                "source": "colabfold",
                "s3_key": m["ref"],
                "model_url": "",
                "proteoform_match": "ours_by_construction",
                "mean_plddt": _fold_plddt(m["ref"]),
                "note": f"folded from the designed protein (untagged native fold, {batch})",
            }
        out.append(row)
        adopted += 1
    out.sort(key=lambda r: (order.get(r["candidate_id"], len(order)), r["is_chosen"] != 1))
    if unmodelled:
        print(f"  {len(unmodelled)} member(s) have no model of their own protein: {unmodelled}")
    return adopted


#: S3 run folders mirrored to the viewer's R2 bucket after the 2026-10-14 runs (the 2026-11-22-v8
#: viewer load, load l1d-bd6f1c9d3a50f968): their models get the R2 key the viewer renders.
R2_MIRRORED = {
    "colabfold/2026-11-10-l1-orthologs": "library-colabfold",
    "colabfold/2026-11-10-l1-ortholog-mismatched": "library-colabfold",
    "tagged/2026-11-12": "library-colabfold-tagged",
    "tagged/2026-11-21-coreceptor": "library-colabfold-tagged",
    "tagged/2026-11-22-natives": "library-colabfold-tagged",
    "tagged/2026-11-23-natives": "library-colabfold-tagged",
}
_S3_RUN = re.compile(r"predicted-structures/((colabfold|tagged)/[^/]+)/([^/]+)$")


def mirrored_r2_key(s3_key: str) -> str:
    """The R2 key of a model in a mirrored S3 run folder, else ''."""
    m = _S3_RUN.search(s3_key or "")
    if not m or m.group(1) not in R2_MIRRORED:
        return ""
    return f"{R2_MIRRORED[m.group(1)]}/{m.group(1).split('/', 1)[1]}/{m.group(3)}"


def write_payload(out: list[dict], members: set[str]) -> None:
    import collections

    for r in out:
        if not r.get("r2_key"):
            r["r2_key"] = mirrored_r2_key(r.get("s3_key", ""))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(out)
    chosen = [r for r in out if r["is_chosen"] == 1]
    print(f"{len(out)} rows for {len({r['candidate_id'] for r in out})} members")
    print(
        f"  chosen   {len(chosen)}  by source {dict(collections.Counter(r['source'] for r in chosen))}"
    )
    alt = [r for r in out if r["is_chosen"] == 0]
    print(
        f"  alternate {len(alt)}  by match {dict(collections.Counter(r['proteoform_match'] for r in alt))}"
    )
    missing = {r["candidate_id"] for r in out} ^ members
    print(f"  members with no row at all: {len(missing)}")
    print(f"wrote {OUT.name}")


def main() -> int:
    # Active members only: the measurement caches keep rows for retired members (dropped
    # from a panel, or duplicates), which the viewer must not offer as library members.
    active = {
        r["candidate_id"]
        for r in rows("candidate_manifest_tranches.tsv")
        if r["library"] == ld.LIBRARY_L1 and ld.is_active(r)
    }
    gap = {
        r["candidate_id"]: r for r in rows("l1_structure_gap.tsv") if r["candidate_id"] in active
    }
    match = {r["candidate_id"]: r for r in rows("proteoform_match.tsv")}
    # Every successful missing-model fold, from the committed summary: the 2026-11-10 run
    # (61 adopted orthologs) was folded from its own FASTA and is not in fold_inputs.tsv.
    folded_no_model = {r["accession"] for r in rows("fold_inputs.tsv")} | {
        r["name"]
        for r in rows("predicted_structures/fold_summary.tsv")
        if r["ok"] in ("True", "true", "1")
    }
    folded_mismatch = {r["name"] for r in rows("fold_inputs_mismatched.tsv")}
    # A fold is "ours by construction" only while the member still orders the protein it was
    # folded from; a later reference fix (SLC22A17's ATG, TSPAN10's swap, 2026-11-12) leaves
    # the fold describing the old protein, so it is offered but not chosen.
    folded_seq = folded_sequences()
    ordered = ordered_proteins()
    plddt = plddt_map()
    length = ordered_lengths()
    offsets = {r["candidate_id"]: r for r in rows("residue_offsets.tsv")}
    # Fold inputs, their structure files and their pLDDT are keyed on the design id (the
    # fold ran before the member ids were assigned, and member ids move until the freeze).
    design = ld.design_ids(HERE / "candidate_manifest_tranches.tsv")
    locations = fold_locations()

    out: list[dict] = []
    for cid, g in gap.items():
        base = {
            "release": RELEASE,
            "candidate_id": cid,
            "uniprot_acc": g["accession"],
            "gene_symbol": g["symbol"],
            "species": g["species"],
            "candidate_set": g["set"],
            "residue_offset": 0,
            "ordered_length_aa": length.get(cid, ""),
        }
        has_model = g["has_alphafold"] in ("True", "true", "1")
        m = match.get(cid)
        verdict = (m or {}).get("verdict", "unknown" if has_model else "no_model")
        model_len = (m or {}).get("model_aa", "")

        if not has_model and g["accession"] in folded_no_model:
            out.append(
                {
                    **base,
                    "source": "colabfold",
                    **fold_keys(locations, g["accession"]),
                    "model_url": "",
                    "proteoform_match": "ours_by_construction",
                    "model_length_aa": base["ordered_length_aa"],
                    "is_chosen": 1,
                    "mean_plddt": plddt.get(g["accession"], ""),
                    "note": "AlphaFold DB has no model for this proteoform; folded here",
                }
            )
            continue

        did = design.get(cid, cid)
        if (
            verdict == "different"
            and did in folded_mismatch
            and folded_seq.get(did) != ordered.get(cid)
        ):
            out.append(
                {
                    **base,
                    "source": "colabfold",
                    **fold_keys(locations, did, stale=True),
                    "model_url": "",
                    "proteoform_match": "stale_fold",
                    "model_length_aa": len(folded_seq.get(did, "")),
                    "is_chosen": 0,
                    "mean_plddt": plddt.get(did, ""),
                    "note": "folded from an earlier designed protein; refold the current one",
                }
            )
            out.append(
                {
                    **base,
                    "source": "alphafold_db",
                    "r2_key": "",
                    "s3_key": "",
                    "model_url": AFDB.format(acc=g["accession"]),
                    "proteoform_match": "different",
                    "model_length_aa": model_len,
                    "is_chosen": 0,
                    "mean_plddt": "",
                    "note": f"AlphaFold DB models a different proteoform "
                    f"({model_len} aa vs {base['ordered_length_aa']} aa ordered)",
                }
            )
            continue
        if verdict == "different" and did in folded_mismatch:
            out.append(
                {
                    **base,
                    "source": "colabfold",
                    **fold_keys(locations, did),
                    "model_url": "",
                    "proteoform_match": "ours_by_construction",
                    "model_length_aa": base["ordered_length_aa"],
                    "is_chosen": 1,
                    "mean_plddt": plddt.get(did, ""),
                    "note": "folded from the designed CDS",
                }
            )
            out.append(
                {
                    **base,
                    "source": "alphafold_db",
                    "r2_key": "",
                    "s3_key": "",
                    "model_url": AFDB.format(acc=g["accession"]),
                    "proteoform_match": "different",
                    "model_length_aa": model_len,
                    "is_chosen": 0,
                    "mean_plddt": "",
                    "note": f"AlphaFold DB models a different proteoform "
                    f"({model_len} aa vs {base['ordered_length_aa']} aa ordered)",
                }
            )
            continue

        if verdict == "contained":
            # The offset must be MEASURED. ``proteoform_match.tsv`` ships 0 on every
            # contained row because nobody computed it there, and reading that 0 as "no
            # shift needed" chose the wrong span for 25 of the 26. Where
            # compute_residue_offsets.py located the ordered protein inside the model, the
            # model is usable and is chosen with the span to render; where it did not, the
            # member is held back rather than rendered at a guessed offset.
            o = offsets.get(cid)
            if o is None:
                out.append(
                    {
                        **base,
                        "source": "alphafold_db",
                        "r2_key": "",
                        "s3_key": "",
                        "model_url": AFDB.format(acc=g["accession"]),
                        "proteoform_match": "contained",
                        "model_length_aa": model_len,
                        "is_chosen": 0,
                        "mean_plddt": "",
                        "note": "ordered protein is a subsequence of the model, but its "
                        "position was never located — run compute_residue_offsets.py",
                    }
                )
                continue
            offset = int(o["residue_offset"])
            out.append(
                {
                    **base,
                    "source": "alphafold_db",
                    "r2_key": "",
                    "s3_key": "",
                    "model_url": AFDB.format(acc=g["accession"]),
                    "proteoform_match": "contained",
                    "residue_offset": offset,
                    "model_length_aa": o["model_length_aa"],
                    "is_chosen": 1,
                    "mean_plddt": "",
                    "note": (
                        f"ordered protein is residues {o['display_first_model_residue']}-"
                        f"{o['display_last_model_residue']} of the model "
                        f"({o['ordered_length_aa']} of {o['model_length_aa']} aa); "
                        + (
                            "a clean C-terminal truncation, so model numbering already "
                            "matches the construct — show only that span"
                            if offset == 0
                            else f"renumber by {offset} and show only that span"
                        )
                    ),
                }
            )
            continue

        out.append(
            {
                **base,
                "source": "alphafold_db",
                "r2_key": "",
                "s3_key": "",
                "model_url": AFDB.format(acc=g["accession"]),
                "proteoform_match": verdict,
                "model_length_aa": model_len,
                "is_chosen": 1 if verdict == "identical" else 0,
                "mean_plddt": "",
                "note": {
                    "identical": "",
                    "different": f"AlphaFold DB models a different proteoform "
                    f"({model_len} aa vs {base['ordered_length_aa']} aa ordered)",
                }.get(verdict, "not compared to the designed CDS"),
            }
        )

    adopted = adopt_native_models(out, gap, ordered)
    write_payload(out, set(gap))
    print(f"  adopted from the native-model index: {adopted}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
