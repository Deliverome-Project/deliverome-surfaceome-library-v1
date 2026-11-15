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
from pathlib import Path

from deliverome_analysis import library_design as ld

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
    from deliverome_analysis.s3 import processed

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
    out: dict[str, float] = {}
    # The committed fold summaries (every run, merged), not a local cache: a fold run on
    # another machine or session is otherwise invisible here.
    for path in (
        HERE / "predicted_structures/fold_summary.tsv",
        HERE / "predicted_structures/fold_summary_mismatched.tsv",
        Path.home() / ".cache/deliverome-analysis/predicted_structures_mismatched/fold_summary.tsv",
    ):
        if not path.exists():
            continue
        with path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("mean_plddt"):
                    out[r["name"]] = float(r["mean_plddt"])
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
                "note": "" if verdict == "identical" else "not compared to the designed CDS",
            }
        )

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(out)

    import collections

    chosen = [r for r in out if r["is_chosen"] == 1]
    print(f"{len(out)} rows for {len({r['candidate_id'] for r in out})} members")
    print(
        f"  chosen   {len(chosen)}  by source {dict(collections.Counter(r['source'] for r in chosen))}"
    )
    alt = [r for r in out if r["is_chosen"] == 0]
    print(
        f"  alternate {len(alt)}  by match {dict(collections.Counter(r['proteoform_match'] for r in alt))}"
    )
    missing = {r["candidate_id"] for r in out} ^ set(gap)
    print(f"  members with no row at all: {len(missing)}")
    print(f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
