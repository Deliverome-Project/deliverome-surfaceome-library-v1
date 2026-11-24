#!/usr/bin/env python3
"""One design of every internal member, under the r8 gate, with a chosen model (v5 or v6).

v5's designs were made in several batches (adopted panel on mask r4, orthologs and start-trimmed
members on r5, the isoform-fixed 11 on r8, routed-in members without a mask). To see what a new
model changes, and only that, this designs every current internal member in one pass with
tag_site_model_v4.design: the same features (every table a v5 design read, the isoform fix as an
override), the same picking and order rule, the r8 membrane-side mask for all, and only the
model swapped. Members with no usable r8 row (both ENTPD2 members, routed in on UniProt topology) are left
to their routed-in design.

    uv run --locked python redesign_all.py --model v6   # writes tag_sites_v6_all_r8.tsv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[3] / "src"))

import membrane_mask  # noqa: E402
import tag_site_model_v4 as V4  # noqa: E402
import tag_site_model_v5 as V5  # noqa: E402
import tag_site_model_v6 as V6  # noqa: E402
from site_count_proposal import designs  # noqa: E402
from v5_region_table import EXTRA, OVERRIDE, PROVENANCE  # noqa: E402

MASK_R8 = HERE / "data" / "mask_r8" / "l1_extracellular_mask_residues.tsv"
MODELS = {"v5": V5, "v6": V6}


def frames(dec: pd.DataFrame) -> dict[str, str]:
    """sequence_sha256 -> the feature key each member's rows are filed under."""
    d = designs()
    frame = (
        d.dropna(subset=["uniprot_acc"])
        .drop_duplicates("sequence_sha256")
        .set_index("sequence_sha256")["uniprot_acc"]
        .to_dict()
    )
    for f in PROVENANCE:
        if (HERE / f).is_file():
            pv = pd.read_csv(HERE / f, sep="\t")
            sha_of = dict(zip(dec["design_id"], dec["sequence_sha256"], strict=False))
            for did, key in zip(pv["design_id"], pv["key"], strict=False):
                if did in sha_of:
                    frame.setdefault(sha_of[did], key)
    return frame


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=sorted(MODELS), required=True)
    args = ap.parse_args(argv)
    m = MODELS[args.model]
    dec = pd.read_csv(HERE / "l1_internal_decisions.tsv", sep="\t")
    frame = frames(dec)
    r8 = pd.read_csv(MASK_R8, sep="\t")
    rows = []
    for r in dec.itertuples(index=False):
        k = frame.get(r.sequence_sha256)
        if k is not None:
            rows.append({"gene": r.gene, "uniprot_acc": k, "member_id": r.member_id,
                         "design_id": r.design_id, "candidate_set": r.candidate_set,
                         "sequence_sha256": r.sequence_sha256})  # fmt: skip
    targets = pd.DataFrame(rows).drop_duplicates("uniprot_acc")
    tfile = f"data/redesign_targets_{args.model}.tsv"
    targets.to_csv(HERE / tfile, sep="\t", index=False)
    key_to_sha = dict(zip(targets["uniprot_acc"], targets["sequence_sha256"], strict=False))
    real_load, real_a2s = membrane_mask.load, membrane_mask.accession_to_sha
    membrane_mask.load = lambda: r8
    membrane_mask.accession_to_sha = lambda: {**real_a2s(), **key_to_sha}
    spec = json.loads(m.SPEC.read_text())
    out = f"tag_sites_{args.model}_all_r8.tsv"
    try:
        V4.design(
            joblib.load(m.MODEL_Z),
            joblib.load(m.MODEL_RAW),
            spec,
            np.array(spec["published_protein_ceilings_oof"]),
            np.array(spec["published_region_ceilings_oof"]),
            tfile,
            out,
            (*V5.EXTRA, *EXTRA),
            use_mask=True,
            skip_refused=True,
            min_region_len=V5.MIN_REGION_LEN,
            override_tables=tuple(OVERRIDE),
        )
    finally:
        membrane_mask.load, membrane_mask.accession_to_sha = real_load, real_a2s
    d = pd.read_csv(HERE / out, sep="\t")
    d["sequence_sha256"] = d["uniprot_acc"].map(key_to_sha)
    d["design_id"] = d["sequence_sha256"].map(
        dict(zip(dec["sequence_sha256"], dec["design_id"], strict=False))
    )
    d["model"] = args.model
    d.to_csv(HERE / out, sep="\t", index=False)
    print(f"{args.model}: {d['design_id'].nunique()} of {len(dec)} members designed -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
