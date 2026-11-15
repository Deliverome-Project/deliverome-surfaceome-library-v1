"""Where a human library reference's form, not its draw, changes the tag lane.

Compares each active L1 human `library_reference` member's designed protein with its gene's
reviewed UniProt canonical, and assigns both the placement lane `l1_tag_paths.classify` would
give them on their own DeepTMHMM v1 string (`topology.TOOL_V1`, matched by sequence hash).
A member whose lane differs from its canonical's is one where the chosen proteoform changes
where the tag can go. So is a member whose CDS starts on a non-ATG codon. The script writes
`reference_form_screen.tsv` (every reference) and prints the flagged rows.

Relation of member to canonical, from the N-terminus:
`identical`; `same_start_differs` (same first 20 residues, differs later);
`starts_downstream` (the member's residues 2-21 sit later in the canonical); `extends_upstream`
(the canonical's residues 2-21 sit later in the member); `divergent_start`; `no_canonical`.

Attribution (Becca's rule, 2026-11-10): Becca asked on 2026-11-12 to "look at where the human
form is the problem". The screen and its rules were written by the Claude Code session
l1-nterm-res1-multipass (PR #280); no person has reviewed the rules. Every flagged row is a
question for Becca, not a change. Nothing here edits a reference or a placement.

Reproducibility: DeepTMHMM v1 by default, against the snapshot pinned in
`topology.SNAPSHOT_KEY`. UniProt is pinned to the release recorded in
`uniprot_reference_screen_snapshot.json`, and the entries used are kept in
`uniprot_reference_screen.tsv` (library genes only). Every row carries the entry and sequence
versions.

    uv run --locked python screen_reference_forms.py            # offline, from the snapshot
    uv run --locked python screen_reference_forms.py --refresh  # re-fetch UniProt
"""

from __future__ import annotations

import argparse
import io
import json
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd
from l1_tag_paths import classify

from deliverome_analysis import gene_resolver, topology
from deliverome_analysis import library_design as ld
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
OUT = HERE / "reference_form_screen.tsv"
SNAPSHOT = HERE / "uniprot_reference_screen.tsv"
META = HERE / "uniprot_reference_screen_snapshot.json"
URL = (
    "https://rest.uniprot.org/uniprotkb/stream?format=tsv"
    "&fields=accession,gene_primary,version,sequence_version,ft_signal,sequence"
    "&query=(organism_id:9606)+AND+(reviewed:true)"
)
#: Residues compared to place the member against its canonical.
K = 20
IMPLEMENTED_BY = "Claude Code session l1-nterm-res1-multipass, PR #280, 2026-11-12"
#: Becca's decisions on the flagged references, 2026-11-12 ("Send all fixes"), keyed on
#: design_id. Applied by the reviewer-package session on PR #283; this table only records them.
DECIDED = 'Becca, 2026-11-12 ("Send all fixes")'
DECISIONS = {
    "l1-p1882": "approve CTG->ATG start replacement (SLC22A17)",
    "l1-p2281": "approve ACG->ATG start replacement (TRPV6)",
    "l1-p1525": "no change; OR8G3 start stays in curation",
    "l1-p0067": "replace with the signal-peptide form Q8IZF5-1 (ADGRF3)",
    "l1-p0215": "replace with the signal-peptide form (CACHD1; approved 2026-11-11)",
    "l1-p0272": "replace with the signal-peptide form (CD200R1L; approved 2026-11-11)",
    "l1-p1136": "replace with the signal-peptide form (LRRC37B; approved 2026-11-11)",
    "l1-p0578": "replace with the form starting at the signal peptide, Q9UF33 (EPHA6)",
    "l1-p1731": "replace with the form starting at the signal peptide, Q6NW40 (RGMB)",
}


def lane(topo: str) -> str:
    if not isinstance(topo, str) or not topo:
        return ""
    rung, _ = classify(topo)
    n_tm = sum(1 for i, c in enumerate(topo) if c == "M" and (i == 0 or topo[i - 1] != "M"))
    if rung == "N at residue 1":
        rung += ", multipass" if n_tm > 1 else ", single-pass (type III)"
    return rung


def relation(member: str, canon: str) -> tuple[str, int | str]:
    if not canon:
        return "no_canonical", ""
    if member == canon:
        return "identical", 0
    if member[:K] == canon[:K]:
        return "same_start_differs", 0
    at = canon.find(member[1 : K + 1])
    if at > 1:
        return "starts_downstream", at - 1
    at = member.find(canon[1 : K + 1])
    if at > 1:
        return "extends_upstream", -(at - 1)
    if member[1 : K + 1] == canon[1 : K + 1]:
        return "same_start_differs", 0  # only residue 1 differs (a non-ATG start read as L/V/T)
    return "divergent_start", ""


def _uniprot(genes: set[str], refresh: bool) -> pd.DataFrame:
    if refresh or not SNAPSHOT.exists():
        with urllib.request.urlopen(URL, timeout=300) as fh:
            text = fh.read().decode()
            release = fh.headers.get("X-UniProt-Release", "")
            rdate = fh.headers.get("X-UniProt-Release-Date", "")
        df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str).fillna("")
        df.columns = [
            "accession",
            "gene_primary",
            "entry_version",
            "sequence_version",
            "ft_signal",
            "sequence",
        ]
        ens = gene_resolver.to_ensembl(sorted(genes | set(df["gene_primary"]) - {""}))
        keep_ens = {e for g in genes for e in ens.get(g, [])}
        df = df[df["gene_primary"].map(lambda g: bool(set(ens.get(g, [])) & keep_ens))]
        df.sort_values("accession").to_csv(SNAPSHOT, sep="\t", index=False)
        META.write_text(
            json.dumps(
                {"release": release, "release_date": rdate, "fetched": str(date.today())}, indent=2
            )
            + "\n"
        )
    return pd.read_csv(SNAPSHOT, sep="\t", dtype=str).fillna("")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args(argv)

    members = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    refs = {
        c: m
        for c, m in members.items()
        if m["candidate_set"] == "library_reference" and m.get("species") == "human"
    }
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str).set_index("candidate_id")
    with open(processed.download(ld.CDS_KEY)) as fh:
        cds = json.load(fh)
    genes = {m["gene_symbol"] for m in refs.values()}
    uni = _uniprot(genes, args.refresh)
    release = json.loads(META.read_text())["release"]
    ens = gene_resolver.to_ensembl(sorted(genes | set(uni["gene_primary"])))
    by_ens: dict[str, list] = {}
    for r in uni.itertuples(index=False):
        for e in ens.get(r.gene_primary, []):
            by_ens.setdefault(e, []).append(r)

    rows, canon_seq = [], {}
    for cid, m in refs.items():
        seq = m["sequence"].rstrip("*")
        hits = {h.accession: h for e in ens.get(m["gene_symbol"], []) for h in by_ens.get(e, [])}

        # Several reviewed entries can share a gene (e.g. PRNP and its alternative-frame
        # AltPrP). Take the one that shares the member's sequence, then the longest.
        def score(h, seq=seq, gene=m["gene_symbol"]):
            return (
                seq[1 : K + 1] in h.sequence or h.sequence[1 : K + 1] in seq,
                h.gene_primary == gene,
                len(h.sequence),
            )

        h = max(hits.values(), key=score) if hits else None
        canon = h.sequence if h else ""
        if h:
            canon_seq[cid] = canon
        rel, off = relation(seq, canon)
        rows.append(
            {
                "design_id": manifest.at[cid, "design_id"],
                "candidate_id": cid,
                "gene_symbol": m["gene_symbol"],
                "transcript": m.get("transcript", ""),
                "start_codon": cds.get(manifest.at[cid, "key"], {}).get("cds", "")[:3],
                "member_aa": len(seq),
                "uniprot_accession": h.accession if h else "",
                "uniprot_release": release if h else "",
                "uniprot_entry_version": h.entry_version if h else "",
                "uniprot_sequence_version": h.sequence_version if h else "",
                "canonical_aa": len(canon),
                "canonical_signal": h.ft_signal.split(";")[0] if h else "",
                "relation": rel,
                "offset_aa": off,
                "sequence_sha256": topology.sequence_sha256(seq),
            }
        )
    out = pd.DataFrame(rows)
    mv1 = ld.topology_report(refs, tool=topology.TOOL_V1).set_index("candidate_id")
    cv1 = topology.resolve(
        canon_seq, tool=topology.TOOL_V1, cohorts=topology.cohorts_for_species("human")
    ).set_index("key")
    out["deeptmhmm_v1_member_lane"] = out["candidate_id"].map(
        lambda c: lane(mv1.at[c, "per_residue_topology"]) if bool(mv1.at[c, "found"]) else ""
    )
    out["deeptmhmm_v1_canonical_lane"] = out["candidate_id"].map(
        lambda c: (
            lane(cv1.at[c, "per_residue_topology"])
            if c in cv1.index and bool(cv1.at[c, "found"])
            else ""
        )
    )
    out["topology_version"] = out["candidate_id"].map(mv1["topology_version"])
    out["lane_differs"] = (
        (out["deeptmhmm_v1_member_lane"] != out["deeptmhmm_v1_canonical_lane"])
        & (out["deeptmhmm_v1_member_lane"] != "")
        & (out["deeptmhmm_v1_canonical_lane"] != "")
    ).astype(int)
    out["non_atg_start"] = (out["start_codon"] != "ATG").astype(int)
    flagged = (out["lane_differs"] == 1) | (out["non_atg_start"] == 1)
    out["decided_by"] = ""
    out.loc[flagged, "decided_by"] = "no decision taken: flagged for Becca"
    decided = out["design_id"].isin(DECISIONS)
    out.loc[decided, "decided_by"] = DECIDED
    out["decision"] = out["design_id"].map(DECISIONS).fillna("")
    out["implemented_by"] = IMPLEMENTED_BY
    out.to_csv(OUT, sep="\t", index=False, lineterminator="\n")

    print(f"{len(out)} human library references; UniProt {release}")
    print(out["relation"].value_counts().to_string())
    print(
        f"lane differs from canonical: {int(out['lane_differs'].sum())}; "
        f"non-ATG start: {int(out['non_atg_start'].sum())}; "
        f"canonical without v1: {int((out['deeptmhmm_v1_canonical_lane'] == '').sum())}"
    )
    cols = [
        "design_id",
        "candidate_id",
        "gene_symbol",
        "transcript",
        "start_codon",
        "member_aa",
        "canonical_aa",
        "canonical_signal",
        "relation",
        "offset_aa",
        "deeptmhmm_v1_member_lane",
        "deeptmhmm_v1_canonical_lane",
    ]
    print(out[flagged][cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
