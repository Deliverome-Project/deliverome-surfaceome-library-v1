"""Native arms with a weak predicted signal peptide, and the Igκ form that would replace each.

Becca, 2026-11-14: "we are relying more on igk now for either weak predicted SP like VSIG10 or
all GPCRs". Definition (thresholds PROPOSED by nterm-reviewer-package, adopted by the
coordinator under Becca's delegation on 2026-11-15, not yet confirmed by Becca):

- W1: the native SignalP 6.0 cleavage-site probability is below ``W_THRESHOLD`` (VSIG10's own
  weakness: 0.35); or
- W2: the built native arm's tagged cleavage-site probability is below ``W_THRESHOLD``.

Rejected: W1 alone, W2 alone, and the native SP probability (none of the 1,090 native arms falls
below 0.5, so it does not separate cases). Would change: Becca setting another cut, or bench
cleavage data.

For each member the Igκ leader replaces residues 1..c, c being the later of the DeepTMHMM v1 and
SignalP 6.0 calls (so it deletes at most a few mature residues rather than leaving signal-peptide
residues), and the cassette stays after the same native residue. Where the two calls differ by
more than ``MAX_CUT_SPREAD`` residues the cut is too uncertain to delete on: no swap, flagged for
Becca. The swap stands only if SignalP keeps the leader's cleavage on the Igκ form and the form
passes refold (or the disposition rule); otherwise the member keeps its native arm, flagged.

    uv run --locked python weak_sp.py fasta     # Igκ forms not yet scored -> weak_sp_new.fasta
    uv run --locked python weak_sp.py table     # weak_sp.tsv
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(S4E9))

import build_constructs as bc  # noqa: E402
import nterm_overnight as no  # noqa: E402
import nterm_placement as npl  # noqa: E402
from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA  # noqa: E402

W_THRESHOLD = 0.5
MAX_CUT_SPREAD = 5
HEAD = len(LEADER_AA) + len(PLUS1_AA)
OUT = HERE / "weak_sp.tsv"
NEW_FASTA = HERE / "weak_sp_new.fasta"


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def rows() -> list[dict]:
    from deliverome_analysis.s3 import processed

    snap = npl.load_snapshot(processed)
    calls = no.prior_calls(processed)
    members, _ = npl.load_members(processed)
    by_design = {m["design_id"]: m for m in members.values()}
    census = {r["design_id"]: r for r in _tsv(HERE / "nterm_members.tsv")}
    by_key: dict[str, list[dict]] = {}
    for (_src, key), r in snap["signalp"].items():
        by_key.setdefault(key, []).append(r)
    out = []
    for d in _tsv(HERE / "nterm_decisions.tsv"):
        if d["arm"] != "native_site":
            continue
        did = d["design_id"]
        seq = by_design[did]["sequence"]
        c = census.get(did, {})
        native = next(
            (
                r
                for r in by_key.get(c.get("source_key", ""), []) + by_key.get(sha(seq), [])
                if r["prediction"] == "SP"
            ),
            None,
        ) or calls.get(sha(seq))
        ncp = float(native["cleavage_probability"] or 0)
        sp_cut = int(native["cleavage_site"])
        v1_cut = c.get("deeptmhmm_v1_cleavage") or by_design[did]["v1_cleavage"] or ""
        tcp = float(d["tagged_cleavage_probability"] or 0)
        w1, w2 = ncp < W_THRESHOLD, tcp < W_THRESHOLD
        if not (w1 or w2):
            continue
        cut = int(d["cleavage_site"])
        n = int(d["native_after_residue"])
        spread = abs(int(v1_cut) - sp_cut) if v1_cut != "" else 0
        parent = LEADER_AA + PLUS1_AA + seq[cut:]
        q = HEAD + (n - cut)
        tagged = parent[:q] + bc.CASSETTE_AA + parent[q:]
        pc, tc = calls.get(sha(parent)), calls.get(sha(tagged))
        verdict = npl.verdict(pc, tc, q) if pc and tc else ""
        out.append(
            {
                "design_id": did,
                "native_tagged_sha256": d["tagged_sha256"],
                "candidate_id": d["candidate_id"],
                "gene_symbol": d["gene_symbol"],
                "w1_native_cleavage_probability_lt": int(w1),
                "w2_tagged_cleavage_probability_lt": int(w2),
                "native_sp_probability": native.get("sp_probability", ""),
                "native_cleavage_probability": ncp,
                "tagged_cleavage_probability": tcp,
                "signalp_cut": sp_cut,
                "deeptmhmm_v1_cut": v1_cut,
                "cut_spread": spread,
                "igk_replaces": f"1-{cut}",
                "native_after_residue": n,
                "insert_after_residue": q,
                "parent_sha256": sha(parent),
                "tagged_sha256": sha(tagged),
                "parent_sequence": parent,
                "tagged_sequence": tagged,
                "signalp_igk_parent": f"{pc['prediction']} {pc['cleavage_site']}" if pc else "",
                "signalp_igk_tagged": (
                    f"{tc['prediction']} {tc['cleavage_site']} p={tc['cleavage_probability']}"
                    if tc
                    else ""
                ),
                "signalp_verdict": verdict,
                "status": "flag_for_becca: the two cut calls differ by more than "
                f"{MAX_CUT_SPREAD} residues"
                if spread > MAX_CUT_SPREAD
                else "",
            }
        )
    return out


def main() -> int:
    rs = rows()
    if sys.argv[1:2] == ["fasta"]:
        from deliverome_analysis.s3 import processed

        calls = no.prior_calls(processed)
        seqs = {}
        for r in rs:
            for k in ("parent", "tagged"):
                if r[f"{k}_sha256"] not in calls:
                    seqs[r[f"{k}_sha256"]] = r[f"{k}_sequence"]
        NEW_FASTA.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        print(f"{len(rs)} weak members; {len(seqs)} sequences to score -> {NEW_FASTA.name}")
        return 0
    cols = [k for k in rs[0] if not k.endswith("_sequence")]
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rs)
    print(f"{len(rs)} weak members -> {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
