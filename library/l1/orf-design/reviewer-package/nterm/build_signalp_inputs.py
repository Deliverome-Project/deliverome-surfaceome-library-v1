"""The sequences the 2026-11-10 SignalP 6.0 run predicts, one FASTA record per distinct protein.

Two sets, both with the p003 51-aa tandem (``GGGGS``-SpyTag003-``GGGGSGGGGS``-P-ALFA-P-``GGGGS``),
the cassette Becca chose for the N-terminal arms on 2026-11-10:

* ``new_member``: the five N-after-SP members added after the sweep (their native designed
  protein, and the tandem after each candidate site from ``nterm_positions.tsv``).
* ``igk_swap``: the members whose native signal peptide fails at every candidate (lost, or
  cleaved at or past the junction; ``nterm_build_options.tsv``). The native signal peptide,
  residues ``1..c``, is replaced by the Igκ leader plus Asp, as ``build_sp_swap_arms.py``
  builds VSIG10's arm, and the tandem goes after the same mature residue as each native
  candidate: native site ``n`` becomes ``len(leader) + 1 + (n - c)``.

Every record is named by the SHA-256 of its sequence; ``signalp_inputs.tsv`` maps each one
back to its member, arm and site. The FASTA keeps a fixed order (member ID, then arm, then
site), which matters because SignalP's numbers depend on batch composition.

**Attribution** (convention of 2026-11-10). The two sets are chosen by the Claude Code
session ``nterm-reviewer-package`` (PR #254), to answer Becca's Igκ question and close the
coverage gap; the cassette is Becca's choice (2026-11-10, "we want the tandem"). The Igκ
construction follows ``build_sp_swap_arms.py`` (decided by Becca, 2026-10-18).

    uv run --locked python build_signalp_inputs.py
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parents[1]
sys.path.insert(0, str(S4E9))
sys.path.insert(0, str(HERE))

import build_constructs as bc  # noqa: E402
import nterm_placement as npl  # noqa: E402
from build_gpcr_rescue_arms import LEADER_AA, PLUS1_AA  # noqa: E402

CASSETTE_AA = bc.CASSETTE_AA
FAILING = ("signal_peptide_lost_everywhere", "cleavage_crosses_junction_everywhere")
FASTA = HERE / "signalp_inputs.fasta"
TABLE = HERE / "signalp_inputs.tsv"


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    from deliverome_analysis.s3 import processed

    members, _ = npl.load_members(processed)
    opts = {r["candidate_id"]: r for r in _tsv(HERE / "nterm_build_options.tsv")}
    positions: dict[str, list[dict]] = {}
    for p in _tsv(HERE / "nterm_positions.tsv"):
        positions.setdefault(p["candidate_id"], []).append(p)

    rows: list[dict] = []

    def add(cid: str, arm: str, seq: str, after: int | str, native_after: int | str) -> None:
        m = members[cid]
        rows.append(
            {
                "candidate_id": cid,
                "design_id": m["design_id"],
                "gene_symbol": m["gene_symbol"],
                "set": "new_member" if arm in ("native", "tagged") else "igk_swap",
                "arm": arm,
                "after_residue": after,
                "native_after_residue": native_after,
                "length_aa": len(seq),
                "sequence_sha256": sha(seq),
                "sequence": seq,
            }
        )

    for cid in sorted(members):
        o = opts[cid]
        native = members[cid]["sequence"]
        pos = sorted(positions.get(cid, []), key=lambda p: int(p["after_residue"]))
        if o["subgroup"] == "not_repredicted":
            add(cid, "native", native, "", "")
            for p in pos:
                n = int(p["after_residue"])
                add(cid, "tagged", native[:n] + CASSETTE_AA + native[n:], n, n)
        elif o["subgroup"] in FAILING:
            c = int(o["cleavage_site"])
            igk = LEADER_AA + PLUS1_AA + native[c:]
            head = len(LEADER_AA) + len(PLUS1_AA)
            add(cid, "igk_native", igk, "", "")
            for p in pos:
                n = int(p["after_residue"])
                q = head + (n - c)
                # the mature residues before the cassette are the native ones, unchanged
                assert igk[head:q] == native[c:n]
                add(cid, "igk_tagged", igk[:q] + CASSETTE_AA + igk[q:], q, n)

    with TABLE.open("w", newline="") as fh:
        cols = [k for k in rows[0] if k != "sequence"]
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)
    seen, n_written = set(), 0
    with FASTA.open("w") as fh:
        for r in rows:
            if r["sequence_sha256"] in seen:
                continue
            seen.add(r["sequence_sha256"])
            fh.write(f">{r['sequence_sha256']}\n{r['sequence']}\n")
            n_written += 1
    by = {}
    for r in rows:
        by.setdefault((r["set"], r["arm"]), set()).add(r["candidate_id"])
    for (s, a), ids in sorted(by.items()):
        print(
            f"{s:11} {a:11} members {len(ids):4}  sequences "
            f"{sum(1 for r in rows if r['set'] == s and r['arm'] == a)}"
        )
    print(f"{len(rows)} rows, {n_written} distinct sequences -> {FASTA.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
