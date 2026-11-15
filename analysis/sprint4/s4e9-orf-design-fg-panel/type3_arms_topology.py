"""Classify DeepTMHMM v1 re-predictions of the three type III arms, and the subgroups.

Input is the 3-line output of the accessible-surfaceome DeepTMHMM 1.0.24 Modal app run on
``type3_arms.fasta`` (written by ``type3_nterm.py``): per member, the native protein, the
cassette after Met1 ("tag"), and Igk leader + Asp + cassette + native from Met1
("igk_tag", leader + Asp + Met1 + cassette + native 2..end). The 3-line file is published to
``processed://sprint4/s4e9-orf-design-fg-panel/type3-nterm/type3_arms.v1.3line`` and
read from there unless ``--three-line`` points at a local copy.

For each arm it records whether the cassette ends up extracellular, whether the native TM
is still called, and the predicted signal-peptide end (Igk is 20 aa). It then assigns each
member to a design subgroup from ``type3_nterm_features.tsv`` (see ``subgroup``) and writes
``type3_nterm_calls.tsv``.

    uv run --locked python type3_arms_topology.py
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
FEATURES = HERE / "type3_nterm_features.tsv"
ARMS_FASTA = HERE / "type3_arms.fasta"
OUT = HERE / "type3_nterm_calls.tsv"
S3_KEY = "sprint4/s4e9-orf-design-fg-panel/type3-nterm/type3_arms_rekey_6bf341f2.v1.3line"
CASSETTE_LEN = 51
IGK_LEN = 20  # leader; +1 Asp follows

#: Native N-tail at or below which the cassette moves the TM from the N-terminal
#: signal-anchor regime into the internal one (Kocik et al. 2012: insertion behaviour
#: changes once the hydrophilic N-domain exceeds ~20 residues).
SHORT_TAIL = 20
#: Native N-tail at or above which the cassette is a minority of what the TM already
#: translocates natively (51 aa is under ~40 % of the tail).
LONG_TAIL = 130


def read_three_line(text: str) -> dict[str, tuple[str, str, str]]:
    """{name: (header label, sequence, topology)} from DeepTMHMM .3line output."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    out = {}
    for i in range(0, len(lines), 3):
        head, seq, topo = lines[i : i + 3]
        name, _, label = head[1:].partition(" | ")
        out[name.split("|")[1] if name.startswith("sp|") else name] = (label, seq, topo)
    return out


def read_fasta(text: str) -> dict[str, str]:
    """{candidate_id__arm: sequence} from ``type3_arms.fasta`` (``>sp|<id>__<arm>|<gene>``)."""
    out, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[1]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


def arm_call(arm: str, topo: str, native_tm1: int) -> dict:
    """Where the cassette lands and whether TM1 survives, for one arm's v1 string."""
    sp = len(topo) - len(topo.lstrip("S"))
    if arm == "native":
        cas = None
    elif arm == "tag":
        cas = topo[1 : 1 + CASSETTE_LEN]
    else:
        cas = topo[IGK_LEN + 2 : IGK_LEN + 2 + CASSETTE_LEN]
    shift = {"native": 0, "tag": CASSETTE_LEN, "igk_tag": IGK_LEN + 1 + CASSETTE_LEN}[arm]
    n_tm = sum(1 for i, c in enumerate(topo) if c == "M" and (i == 0 or topo[i - 1] != "M"))
    tm1 = native_tm1 - 1 + shift
    win = topo[max(0, tm1 - 5) : tm1 + 15]
    if cas is None:
        cassette = ""
    elif cas.count("O") >= 0.8 * len(cas):
        cassette = "out"
    elif cas.count("I") >= 0.8 * len(cas):
        cassette = "in"
    else:
        cassette = "mixed"
    return {
        f"{arm}_label": "",
        f"{arm}_sp_aa": sp,
        f"{arm}_n_tm": n_tm,
        f"{arm}_native_tm_kept": int("M" in win),
        f"{arm}_cassette": cassette,
        f"{arm}_n_term_state": topo[sp] if sp < len(topo) else "",
    }


def subgroup(f: dict) -> str:
    """Design subgroup, in order of precedence.

    * ``no_tm``: v1 calls no TM at all; the rung admits n_tm == 0 by construction.
    * ``sp_in_member``: the member carries its gene's UniProt signal peptide from residue 1
      (v1 did not call it), so it belongs on the N-after-SP rung.
    * ``sp_missing``: the gene has a UniProt signal peptide and the ordered form starts
      downstream of it (alternative start, ortholog model) or replaces it.
    * ``dual_topology``: MRAP / MRAP2, whose function rests on both orientations.
    * ``short_tail`` / ``mid_tail`` / ``long_tail``: genuine type III by native N-tail.
    """
    if f["n_tm_v1"] == "0":
        return "no_tm"
    if f.get("uniprot_signal"):
        return "sp_in_member" if f.get("member_has_canonical_sp") == "1" else "sp_missing"
    if f["gene_symbol"] in {"MRAP", "MRAP2"}:
        return "dual_topology"
    tail = int(f["n_tail_native_aa"])
    if tail <= SHORT_TAIL:
        return "short_tail"
    if tail >= LONG_TAIL:
        return "long_tail"
    return "mid_tail"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--three-line", type=Path, default=None)
    args = ap.parse_args(argv)
    path = args.three_line or Path(processed.download(S3_KEY))
    # Match on SEQUENCE, not on the FASTA name: member IDs were renumbered (l1-pNNNN ->
    # l1-mNNNN) after the 3-line file was produced, and the proteins did not change.
    by_seq = {
        seq: (label, seq, topo) for label, seq, topo in read_three_line(path.read_text()).values()
    }
    calls = {}
    for name, seq in read_fasta(ARMS_FASTA.read_text()).items():
        assert seq in by_seq, f"{name}: arm sequence has no v1 prediction"
        calls[name] = by_seq[seq]

    with FEATURES.open() as fh:
        feats = list(csv.DictReader(fh, delimiter="\t"))
    rows = []
    for f in feats:
        cid = f["candidate_id"]
        row = {
            "candidate_id": cid,
            "design_id": f.get("design_id", ""),
            "candidate_set": f["candidate_set"],
            "species": f["species"],
            "gene_symbol": f["gene_symbol"],
            "subgroup": subgroup(f),
            "n_tail_native_aa": f.get("n_tail_native_aa", ""),
            "delta_charge_native": f.get("delta_charge_native", ""),
            "tm1_len": f.get("tm1_len", ""),
            "tm1_kd_mean": f.get("tm1_kd_mean", ""),
            "v2_label": f.get("v2_label", ""),
            "v2_n_term": f.get("v2_n_term", ""),
        }
        tm1 = int(f["tm1_start"]) if f.get("tm1_start") else 1
        for arm in ("native", "tag", "igk_tag"):
            label, _seq, topo = calls[f"{cid}__{arm}"]
            c = arm_call(arm, topo, tm1)
            c[f"{arm}_label"] = label.strip()
            row.update(c)
        rows.append(row)

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT.name}: {len(rows)} members")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
