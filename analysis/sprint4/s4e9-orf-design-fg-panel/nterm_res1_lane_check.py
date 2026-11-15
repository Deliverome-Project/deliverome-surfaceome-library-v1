"""Is "N at residue 1, multipass" the right lane for each member on it?

`l1_tag_paths.py` puts a member on this rung from its own DeepTMHMM v1 string (`deeptmhmm_v1`,
matched by sequence hash): v1 says residue 1 is extracellular, calls no signal peptide, and calls
more than one helix. **DeepTMHMM v1 is the topology of record and the default lane call**
(Becca, 2026-11-11). Every other source here is a cross-check: it can only flag a disagreement
with v1 (`crosscheck_disagrees`, naming the source) and send the member to literature review.
No cross-check replaces the v1 call in any column. The review's verdicts are a committed input
(`nterm_res1_lane_review.tsv`), not computed here, and they are recommendations until Becca
confirms them.

A member is in the right lane only if its ordered protein has, at residue 1:

1. an extracellular N-terminus;
2. no cleavable signal peptide; and
3. at least two TM helices.

**Cross-checks, per condition.**

| Source | Conditions | When it is read |
|---|---|---|
| UniProt topological domains and SIGNAL feature, with evidence codes | 1, 2, 3 | 1 and 2 only when the member's first 20 residues are the canonical's; 3 only when the whole member is the canonical |
| DeepTMHMM v2 on the member's own sequence (sha256 match) | 1, 2, 3 | always, when a v2 row exists |
| SignalP 6 on the member's first 70 residues | 2 | when an existing run saw the same 70 residues (SignalP truncates every input to 70, so that prefix decides its call exactly) |
| The member's chosen AlphaFold DB / ColabFold model | 1, 3 | always, when a model is chosen: membrane crossings counted on the model, and the cytoplasmic face from the positive-inside rule |

**v2 alphabet.** `nterm_res1_multipass.py` read only v2's leading `S`, because v2's other state
letters are not documented in this repo. The topology snapshot carries the documentation
implicitly: every v2 row also stores `per_residue_topology`, v2 projected onto the v1 alphabet
(I / O / M / S / B). `v2_alphabet()` derives the letter-by-letter projection from all v2 rows,
refuses to run unless each letter projects to exactly one v1 state, and checks the derived
helix count against the snapshot's own `tm_helix_count`. Sides are read through that
projection. Three v2 letters project to `M`; only one is a TM helix (it reproduces
`tm_helix_count` exactly), the others are reentrant and interfacial segments such as GPCR
helix 8.

**Structure check.** A planar membrane is fitted to the member's model: the normal is the mean
axis of the helices that v1 or v2 places (positions only, never their sides), using C-alpha atoms
at pLDDT >= 70; the mid-plane passes through their centroid. The same construction as
`deliverome_analysis.membrane_side.fit_helix_axis` on the open PR #291; it is written out here
because that module is not on `main`. Switch to it once it lands. Then, independently of any
topology string:

* **crossings**: walking the confident residues (pLDDT >= 50) in order, the number of times the
  chain passes from one side of the slab (|z| > 10 A) to the other. Reentrant loops return to
  the side they left and are not counted.
* **cytoplasmic face**: Lys + Arg at the helix flanks (within 10 residues of a residue in the
  membrane core, 5-30 A from the mid-plane) on each face; the face with more is cytoplasmic
  (positive-inside rule). The N-terminus is out when it sits on the other face.

**Pinned inputs** (Becca, 2026-11-11: pin every input so we can tell when a rerun is needed).
The script runs offline from snapshots committed beside it, and writes each member's
fingerprints into every row of the output:

| Input | Snapshot | Fingerprint columns |
|---|---|---|
| UniProt entries | `nterm_res1_lane_uniprot.tsv` (`--refresh-uniprot` re-fetches) | `uniprot_release` (REST `X-UniProt-Release`), `uniprot_entry_version`, `uniprot_sequence_version`, `uniprot_sequence_sha256` |
| Topology | `topology.SNAPSHOT_KEY` (pinned parquet) | `topology_snapshot_key`, `deeptmhmm_v1_topology_version`, `deeptmhmm_v1_tool_version`, `deeptmhmm_v2_topology_version`, `deeptmhmm_v2_tool_version`, `sequence_sha256` |
| SignalP 6 | `nterm_res1_signalp6.tsv` (`--signalp-dir` rebuilds) | `signalp6_version`, `signalp6_source` (run/shard/input id) |
| Structure models | `library_structure_payload.tsv` (chosen model per member) | `model_version` (AFDB model version or ColabFold release), `model_sha256` |

`--check-drift` re-fetches only the UniProt fingerprints, recomputes each member's sequence hash
and topology version from the manifest and the pinned snapshot, lists the members whose UniProt
entry or sequence version, sequence hash, topology version or rung membership changed, and exits
nonzero if any did.

**Attribution** (Becca's rule, 2026-11-10: note the agent behind every manual fix and its
reasoning). Nothing here overrides a value or carries a hand-picked member list: the members
come from `l1_tag_paths.tsv` and the manifest, so a rerun after the Compara 116 ortholog panel
(PR #283) lands checks the new members. The flag rules and thresholds below are analysis
choices of the Claude Code session l1-nterm-lane-check (branch
`claude/l1-nterm-lane-check-98e494`), not yet reviewed by a person; each constant's comment
gives its reason. Verdicts are recommendations until Becca confirms them.

    uv run --locked python nterm_res1_lane_check.py                    # offline, from snapshots
    uv run --locked python nterm_res1_lane_check.py --check-drift      # which members to rerun
    uv run --locked python nterm_res1_lane_check.py --refresh-uniprot  # re-pin UniProt
    # rebuild the SignalP 6 snapshot from run outputs pulled off the signalp6-runs volume:
    uv run --locked python nterm_res1_lane_check.py --signalp-dir <dir>
"""

from __future__ import annotations

import argparse
import collections
import datetime
import hashlib
import io
import itertools
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from deliverome_analysis import library_design as ld
from deliverome_analysis import topology
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
TAG_PATHS = HERE / "l1_tag_paths.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
#: Per-member features from nterm_res1_multipass.py (GPCR status, Tedman class and intensity,
#: UniProt accession). Must cover exactly the rung's members; rerun that script first if not.
FEATURES = HERE / "nterm_res1_multipass.tsv"
#: UniProt entries for the rung's accessions, with the release and per-entry versions they were
#: fetched at. Written by --refresh-uniprot; read offline otherwise.
UNIPROT_SNAPSHOT = HERE / "nterm_res1_lane_uniprot.tsv"
UNIPROT_FIELDS = "accession,gene_primary,ft_signal,ft_topo_dom,version,sequence_version,sequence"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search?format=tsv&size=500&fields="
STRUCTURES = HERE / "library_structure_payload.tsv"
#: SignalP 6 calls for the rung's members, keyed on their first 70 residues. Built by
#: --signalp-dir from the signalp6-runs Modal volume (runs sp6_2026_09_28,
#: sp6_ortho_2026_09_28, sp6_2026_10_18); committed so a rerun needs no Modal access.
SIGNALP = HERE / "nterm_res1_signalp6.tsv"
#: Literature verdicts for the flagged members (reviewed input, not computed).
REVIEW = HERE / "nterm_res1_lane_review.tsv"
SCREEN_OUT = HERE / "nterm_res1_lane_screen.tsv"
OUT = HERE / "nterm_res1_lane_check.tsv"
RUNG = "N at residue 1, multipass"

AF_CACHE = Path.home() / ".cache" / "deliverome-analysis" / "afdb"
COLABFOLD_PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/"

IMPLEMENTED_BY = "Claude Code session l1-nterm-lane-check, branch claude/l1-nterm-lane-check-98e494"
DECIDED_BY_SCREEN = (
    "no decision taken: DeepTMHMM v1 call of record (Becca, 2026-11-11), no cross-check "
    "disagrees; not yet reviewed by Becca"
)

#: SignalP truncates every input to its first 70 residues (signalp/utils.py), so two
#: sequences with the same first 70 residues get the same call.
SIGNALP_WINDOW = 70
#: Residues compared to decide whether a member starts where its UniProt canonical does
#: (same as nterm_res1_multipass.N_START).
N_START = 20
#: SignalP 6 SP probability above which a signal peptide is called (SignalP's own cutoff).
SP_PROB = 0.5
#: pLDDT for C-alpha atoms that define the helix axis (as membrane_side.FIT_MIN_PLDDT).
FIT_PLDDT = 70.0
#: pLDDT below which a residue is ignored when counting crossings: very-low-confidence loops
#: are often modelled as extended ribbons through the slab.
WALK_PLDDT = 50.0
#: A helix contributes to the axis only with this many confident C-alpha atoms.
MIN_HELIX_CA = 10
#: Distance from the mid-plane (A) that counts as being on a side when walking the chain.
#: The hydrocarbon core is ~30 A thick, so +/-10 A is well inside each face.
SIDE_Z = 10.0
#: Membrane core used to find helix flanks for the positive-inside count (A).
CORE_Z = 12.0
#: Flank band either side of the mid-plane for the K/R count (A).
FLANK_Z = (5.0, 30.0)
#: Sequence distance from a core residue within which a K/R counts as a flank.
FLANK_RES = 10
#: K/R difference between faces below which the structure does not call a side.
KR_MARGIN = 3


# --- v2 ------------------------------------------------------------------------------


def v2_alphabet() -> tuple[dict[str, str], str]:
    """(letter -> v1 state, TM-helix letter), derived from the snapshot and checked.

    Every v2 row carries both its own string and the projection onto the v1 alphabet. A letter
    that projects to two different states would make the projection ambiguous; stop then.
    Three letters project to v1's ``M``; the TM-helix letter is the one whose run count
    reproduces the snapshot's ``tm_helix_count``. The others are reentrant and interfacial
    segments (GPCR helix 8 is one), which do not cross the membrane.
    """
    snap = topology.load_snapshot()
    v2 = snap[snap["tool_version"].astype(str).str.startswith(topology.V2_TOOL_PREFIX)]
    v2 = v2.dropna(subset=["dtm2_topology_string", "per_residue_topology"])
    seen: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for a, b in zip(v2["dtm2_topology_string"], v2["per_residue_topology"], strict=True):
        if len(a) == len(b):
            for x, y in zip(a, b, strict=True):
                seen[x][y] += 1
    ambiguous = {k: dict(c) for k, c in seen.items() if len(c) > 1}
    if ambiguous:
        raise SystemExit(f"v2 letters project to more than one v1 state: {ambiguous}")
    proj = {k: next(iter(c)) for k, c in seen.items()}
    sample = v2.dropna(subset=["tm_helix_count"])
    agree = {}
    for letter in (k for k, v in proj.items() if v == "M"):
        runs = sample["dtm2_topology_string"].map(lambda s, L=letter: len(re.findall(f"{L}+", s)))
        agree[letter] = float((runs == sample["tm_helix_count"]).mean())
    helix = max(agree, key=agree.get)
    print(
        f"v2 alphabet: {len(proj)} letters; TM-helix letter '{helix}' reproduces tm_helix_count "
        f"for {agree[helix]:.2%} of {len(sample)} v2 rows; other M-projecting letters "
        + ", ".join(f"'{k}' {v:.1%}" for k, v in agree.items() if k != helix)
    )
    if agree[helix] < 0.99:
        raise SystemExit("v2 alphabet decoding does not reproduce the snapshot's helix counts")
    return proj, helix


def v2_features(s: str, proj: dict[str, str], helix_letter: str) -> dict:
    v1 = "".join(proj[c] for c in s)
    sp = len(v1) - len(v1.lstrip("S"))
    helix = "".join("H" if c == helix_letter else "-" for c in s)
    n_side = next((c for c in v1[sp:] if c in "IO"), "")
    c_side = next((c for c in reversed(v1) if c in "IO"), "")
    return {
        "v2_signal_peptide_aa": sp,
        "v2_n_side": {"O": "out", "I": "in"}.get(n_side, ""),
        "v2_tm": len(re.findall("H+", helix)),
        "v2_c_side": {"O": "out", "I": "in"}.get(c_side, ""),
        "_v2_helix_mask": helix,
    }


# --- SignalP 6 -------------------------------------------------------------------------


def build_signalp_snapshot(run_dir: Path, members: dict[str, str]) -> None:
    """Collect SignalP 6 calls made on the same first 70 residues as a member.

    ``run_dir`` holds <run>/<shard>/prediction_results.txt as pulled from the signalp6-runs
    Modal volume. The runs key their inputs two ways and keep only the mature sequence of
    predicted signal peptides, so the input sequence is recovered per key:

    * a UniProt accession (sp6_2026_09_28): the input was that entry's sequence. Only
      accessions in the UniProt snapshot are resolvable, and the call is used only if that
      sequence's first 70 residues are the member's.
    * a sha256 (sp6_2026_10_18): matched to the sha256 of the member's ordered sequence.
    """
    up = load_uniprot()
    acc_prefix = dict(zip(up["accession"], up["sequence"].str[:SIGNALP_WINDOW], strict=True))
    sha_prefix = {topology.sequence_sha256(s): s[:SIGNALP_WINDOW] for s in members.values()}
    want: dict[str, list[str]] = collections.defaultdict(list)
    for cid, seq in members.items():
        want[seq[:SIGNALP_WINDOW]].append(cid)
    by_prefix: dict[str, dict] = {}
    for pred in sorted(run_dir.glob("*/*/prediction_results.txt")):
        shard = pred.parent
        lines = pred.read_text().splitlines()
        version = re.search(r"SignalP-([\d.]+)", lines[0]) if lines else None
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            prefix = acc_prefix.get(f[0]) or sha_prefix.get(f[0])
            if prefix in want and prefix not in by_prefix:
                cs = re.search(r"CS pos: (\d+)-", f[4]) if len(f) > 4 else None
                by_prefix[prefix] = {
                    "signalp6_prediction": f[1],
                    "signalp6_sp_prob": float(f[3]),
                    "signalp6_cs_after": int(cs.group(1)) if cs else "",
                    "signalp6_version": f"SignalP {version.group(1)}" if version else "",
                    "signalp6_source": f"{shard.parent.name}/{shard.name}/{f[0]}",
                }
    rows = [
        {
            "design_id": design_ids()[cid],
            "candidate_id": cid,
            "prefix70_sha256": topology.sequence_sha256(p),
            **by_prefix[p],
        }
        for p, cids in want.items()
        if p in by_prefix
        for cid in cids
    ]
    pd.DataFrame(rows).sort_values("candidate_id").to_csv(
        SIGNALP, sep="\t", index=False, lineterminator="\n"
    )
    print(f"wrote {SIGNALP.name}: {len(rows)} of {len(members)} members")


def design_ids() -> dict[str, str]:
    """candidate_id -> design_id from the manifest.

    design_id is stable across the candidate_id renumberings (l1-pNNNN -> l1-mNNNN, and PR
    #283's ortholog adoption), so every table this script writes or reads carries it
    (Becca, 2026-11-12).
    """
    man = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    return dict(zip(man["candidate_id"], man["design_id"], strict=True))


# --- UniProt ---------------------------------------------------------------------------


def fetch_uniprot(accessions: list[str]) -> pd.DataFrame:
    """Entries for these accessions, with the UniProt release they were served from."""
    frames, releases = [], set()
    for i in range(0, len(accessions), 100):
        q = urllib.parse.quote(" OR ".join(f"accession:{a}" for a in accessions[i : i + 100]))
        with urllib.request.urlopen(
            f"{UNIPROT_SEARCH}{UNIPROT_FIELDS}&query={q}", timeout=120
        ) as fh:
            releases.add(fh.headers.get("X-UniProt-Release", ""))
            frames.append(pd.read_csv(io.StringIO(fh.read().decode()), sep="\t", dtype=str))
    if len(releases) != 1:
        raise SystemExit(f"UniProt release changed during the fetch: {sorted(releases)}")
    df = pd.concat(frames).rename(
        columns={
            "Entry": "accession",
            "Gene Names (primary)": "gene_primary",
            "Signal peptide": "ft_signal",
            "Topological domain": "ft_topo_dom",
            "Entry version": "entry_version",
            "Sequence version": "sequence_version",
            "Sequence": "sequence",
        }
    )
    df = df.fillna("").drop_duplicates("accession").sort_values("accession")
    df["sequence_sha256"] = df["sequence"].map(topology.sequence_sha256)
    df["uniprot_release"] = releases.pop()
    df["fetched"] = datetime.date.today().isoformat()
    missing = sorted(set(accessions) - set(df["accession"]))
    if missing:
        print(f"UniProt returned no entry for {len(missing)} accession(s): {missing}")
    return df


def rung_accessions() -> list[str]:
    feats = pd.read_csv(FEATURES, sep="\t", dtype=str).fillna("")
    return sorted({a for a in feats["uniprot_accession"] if a})


def load_uniprot() -> pd.DataFrame:
    if not UNIPROT_SNAPSHOT.exists():
        raise SystemExit(f"{UNIPROT_SNAPSHOT.name} missing; run with --refresh-uniprot")
    return pd.read_csv(UNIPROT_SNAPSHOT, sep="\t", dtype=str).fillna("")


#: A single-residue domain is written ``TOPO_DOM 1;`` (no range), and some domains carry no
#: evidence code; both must parse, or the second domain is mistaken for the first.
TOPO = re.compile(r'TOPO_DOM (\d+)(?:\.\.(\d+))?; /note="([^"]+)"(?:; /evidence="([^"]+)")?')


def uniprot_features(seq: str, rec: dict | None, helix_letter: str) -> dict:
    """UniProt canonical features, read only where they describe this member's residues.

    Canonical TM1 locates an isoform's start relative to the first helix. UniProt gives it
    as the end of the first topological domain; where the entry has none (GPM6B), it is
    taken from DeepTMHMM v2 on the canonical's own sequence, recorded in
    ``up_canonical_tm1_source``.
    """
    out = {
        "up_same_start": 0,
        "up_identical": 0,
        "up_start_in_canonical": "",
        "up_canonical_tm1_start": "",
        "up_canonical_tm1_source": "",
        "up_signal_end": "",
        "up_n_side": "",
        "up_n_evidence": "",
        "up_tm": "",
        "up_c_side": "",
        "uniprot_release": "",
        "uniprot_entry_version": "",
        "uniprot_sequence_version": "",
        "uniprot_sequence_sha256": "",
    }
    if not rec:
        return out
    out["uniprot_release"] = rec["uniprot_release"]
    out["uniprot_entry_version"] = rec["entry_version"]
    out["uniprot_sequence_version"] = rec["sequence_version"]
    out["uniprot_sequence_sha256"] = rec["sequence_sha256"]
    can = rec["sequence"]
    doms = [
        (
            int(a),
            int(b or a),
            note,
            "experimental" if "ECO:0000269" in ev else ("inferred" if ev else "no evidence code"),
        )
        for a, b, note, ev in TOPO.findall(rec["ft_topo_dom"])
    ]
    m = re.search(r"SIGNAL 1\.\.(\d+)", rec["ft_signal"])
    can_sp = int(m.group(1)) if m else 0
    pos = can.find(seq[:N_START])
    out["up_start_in_canonical"] = pos + 1 if pos >= 0 else ""
    out["up_signal_end"] = can_sp
    if doms:
        out["up_canonical_tm1_start"] = doms[0][1] + 1
        out["up_canonical_tm1_source"] = "UniProt TOPO_DOM"
    else:
        hit = topology.for_sequence(can, tool=topology.TOOL_V2)
        if hit is not None and helix_letter in (hit.get("dtm2_topology_string") or ""):
            out["up_canonical_tm1_start"] = hit["dtm2_topology_string"].index(helix_letter) + 1
            out["up_canonical_tm1_source"] = "DeepTMHMM v2 on the canonical"
    same = can.startswith(seq[:N_START])
    out["up_same_start"] = int(same)
    out["up_identical"] = int(can == seq)
    if same and doms:
        first = doms[0]
        # A domain starting after the signal peptide still describes residue 1's fate only if
        # there is no signal peptide; record its side either way and let the SP column speak.
        out["up_n_side"] = side_word(first[2])
        out["up_n_evidence"] = first[3]
    if out["up_identical"] and doms:
        sides = [side_word(d[2]) for d in doms]
        out["up_tm"] = sum(1 for a, b in itertools.pairwise(sides) if a and b and a != b)
        out["up_c_side"] = sides[-1]
    return out


def side_word(note: str) -> str:
    n = note.lower()
    if n.startswith("cytoplasmic"):
        return "in"
    if n.startswith(("extracellular", "lumenal", "vacuolar")):
        return "out"
    return ""


# --- structure -------------------------------------------------------------------------


def read_ca(text: str) -> tuple[str, np.ndarray, np.ndarray]:
    three = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
        "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M",
    }  # fmt: skip
    seq, xyz, b = [], [], []
    for line in text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in " A":
            seq.append(three.get(line[17:20], "X"))
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            b.append(float(line[60:66]))
    return "".join(seq), np.array(xyz), np.array(b)


def model_version(row: dict) -> str:
    """AFDB model version from the file name (``..._v6.pdb``) or the ColabFold release."""
    if row["source"] == "alphafold_db" and row["model_url"]:
        m = re.search(r"model_(v\d+)\.pdb$", row["model_url"])
        return f"AFDB {m.group(1)}" if m else "AFDB"
    return f"ColabFold {row['release']}" if row.get("release") else row["source"]


def model_text(row: dict) -> str | None:
    if row["source"] == "alphafold_db" and row["model_url"]:
        AF_CACHE.mkdir(parents=True, exist_ok=True)
        path = AF_CACHE / row["model_url"].rsplit("/", 1)[-1]
        for attempt in range(6):
            if path.exists():
                break
            try:
                with urllib.request.urlopen(row["model_url"], timeout=60) as fh:
                    path.write_bytes(fh.read())
            except urllib.error.HTTPError as e:
                # AlphaFold DB rate-limits bulk fetches with 429; back off and retry.
                if e.code != 429 or attempt == 5:
                    raise
                time.sleep(5 * 2**attempt)
        return path.read_text()
    if row["r2_key"]:
        key = COLABFOLD_PREFIX + row["r2_key"].split("/", 1)[1]
        # A chosen model can be listed before its file is published (MC4R's adopted cyno
        # ortholog l1-p3001 on 2026-11-12); that member gets no structure cross-check.
        return processed.download(key).read_text() if processed.exists(key) else None
    return None


def structure_features(seq: str, helix_mask: str, text: str | None) -> dict:
    out = {
        "af_coverage": "",
        "af_mean_plddt": "",
        "af_crossings": "",
        "af_n_face": "",
        "af_c_face": "",
        "af_kr_face_plus": "",
        "af_kr_face_minus": "",
        "af_n_side": "",
        "af_note": "",
    }
    if text is None:
        out["af_note"] = "no chosen model"
        return out
    mseq, xyz, plddt = read_ca(text)
    # Member residue i -> model residue i - off. Exact substring either way.
    if (k := seq.find(mseq)) >= 0:
        off = k
    elif (k := mseq.find(seq)) >= 0:
        off = -k
    else:
        out["af_note"] = "model sequence does not align exactly to the member"
        return out
    idx = np.array([i - off for i in range(len(seq))])
    ok = (idx >= 0) & (idx < len(mseq))
    out["af_coverage"] = round(float(ok.mean()), 3)
    z = np.full(len(seq), np.nan)
    p = np.full(len(seq), np.nan)
    pos = np.full((len(seq), 3), np.nan)
    pos[ok] = xyz[idx[ok]]
    p[ok] = plddt[idx[ok]]
    out["af_mean_plddt"] = round(float(np.nanmean(p)), 1)

    # Helix axis from v1 + v2 helix positions (positions only).
    axes, pts = [], []
    for m in re.finditer("H+", helix_mask):
        sel = [i for i in range(m.start(), m.end()) if ok[i] and p[i] >= FIT_PLDDT]
        if len(sel) < MIN_HELIX_CA:
            continue
        a = pos[sel[-4:]].mean(0) - pos[sel[:4]].mean(0)
        a /= np.linalg.norm(a)
        if axes and np.dot(a, axes[0]) < 0:
            a = -a
        axes.append(a)
        pts.append(pos[sel])
    if len(axes) < 2:
        out["af_note"] = "fewer than two confident helices to fit a membrane"
        return out
    n = np.mean(axes, axis=0)
    # One refit without helices more than 45 degrees off the mean (amphipathic or kinked).
    keep = [i for i, a in enumerate(axes) if np.dot(a, n / np.linalg.norm(n)) > np.cos(np.pi / 4)]
    if len(keep) >= 2:
        n = np.mean([axes[i] for i in keep], axis=0)
        pts = [pts[i] for i in keep]
    n /= np.linalg.norm(n)
    center = np.vstack(pts).mean(0)
    z[ok] = (pos[ok] - center) @ n

    conf = ok & (p >= WALK_PLDDT)
    state, crossings = 0, 0
    for i in np.where(conf)[0]:
        s = 1 if z[i] > SIDE_Z else (-1 if z[i] < -SIDE_Z else 0)
        if s and state and s != state:
            crossings += 1
        if s:
            state = s
    out["af_crossings"] = crossings

    first_h = helix_mask.find("H")
    last_h = helix_mask.rfind("H")
    n_res = [i for i in range(max(first_h, 1)) if ok[i]] or [i for i in range(3) if ok[i]]
    c_res = [i for i in range(last_h + 1, len(seq)) if ok[i]] or [
        i for i in range(len(seq) - 3, len(seq)) if ok[i]
    ]
    nz = float(np.nanmean(z[n_res])) if n_res else float("nan")
    cz = float(np.nanmean(z[c_res])) if c_res else float("nan")
    out["af_n_face"] = "" if np.isnan(nz) else ("+" if nz > 0 else "-")
    out["af_c_face"] = "" if np.isnan(cz) else ("+" if cz > 0 else "-")

    core = np.where(ok & (np.abs(np.nan_to_num(z, nan=99)) < CORE_Z))[0]
    near = np.zeros(len(seq), bool)
    for i in core:
        near[max(0, i - FLANK_RES) : i + FLANK_RES + 1] = True
    kr = np.array([c in "KR" for c in seq])
    band = ok & near & kr & conf
    zz = np.nan_to_num(z, nan=0.0)
    plus = int((band & (zz > FLANK_Z[0]) & (zz < FLANK_Z[1])).sum())
    minus = int((band & (zz < -FLANK_Z[0]) & (zz > -FLANK_Z[1])).sum())
    out["af_kr_face_plus"], out["af_kr_face_minus"] = plus, minus
    if abs(plus - minus) >= KR_MARGIN and out["af_n_face"]:
        cyto = "+" if plus > minus else "-"
        out["af_n_side"] = "in" if out["af_n_face"] == cyto else "out"
    return out


# --- flags ------------------------------------------------------------------------------

#: GPCR classes whose members normally carry a cleavable signal peptide. A leaderless member
#: of one suggests a truncated or alternative start. Tedman class letters as in the screen.
SP_CLASSES = {"B", "C", "F", "Adhesion"}


def flags(r: dict) -> list[str]:
    f = []
    # Condition 1: extracellular N-terminus.
    if r["up_n_side"] == "in" and not r["up_signal_end"]:
        f.append(f"C1: UniProt canonical N-terminus cytoplasmic ({r['up_n_evidence']})")
    if r["v2_n_side"] == "in":
        f.append("C1: DeepTMHMM v2 N-terminus cytoplasmic")
    if r["af_n_side"] == "in":
        f.append(
            "C1: model + positive-inside put the N-terminus on the cytoplasmic face "
            f"(K/R {r['af_kr_face_plus']} vs {r['af_kr_face_minus']})"
        )
    # Condition 2: no cleavable signal peptide.
    if r["up_same_start"] and r["up_signal_end"]:
        f.append(f"C2: UniProt SIGNAL 1..{r['up_signal_end']}")
    if r["v2_signal_peptide_aa"]:
        f.append(f"C2: DeepTMHMM v2 signal peptide 1..{r['v2_signal_peptide_aa']}")
    if r["signalp6_sp_prob"] != "" and float(r["signalp6_sp_prob"]) >= SP_PROB:
        f.append(f"C2: SignalP 6 SP p={float(r['signalp6_sp_prob']):.2f}")
    # Condition 3: at least two helices.
    if r["v2_tm"] != "" and int(r["v2_tm"]) < 2:
        f.append(f"C3: DeepTMHMM v2 calls {r['v2_tm']} helix(es)")
    if r["af_crossings"] != "" and int(r["af_crossings"]) < 2 and r["af_coverage"] == 1.0:
        f.append(f"C3: model crosses the membrane {r['af_crossings']} time(s)")
    if r["up_tm"] != "" and int(r["up_tm"]) < 2:
        f.append(f"C3: UniProt canonical has {r['up_tm']} TM")
    # Start of the ordered protein.
    if not r["sequence"].startswith("M"):
        f.append("start: no Met at residue 1")
    tm1 = r["up_canonical_tm1_start"]
    start = r["up_start_in_canonical"]
    if start != "" and tm1 != "" and int(start) >= int(tm1):
        f.append(f"start: begins at canonical residue {start}, at or after canonical TM1 ({tm1})")
    if r["up_signal_end"] and not r["up_same_start"]:
        f.append(
            f"start: canonical carries SIGNAL 1..{r['up_signal_end']}, but the member starts elsewhere"
        )
    if r["gpcr"] and r["tedman_class"] in SP_CLASSES:
        f.append(f"family: class {r['tedman_class']} GPCR without a signal peptide")
    # C-terminal side: v1's helix count fixes it; GPCRs end inside.
    v1_c = "out" if r["n_tm"] % 2 == 0 else "in"
    if r["gpcr"] and v1_c == "out":
        f.append(f"family: v1 calls {r['n_tm']} helices, so the GPCR C-terminus lands outside")
    if r["v2_tm"] != "" and r["v2_n_side"] == "out" and (int(r["v2_tm"]) - r["n_tm"]) % 2:
        f.append(f"count: v2 {r['v2_tm']} helices vs v1 {r['n_tm']}, C-terminal side differs")
    if (
        r["af_crossings"] != ""
        and r["af_coverage"] == 1.0
        and (int(r["af_crossings"]) - r["n_tm"]) % 2
    ):
        f.append(
            f"count: model {r['af_crossings']} crossings vs v1 {r['n_tm']} helices, "
            "C-terminal side differs"
        )
    return f


def support(r: dict) -> list[str]:
    """Independent sources that affirm all three conditions for this member."""
    s = []
    if r["up_n_side"] == "out" and not r["up_signal_end"]:
        s.append(f"UniProt N-out ({r['up_n_evidence']})")
    if r["v2_n_side"] == "out" and not r["v2_signal_peptide_aa"] and int(r["v2_tm"] or 0) >= 2:
        s.append(f"v2 N-out, no SP, {r['v2_tm']} TM")
    if r["signalp6_prediction"] == "OTHER":
        s.append(f"SignalP 6 no SP (p_SP={float(r['signalp6_sp_prob']):.3f})")
    if r["af_n_side"] == "out" and int(r["af_crossings"] or 0) >= 2:
        s.append(
            f"model N-out, {r['af_crossings']} crossings, pLDDT {r['af_mean_plddt']} "
            f"(K/R {r['af_kr_face_plus']}/{r['af_kr_face_minus']})"
        )
    if r["in_tedman"] and r["tedman_au"] != "":
        s.append(f"Tedman HA-after-Met1 {float(r['tedman_au']):,.0f} a.u.")
    return s


def screen() -> pd.DataFrame:
    tp = pd.read_csv(TAG_PATHS, sep="\t", dtype=str)
    rung = tp[tp["placement"] == RUNG].copy()
    rung["n_tm"] = rung["n_tm"].astype(int)
    feats = pd.read_csv(FEATURES, sep="\t", dtype=str).fillna("")
    if set(feats["candidate_id"]) != set(rung["candidate_id"]):
        raise SystemExit(
            f"{FEATURES.name} does not cover exactly the rung's {len(rung)} members; "
            "rerun nterm_res1_multipass.py first"
        )
    feats = feats.set_index("candidate_id")
    members = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    seqs = {cid: members[cid]["sequence"].rstrip("*") for cid in rung["candidate_id"]}
    for cid, sha in zip(rung["candidate_id"], rung["sequence_sha256"], strict=True):
        if topology.sequence_sha256(seqs[cid]) != sha:
            raise SystemExit(f"{cid}: designed sequence no longer matches l1_tag_paths.tsv")

    proj, helix_letter = v2_alphabet()
    v1 = ld.topology_report({c: members[c] for c in seqs}, tool=topology.TOOL_V1)
    v2 = ld.topology_report({c: members[c] for c in seqs}, tool=topology.TOOL_V2)
    v1 = v1.set_index("candidate_id")
    v2 = v2.set_index("candidate_id")
    up = load_uniprot().set_index("accession").to_dict("index")
    sp6 = (
        # Keyed on design_id: candidate ids are renumbered between library passes.
        pd.read_csv(SIGNALP, sep="\t", dtype=str).set_index("design_id").to_dict("index")
        if SIGNALP.exists()
        else {}
    )
    dids = design_ids()
    st = pd.read_csv(STRUCTURES, sep="\t", dtype=str).fillna("")
    st = st[st["is_chosen"] == "1"].set_index("candidate_id").to_dict("index")

    rows = []
    for r in rung.to_dict("records"):
        cid = r["candidate_id"]
        seq = seqs[cid]
        fe = feats.loc[cid]
        v1s = v1.at[cid, "per_residue_topology"]
        row = {
            "design_id": dids[cid],
            "candidate_id": cid,
            "gene_symbol": r["gene_symbol"],
            "candidate_set": r["candidate_set"],
            "species": r["species"],
            "gpcr": int(fe["gpcr"]),
            "tedman_class": fe["tedman_class"],
            "in_tedman": int(fe["in_tedman"]),
            "tedman_au": fe["tedman_au"],
            "length_aa": len(seq),
            "deeptmhmm_v1_lane": r["placement"],
            "deeptmhmm_v1_n_side": "out",
            "deeptmhmm_v1_signal_peptide_aa": int(r["signal_peptide_aa"]),
            "n_tm": r["n_tm"],
            "v1_n_tail_aa": int(fe["n_tail_aa"]),
            "first_residues": seq[:12],
            "uniprot_accession": fe["uniprot_accession"],
        }
        if len(v1s) != len(seq):
            raise SystemExit(f"{cid}: v1 string {len(v1s)} aa, sequence {len(seq)} aa")
        mask_v1 = "".join("H" if c == "M" else "-" for c in v1s)
        v2s = v2.at[cid, "dtm2_topology_string"] if bool(v2.at[cid, "found"]) else None
        # A v2 string whose length differs from the sequence it was matched to cannot be
        # read position by position; treat it as missing rather than misalign it.
        row["v2_note"] = "" if v2s is None or len(v2s) == len(seq) else "v2 string length differs"
        if v2s is not None and len(v2s) == len(seq):
            v2f = v2_features(v2.at[cid, "dtm2_topology_string"], proj, helix_letter)
        else:
            v2f = {
                "v2_signal_peptide_aa": "",
                "v2_n_side": "",
                "v2_tm": "",
                "v2_c_side": "",
                "_v2_helix_mask": "-" * len(seq),
            }
        mask = "".join(
            "H" if a == "H" or b == "H" else "-"
            for a, b in zip(mask_v1, v2f.pop("_v2_helix_mask"), strict=True)
        )
        row.update(v2f)
        s = sp6.get(dids[cid], {})
        if s and s["prefix70_sha256"] != topology.sequence_sha256(seq[:SIGNALP_WINDOW]):
            s = {}  # the member's first 70 residues changed since the SignalP run
        row["signalp6_prediction"] = s.get("signalp6_prediction", "")
        row["signalp6_sp_prob"] = s.get("signalp6_sp_prob", "")
        row["signalp6_source"] = s.get("signalp6_source", "")
        row["signalp6_version"] = s.get("signalp6_version", "")
        row.update(uniprot_features(seq, up.get(fe["uniprot_accession"]), helix_letter))
        srow = st.get(cid)
        row["model_source"] = srow["source"] if srow else ""
        row["model_ref"] = (srow["model_url"] or srow["r2_key"]) if srow else ""
        text = model_text(srow) if srow else None
        row["model_version"] = model_version(srow) if srow else ""
        row["model_sha256"] = hashlib.sha256(text.encode()).hexdigest() if text else ""
        row.update(structure_features(seq, mask, text))
        if srow and text is None:
            row["af_note"] = "chosen model not published in S3"
        row["sequence"] = seq
        fl = flags(row)
        row["screen_flags"] = "; ".join(fl)
        row["screen_support"] = "; ".join(support(row))
        row["flagged"] = int(bool(fl))
        row["topology_snapshot_key"] = topology.SNAPSHOT_KEY
        row["deeptmhmm_v1_topology_version"] = v1.at[cid, "topology_version"]
        row["deeptmhmm_v1_tool_version"] = v1.at[cid, "tool_version"]
        row["deeptmhmm_v2_topology_version"] = v2.at[cid, "topology_version"] if v2s else ""
        row["deeptmhmm_v2_tool_version"] = v2.at[cid, "tool_version"] if v2s else ""
        if row["deeptmhmm_v1_topology_version"] != r["topology_version"]:
            raise SystemExit(f"{cid}: v1 topology_version differs from l1_tag_paths.tsv")
        row["sequence_sha256"] = r["sequence_sha256"]
        rows.append(row)
    out = pd.DataFrame(rows).drop(columns=["sequence"]).rename(columns=SCREEN_RENAME)
    out.to_csv(SCREEN_OUT, sep="\t", index=False, lineterminator="\n")
    return out


#: Output names: v1 columns say deeptmhmm_v1; cross-check columns say which source they are.
SCREEN_RENAME = {
    "n_tm": "deeptmhmm_v1_tm",
    "v1_n_tail_aa": "deeptmhmm_v1_n_tail_aa",
    "v2_signal_peptide_aa": "deeptmhmm_v2_signal_peptide_aa",
    "v2_n_side": "deeptmhmm_v2_n_side",
    "v2_tm": "deeptmhmm_v2_tm",
    "v2_c_side": "deeptmhmm_v2_c_side",
    "v2_note": "deeptmhmm_v2_note",
    "screen_flags": "crosscheck_disagrees",
    "screen_support": "crosscheck_agrees",
    "flagged": "sent_to_review",
}

#: Fingerprints carried into every row of the verdict table, so each row says what it rests on.
FINGERPRINTS = [
    "uniprot_accession",
    "uniprot_release",
    "uniprot_entry_version",
    "uniprot_sequence_version",
    "uniprot_sequence_sha256",
    "topology_snapshot_key",
    "deeptmhmm_v1_topology_version",
    "deeptmhmm_v1_tool_version",
    "deeptmhmm_v2_topology_version",
    "deeptmhmm_v2_tool_version",
    "signalp6_version",
    "signalp6_source",
    "model_version",
    "model_ref",
    "model_sha256",
    "sequence_sha256",
]


# --- verdicts ---------------------------------------------------------------------------

VERDICT_COLS = [
    "design_id",
    "candidate_id",
    "gene_symbol",
    "candidate_set",
    "species",
    "deeptmhmm_v1_lane",
    "deeptmhmm_v1_tm",
    "crosscheck_disagrees",
    "sent_to_review",
    "lane_ok",
    "recommended_lane",
    "confidence",
    "evidence",
    "pmids",
    "decided_by",
    "implemented_by",
    "reason",
    *FINGERPRINTS,
]


def verdicts(scr: pd.DataFrame) -> pd.DataFrame:
    """v1's call stands unless a cross-check disagrees; then the review's verdict applies.

    A member with a review row takes it even when no cross-check disagrees any more (the
    review read the literature for it, and that evidence should travel with the row).
    """
    review = (
        pd.read_csv(REVIEW, sep="\t", dtype=str).fillna("").set_index("design_id")
        if REVIEW.exists()
        else pd.DataFrame()
    )
    rows = []
    for r in scr.fillna("").to_dict("records"):
        cid = r["candidate_id"]
        base = {
            k: r[k]
            for k in (
                "design_id",
                "candidate_id",
                "gene_symbol",
                "candidate_set",
                "species",
                "deeptmhmm_v1_lane",
                "deeptmhmm_v1_tm",
                "crosscheck_disagrees",
                "sent_to_review",
                *FINGERPRINTS,
            )
        }
        did = r["design_id"]
        if did not in review.index:
            if r["sent_to_review"]:
                raise SystemExit(f"{cid}: a cross-check disagrees but {REVIEW.name} has no row")
            agrees = r["crosscheck_agrees"]
            n = len(agrees.split("; ")) if agrees else 0
            exp = "experimental" in agrees or "Tedman" in agrees
            rows.append(
                {
                    **base,
                    "lane_ok": "yes",
                    "recommended_lane": RUNG,
                    "confidence": "high" if n >= 3 and exp else ("medium" if n >= 2 else "low"),
                    "evidence": "DeepTMHMM v1 call of record; cross-checks agreeing: "
                    + (agrees or "none available"),
                    "pmids": "",
                    "decided_by": DECIDED_BY_SCREEN,
                    "implemented_by": IMPLEMENTED_BY,
                    "reason": "no cross-check disagrees with the v1 call on any of the three conditions",
                }
            )
            continue
        v = review.loc[did]
        if v["reviewed_sequence_sha256"] != r["sequence_sha256"]:
            raise SystemExit(
                f"{did}: ordered sequence changed since its review; review it again in {REVIEW.name}"
            )
        if v["candidate_id"] != cid:
            raise SystemExit(
                f"{did}: {REVIEW.name} names it {v['candidate_id']}, the manifest {cid}; "
                "run --rekey-review"
            )
        rows.append(
            {
                **base,
                "lane_ok": v["lane_ok"],
                "recommended_lane": v["recommended_lane"],
                "confidence": v["confidence"],
                "evidence": f"cross-checks disagreeing: {r['crosscheck_disagrees'] or 'none'}"
                f" | review: {v['evidence']}",
                "pmids": v["pmids"],
                "decided_by": v["decided_by"],
                "implemented_by": v["implemented_by"],
                "reason": v["reason"],
            }
        )
    out = pd.DataFrame(rows, columns=VERDICT_COLS)
    extra = set(review.index) - set(scr.loc[scr["sent_to_review"] == 1, "design_id"])
    if extra:
        print(
            f"note: {len(extra)} reviewed members no cross-check now disagrees with: {sorted(extra)}"
        )
    return out


def rekey_review() -> None:
    """Bring the review table's candidate_id up to the manifest's numbering, by design_id.

    The review is keyed by design_id, which survives renumbering; candidate_id is carried for
    reading. A design_id no longer in the manifest is marked ``retired`` and kept, so the
    record of what was reviewed is not lost.
    """
    rev = pd.read_csv(REVIEW, sep="\t", dtype=str).fillna("")
    man = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
    # A retired member stays in the manifest, under its design_id, with a retired_reason.
    live = man[man["retired_reason"] == ""]
    now = dict(zip(live["design_id"], live["candidate_id"], strict=True))
    reason = dict(zip(man["design_id"], man["retired_reason"], strict=True))
    old = rev["candidate_id"].copy()
    rev["candidate_id"] = rev["design_id"].map(lambda d: now.get(d, "retired"))
    rev["retired_reason"] = rev["design_id"].map(lambda d: "" if d in now else reason.get(d, ""))
    changed = int((old != rev["candidate_id"]).sum())
    retired = sorted(rev.loc[rev["candidate_id"] == "retired", "design_id"])
    rev.to_csv(REVIEW, sep="\t", index=False, lineterminator="\n")
    print(f"{REVIEW.name}: {changed} candidate_id(s) rekeyed; retired: {retired}")


def check_drift() -> int:
    """Members whose pinned inputs have moved since `nterm_res1_lane_check.tsv` was written.

    Re-fetches only UniProt's per-entry versions (not the features), and recomputes each
    member's sequence hash and v1 topology version from the manifest and the pinned snapshot.
    """
    if not OUT.exists():
        raise SystemExit(f"{OUT.name} missing; run the check first")
    old = pd.read_csv(OUT, sep="\t", dtype=str).fillna("").set_index("candidate_id")
    changes: list[str] = []

    tp = pd.read_csv(TAG_PATHS, sep="\t", dtype=str)
    rung = tp[tp["placement"] == RUNG].set_index("candidate_id")
    for cid in sorted(set(rung.index) ^ set(old.index)):
        changes.append(f"{cid}: {'joined' if cid in rung.index else 'left'} the rung")

    members = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    here = {c: members[c] for c in old.index if c in members}
    v1 = ld.topology_report(here, tool=topology.TOOL_V1).set_index("candidate_id")
    for cid, r in old.iterrows():
        if cid not in members:
            changes.append(f"{cid}: no longer in the manifest")
            continue
        sha = topology.sequence_sha256(members[cid]["sequence"].rstrip("*"))
        if sha != r["sequence_sha256"]:
            changes.append(f"{cid}: ordered sequence hash changed")
        tv = v1.at[cid, "topology_version"] if bool(v1.at[cid, "found"]) else ""
        if tv != r["deeptmhmm_v1_topology_version"]:
            changes.append(
                f"{cid}: v1 topology_version {r['deeptmhmm_v1_topology_version']} -> {tv}"
            )
    if topology.SNAPSHOT_KEY != old["topology_snapshot_key"].iloc[0]:
        changes.append(f"topology snapshot key changed to {topology.SNAPSHOT_KEY}")

    accs = sorted({a for a in old["uniprot_accession"] if a})
    now = fetch_uniprot(accs).set_index("accession")
    release = now["uniprot_release"].iloc[0]
    pinned = sorted(set(old["uniprot_release"]) - {""})
    print(f"UniProt release now {release}; pinned {', '.join(pinned)}")
    for cid, r in old.iterrows():
        a = r["uniprot_accession"]
        if not a:
            continue
        if a not in now.index:
            changes.append(f"{cid}: UniProt {a} no longer returned")
            continue
        n = now.loc[a]
        for col, new in (
            ("uniprot_entry_version", n["entry_version"]),
            ("uniprot_sequence_version", n["sequence_version"]),
            ("uniprot_sequence_sha256", n["sequence_sha256"]),
        ):
            if r[col] != new:
                changes.append(f"{cid}: {a} {col} {r[col][:12]} -> {new[:12]}")

    for c in changes:
        print("  " + c)
    print(
        f"{len(changes)} change(s)" + ("; rerun these members" if changes else "; nothing to rerun")
    )
    return 1 if changes else 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--signalp-dir", type=Path, help="rebuild the SignalP 6 snapshot first")
    ap.add_argument("--screen-only", action="store_true", help="stop after the screen table")
    ap.add_argument("--refresh-uniprot", action="store_true", help="re-fetch the UniProt snapshot")
    ap.add_argument(
        "--check-drift",
        action="store_true",
        help="re-fetch fingerprints only; list members to rerun; exit 1 if any changed",
    )
    ap.add_argument(
        "--rekey-review",
        action="store_true",
        help="rewrite the review table's candidate_id from the manifest, by design_id",
    )
    args = ap.parse_args(argv)
    if args.rekey_review:
        rekey_review()
    if args.check_drift:
        return check_drift()
    if args.refresh_uniprot:
        df = fetch_uniprot(rung_accessions())
        df.to_csv(UNIPROT_SNAPSHOT, sep="\t", index=False, lineterminator="\n")
        print(
            f"wrote {UNIPROT_SNAPSHOT.name}: {len(df)} entries, release {df['uniprot_release'].iloc[0]}"
        )
    if args.signalp_dir:
        tp = pd.read_csv(TAG_PATHS, sep="\t", dtype=str)
        cids = set(tp.loc[tp["placement"] == RUNG, "candidate_id"])
        m = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
        build_signalp_snapshot(
            args.signalp_dir, {c: m[c]["sequence"].rstrip("*") for c in sorted(cids)}
        )

    scr = screen()
    n_rev = int(scr["sent_to_review"].sum())
    print(f"\nwrote {SCREEN_OUT.name}: {len(scr)} members, {n_rev} with a cross-check disagreeing")
    cov = {
        "v2": int((scr["deeptmhmm_v2_tm"].astype(str) != "").sum()),
        "SignalP 6": int((scr["signalp6_prediction"].fillna("") != "").sum()),
        "UniProt same start": int(scr["up_same_start"].sum()),
        "model membrane fit": int((scr["af_crossings"].astype(str) != "").sum()),
        "model side call": int((scr["af_n_side"].fillna("") != "").sum()),
    }
    print("coverage:", cov)
    kinds = collections.Counter(
        f.split(":")[0] + ":" + re.sub(r"[\d.,=()/]+", "#", f.split(":", 1)[1]).strip()
        for fl in scr["crosscheck_disagrees"].fillna("")
        for f in fl.split("; ")
        if f
    )
    print("flag kinds:")
    for k, n in kinds.most_common():
        print(f"  {n:>4}  {k}")
    print(
        pd.crosstab(scr["sent_to_review"], scr["gpcr"].map({1: "GPCR", 0: "non-GPCR"})).to_string()
    )
    if args.screen_only:
        return 0

    out = verdicts(scr)
    out.to_csv(OUT, sep="\t", index=False, lineterminator="\n")
    print(f"\nwrote {OUT.name}: {len(out)} rows")
    print(out.groupby(["lane_ok", "confidence"]).size().to_string())
    moves = out[out["recommended_lane"] != RUNG]
    print(f"\nrecommended for another lane ({len(moves)}):")
    for r in moves.itertuples():
        print(
            f"  {r.design_id} {r.candidate_id} {r.gene_symbol:10s} -> {r.recommended_lane} ({r.confidence})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
