"""Biological-correspondence review of the 200-member CDS-eligible ortholog review panel.

For every selected mouse/cyno form, compare the exact selected protein (translated from
the archived selected CDS) with the exact human L1 library protein (translated designed
CDS), both matched by sequence hash, never accession. Evidence per pair:

* global alignment (BLOSUM62, BLAST gap costs 11/1, free end gaps): coverage, identity,
  human-only and ortholog-only segments classed N-terminal / internal / C-terminal;
* DeepTMHMM v1 (``topology.TOOL_V1``) by exact sequence from the pinned 2026-11-06
  snapshot: label, signal peptide, TM helices (each human helix mapped through the
  alignment), N-terminal orientation;
* GPI: the human class comes from the frozen GPI set, the ortholog call from NetGPI 1.1
  on the exact ortholog sequence (``anchor_call`` recorded in the panel);
* Pfam 39.0 domains on both proteins with each family's curated gathering threshold
  (hmmsearch via pyhmmer). A human domain is conserved when the ortholog has a hit of
  the same Pfam accession that overlaps the residues aligned to it.

No new numeric cutoff is introduced: every flag is a categorical presence/absence call.
Recommendations follow fixed rules (see ``recommend``). For strata holding a form
recommended for replacement, the archived seeded draw is replayed exactly and the next
unused candidates in shuffle order are reviewed with the same evidence, so the proposed
replacement is the first architecture-conserving candidate after the drawn ones.

Manual overrides (attribution per Becca Carlson's convention of 2026-11-10: "for any
manual fixes please note the agent who made the fix / adjustment and reasoning"): the one
reviewer departure from the rule, mouse OR10G6 replace -> ask, was an agent judgement made
on 2026-11-07 and implemented by the Claude Code session ortholog-correspondence-review
(started by l1-reviewer-package), PR 273. Becca kept all 6 asks as drawn on 2026-11-10
(PR 283). Its reason is in ``ortholog_correspondence_judgements.tsv``. Every output row
carries ``decided_by`` / ``implemented_by`` / ``recommendation_basis``; rows decided by the
mechanical rule say so. The override is applied by code from that table, so all outputs
regenerate.

Review only: no manifest, construct table, pipeline, D1 record or existing S3 key changes.

    uv run --locked --group modal --group domains python \
        analysis/sprint4/s4e9-orf-design-fg-panel/review_ortholog_correspondence.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from deliverome_analysis import s3, seqtools, topology  # noqa: E402

RULE_VERSION = "ortholog-correspondence-review-v1"
FINAL = (
    "sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-cds-eligibility/2026-11-05/final-v2/"
)
REBUILD = "sprint4/s4e9-orf-design-fg-panel/library-design/ortholog-rebuild-inputs/2026-11-05/"

#: name -> (bucket, key, sha256 of stored object, sha256 after gunzip or None)
INPUTS = {
    "selected_cds": (
        "processed",
        FINAL + "results/panel/selected_cds.json.gz",
        "4f83238205b7f4fa1156d0129d40050300a794cfb9881ace35e1102a2162b222",
        "ccfc8f2d16f8e03e264428159323c1c801c1f2ba5fd7e019d2a483578a5d2a37",
    ),
    "human_cds": (
        "processed",
        "sprint4/s4e9-orf-design-fg-panel/library-candidates/candidate_cds.json",
        "03c3413d70e15612868133d15839e092dba20fb05f3e8318d17f48928c91c662",
        None,
    ),
    "topology_review": (
        "external",
        "accessible-surfaceome/topology/2026-11-06-ortholog-cds-policy/topology.parquet",
        "7683d75a3e04228d95fcfb8fe9b7e9b46b38c0ea254ca42550a666f5d8bbf216",
        None,
    ),
    "pfam": (
        "external",
        "pfam/39.0/Pfam-A.hmm.gz",
        "277356d1ceeba7ed3fc554bfe0ad58fdd87e52c26c1cf490cc90610dbc098f0e",
        "4f86f23f7d673cf64e35a786dd89c0159fb6cafafd5301615162af9927b3c243",
    ),
    "pfam_clans": (
        "external",
        "pfam/39.0/Pfam-A.clans.tsv.gz",
        "69b2cfad70b4b0e201826b32d97a27bffa6f1cdf7b33d951a9e3a98fd020be20",
        None,
    ),
    # Exact inputs of the archived draw (workflow stage 3), for the replacement queue.
    "draw_constructs": (
        "processed",
        FINAL + "inputs/candidate_constructs.json.gz",
        "11d645d255037f0cb2ab855fc438df16178f04c5d90cee6fd56e22cc176e1358",
        "cdc5fa5e119e7590d910b44db03138a829f7a59e6d8ec430d8ebbe9f21008e3f",
    ),
    "draw_manifest": (
        "processed",
        REBUILD + "candidate_manifest_tranches.tsv",
        "7ad0538b831b521a10526f290f3b40eb9ec1f2333f0d674f88a831092c13cdf4",
        None,
    ),
    "draw_gpi": (
        "processed",
        REBUILD + "gpi_set.tsv",
        "dbb1db6e68bc3e1bdb3e69fee365e963a5a36f90434f88dd9e8d42a0df6f1b20",
        None,
    ),
    "draw_topology": (
        "external",
        "accessible-surfaceome/topology/2026-11-05-l1-compara116/topology.parquet",
        "604521dcbc920f520ea566aee3da638dad6540772c0f38d1792a631ca301e048",
        None,
    ),
    "draw_eligible": (
        "processed",
        FINAL + "results/eligibility/eligible_proteins.json.gz",
        "2c301da225cd910aa74ae1a2aefe42815e4978cea98ca3d60da7c489461ee0d9",
        "1cabf363947d346cda1bae8b8863477bd0c9f5197cbe3efe6a1ea7f93d7fa76d",
    ),
    "draw_picks": (
        "processed",
        FINAL + "results/scored/ortholog_provisional_picks.tsv.gz",
        "760b3a7d3acd4525df05d88c64849502747a03e5d8ca7538cdab79a484717671",
        "58def62b49dde1973a6721a20957e88371521a747da47f00d07450d48d2322f5",
    ),
    "draw_score_summary": (
        "processed",
        FINAL + "results/scored/ortholog_scoring_summary.json",
        "4803d67379fd992a185a9b2bde2fae0ce223ceb6df58798520c4980b97c2a8b5",
        None,
    ),
    "draw_proteins": (
        "processed",
        FINAL + "results/scored/ortholog_picks.fasta",
        "4e524e85b8be39ed170dbf5c20d32de148ff722639392d83c297d64a03f971ce",
        None,
    ),
}
PANEL = HERE / "ortholog_panel_cds_eligible.tsv"
#: Reviewer departures from the mechanical rule: form_id, rule_recommendation,
#: recommendation, decided_by, implemented_by, reason. Each row must name the rule result it
#: overrides and who decided and implemented it (Becca's attribution convention, 2026-11-10).
JUDGEMENTS = HERE / "ortholog_correspondence_judgements.tsv"
JUDGEMENT_COLUMNS = (
    "form_id",
    "rule_recommendation",
    "recommendation",
    "decided_by",
    "implemented_by",
    "reason",
)
RULE_DECIDED_BY = "no decision taken (mechanical rule)"
RULE_IMPLEMENTED_BY = (
    "review_ortholog_correspondence.py recommend(), Claude Code session "
    "ortholog-correspondence-review, PR 273"
)
OLD_PANEL = HERE / "ortholog_panel_compara116.tsv"
OLD_SUMMARY = HERE / "ortholog_panel_compara116_summary.json"
NAMED_CASES = {("cyno", "SLC8A1"), ("cyno", "IL31RA"), ("cyno", "MS4A10")}
ALTERNATES_PER_STRATUM = 6

#: Flags that mean the human architecture is not reproduced. All categorical.
ARCHITECTURE_FLAGS = (
    "signal_peptide_lost",
    "signal_peptide_gained",
    "tm_helix_lost",
    "tm_helix_gained",
    "topology_label_changed",
    "n_terminal_orientation_changed",
    "gpi_anchor_lost",
    "pfam_domain_lost",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(name: str, cache: Path) -> Path:
    """Download one pinned input, verify its checksum, gunzip when a decompressed hash is set."""
    bucket, key, digest, inner = INPUTS[name]
    stored = cache / key
    if not stored.exists() or sha256_file(stored) != digest:
        stored.parent.mkdir(parents=True, exist_ok=True)
        getattr(s3, bucket).download(key, stored)
    if sha256_file(stored) != digest:
        raise ValueError(f"{bucket}://{key} checksum mismatch")
    if inner is None:
        return stored
    plain = stored.with_suffix("")
    if not plain.exists() or sha256_file(plain) != inner:
        with gzip.open(stored, "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    if sha256_file(plain) != inner:
        raise ValueError(f"{bucket}://{key} decompressed checksum mismatch")
    return plain


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def runs(mask: list[bool]) -> list[tuple[int, int]]:
    """Half-open [start, end) runs where mask is True."""
    out, start = [], None
    for i, value in enumerate([*mask, False]):
        if value and start is None:
            start = i
        elif not value and start is not None:
            out.append((start, i))
            start = None
    return out


def char_runs(text: str, char: str) -> list[tuple[int, int]]:
    return runs([c == char for c in text])


# --- evidence -------------------------------------------------------------------------


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.open_end_gap_score = 0
    aligner.extend_end_gap_score = 0
    return aligner


def align_map(aligner: PairwiseAligner, human: str, ortholog: str) -> list[int | None]:
    """human residue index -> aligned ortholog residue index (or None)."""
    alignment = aligner.align(human, ortholog)[0]
    mapping: list[int | None] = [None] * len(human)
    for (h0, h1), (o0, _o1) in zip(*alignment.aligned, strict=True):
        for k in range(h1 - h0):
            mapping[h0 + k] = o0 + k
    return mapping


def segments(length: int, matched: set[int]) -> list[dict]:
    out = []
    for start, end in runs([i not in matched for i in range(length)]):
        where = "N" if start == 0 else "C" if end == length else "internal"
        out.append({"start": start + 1, "end": end, "length": end - start, "where": where})
    return out


def text(value: str | bytes) -> str:
    return value.decode() if isinstance(value, bytes) else value


def load_clans(path: Path) -> dict[str, str]:
    """Pfam family accession -> clan accession (families outside a clan are omitted)."""
    out = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            family, clan, *_ = line.rstrip("\n").split("\t")
            if clan:
                out[family] = clan
    return out


def pfam_scan(
    sequences: dict[str, str], hmm_path: Path, clans: dict[str, str]
) -> dict[str, list[dict]]:
    """Pfam domains per key, down to each family's noise cutoff (NC), tiered against GA.

    ``tier`` is ``gathering`` when both the sequence and domain scores reach the family's
    gathering threshold (GA, the Pfam annotation cutoff) and ``noise`` when they reach only
    its noise cutoff (NC, the highest score Pfam curators saw for a known false positive).
    Overlapping hits within one clan (or one family outside a clan) are then resolved to
    the best-scoring hit, as Pfam's own pfam_scan.pl does, separately per tier.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    keys = sorted(sequences)
    digital = [
        pyhmmer.easel.TextSequence(name=k, sequence=sequences[k]).digitize(alphabet) for k in keys
    ]
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, digital)
    hits: dict[str, list[dict]] = defaultdict(list)
    with pyhmmer.plan7.HMMFile(hmm_path) as hmms:
        for top in pyhmmer.hmmsearch(hmms, block, bit_cutoffs="noise"):
            accession = text(top.query.accession)
            family = accession.split(".")[0]
            name = text(top.query.name)
            ga_seq, ga_dom = top.query.cutoffs.gathering
            for hit in top:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    if not domain.included:
                        continue
                    gathering = hit.score >= ga_seq and domain.score >= ga_dom
                    hits[text(hit.name)].append(
                        {
                            "pfam": family,
                            "pfam_version": accession,
                            "name": name,
                            "clan": clans.get(family, ""),
                            "start": int(domain.env_from),
                            "end": int(domain.env_to),
                            "score": round(float(domain.score), 1),
                            "tier": "gathering" if gathering else "noise",
                        }
                    )
    out = {}
    for k in keys:
        resolved = []
        for tier in ("gathering", "noise"):
            resolved += resolve_clan_overlaps([d for d in hits.get(k, []) if d["tier"] == tier])
        out[k] = sorted(resolved, key=lambda d: (d["start"], d["pfam"]))
    return out


def clan_key(d: dict) -> str:
    return d["clan"] or d["pfam"]


def resolve_clan_overlaps(hits: list[dict]) -> list[dict]:
    """Keep the best-scoring hit among overlapping hits sharing a clan (pfam_scan.pl rule)."""
    kept: list[dict] = []
    for d in sorted(hits, key=lambda d: (-d["score"], d["start"], d["pfam"])):
        if not any(
            clan_key(k) == clan_key(d) and k["start"] <= d["end"] and d["start"] <= k["end"]
            for k in kept
        ):
            kept.append(d)
    return kept


def topology_rows(sequences: dict[str, str], species: str, snapshot: pd.DataFrame) -> dict:
    report = topology.resolve(
        sequences,
        tool=topology.TOOL_V1,
        cohorts=topology.cohorts_for_species(species),
        snapshot=snapshot,
    )
    missing = report.loc[~report["found"], "key"].tolist()
    if missing:
        raise ValueError(f"No exact-sequence {topology.TOOL_V1} row for {species}: {missing}")
    return {r["key"]: r for r in report.to_dict("records")}


def one_to_one(
    human: list[tuple[int, int]],
    ortholog: list[tuple[int, int]],
    mapping: list[int | None],
    pairs: bool = False,
) -> tuple[list, set[int]]:
    """Pair human and ortholog intervals one-to-one through the alignment.

    Overlap = number of human residues in the human interval whose aligned partner lies in
    the ortholog interval. Pairs maximise total overlap (Hungarian assignment); a human
    interval counts as present only if its partner overlaps by at least one residue, so two
    human helices can never both be satisfied by one ortholog helix.
    """
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    if not human:
        return [], set()
    if not ortholog:
        return [None if pairs else False] * len(human), set()
    overlap = np.zeros((len(human), len(ortholog)), dtype=int)
    for i, (a, b) in enumerate(human):
        for k in range(a, b):
            j = mapping[k]
            if j is None:
                continue
            for n, (c, d) in enumerate(ortholog):
                if c <= j < d:
                    overlap[i, n] += 1
    rows, cols = linear_sum_assignment(-overlap)
    matched: list = [None if pairs else False] * len(human)
    used = set()
    for i, n in zip(rows, cols, strict=True):
        if overlap[i, n] > 0:
            matched[i] = int(n) if pairs else True
            used.add(int(n))
    return matched, used


def first_orientation(per_residue: str) -> str:
    return next((c for c in per_residue if c in "IO"), "")


def compare(
    human: str,
    ortholog: str,
    human_topo: dict,
    ortholog_topo: dict,
    human_pfam: list[dict],
    ortholog_pfam: list[dict],
    human_gpi: bool,
    anchor_call: str,
    aligner: PairwiseAligner,
) -> dict:
    mapping = align_map(aligner, human, ortholog)
    matched_h = {i for i, j in enumerate(mapping) if j is not None}
    matched_o = {j for j in mapping if j is not None}
    identical = sum(1 for i, j in enumerate(mapping) if j is not None and human[i] == ortholog[j])
    h_gaps = segments(len(human), matched_h)
    o_gaps = segments(len(ortholog), matched_o)

    ht, ot = human_topo["per_residue_topology"], ortholog_topo["per_residue_topology"]
    if len(ht) != len(human) or len(ot) != len(ortholog):
        raise ValueError("Topology string length differs from its sequence")
    h_sp = int(human_topo["signal_peptide_length"] or 0)
    o_sp = int(ortholog_topo["signal_peptide_length"] or 0)
    h_tm, o_tm = char_runs(ht, "M"), char_runs(ot, "M")
    tm_matched, o_tm_hit = one_to_one(h_tm, o_tm, mapping)

    features = []
    if h_sp:
        features.append(("SP", 0, h_sp))
    features += [(f"TM{n + 1}", a, b) for n, (a, b) in enumerate(h_tm)]
    human_pfam = [d for d in human_pfam if d["tier"] == "gathering"]
    lost_domains, partial_domains, weak_domains, disrupted_domains = [], [], [], []
    for key in sorted({clan_key(d) for d in human_pfam}):
        mine = [d for d in human_pfam if clan_key(d) == key]
        theirs = [e for e in ortholog_pfam if clan_key(e) == key]
        present, _ = one_to_one(
            [(d["start"] - 1, d["end"]) for d in mine],
            [(e["start"] - 1, e["end"]) for e in theirs],
            mapping,
            pairs=True,
        )
        for d, partner in zip(mine, present, strict=True):
            a, b = d["start"] - 1, d["end"]
            label = f"{d['name']}({d['pfam']}):{d['start']}-{d['end']}"
            unaligned = sum(mapping[i] is None for i in range(a, b))
            if partner is None:
                lost_domains.append(label)
                mapped = [mapping[i] for i in range(a, b) if mapping[i] is not None]
                inserted = bool(mapped) and any(
                    j not in matched_o for j in range(min(mapped), max(mapped) + 1)
                )
                if unaligned or inserted or not mapped:
                    disrupted_domains.append(label)
                continue
            e = theirs[partner]
            if e["tier"] == "noise" or e["pfam"] != d["pfam"]:
                weak_domains.append(f"{label}->{e['name']}({e['pfam']},{e['tier']})")
            if unaligned:
                partial_domains.append(f"{label}[{unaligned} unaligned]")
    features += [
        (f"{d['name']}({d['pfam']}):{d['start']}-{d['end']}", d["start"] - 1, d["end"])
        for d in human_pfam
    ]
    human_set = {clan_key(d) for d in human_pfam}
    gained_domains = sorted(
        {
            f"{e['name']}({e['pfam']})"
            for e in ortholog_pfam
            if e["tier"] == "gathering" and clan_key(e) not in human_set
        }
    )

    for seg in h_gaps:
        a, b = seg["start"] - 1, seg["end"]
        seg["features"] = [name for name, fa, fb in features if fa < b and fb > a]
    h_label, o_label = human_topo["deeptmhmm_label"], ortholog_topo["deeptmhmm_label"]
    h_orient, o_orient = first_orientation(ht), first_orientation(ot)
    flags = {
        "signal_peptide_lost": bool(h_sp) and not o_sp,
        "signal_peptide_gained": not h_sp and bool(o_sp),
        "tm_helix_lost": not all(tm_matched),
        "tm_helix_gained": len(o_tm) > len(o_tm_hit),
        "topology_label_changed": h_label != o_label,
        "n_terminal_orientation_changed": h_orient != o_orient,
        "gpi_anchor_lost": human_gpi and anchor_call != "anchored",
        "pfam_domain_lost": bool(lost_domains),
    }
    indel_loss = bool(disrupted_domains)
    in_feature = [s for s in h_gaps if s["features"]]
    return {
        "human_length": len(human),
        "ortholog_length": len(ortholog),
        "length_ratio": round(len(ortholog) / len(human), 4),
        "human_aligned_fraction": round(len(matched_h) / len(human), 4),
        "aligned_identity": round(identical / len(matched_h), 4) if matched_h else 0.0,
        "human_unmatched_n": sum(s["length"] for s in h_gaps if s["where"] == "N"),
        "human_unmatched_internal": sum(s["length"] for s in h_gaps if s["where"] == "internal"),
        "human_unmatched_c": sum(s["length"] for s in h_gaps if s["where"] == "C"),
        "ortholog_only_n": sum(s["length"] for s in o_gaps if s["where"] == "N"),
        "ortholog_only_internal": sum(s["length"] for s in o_gaps if s["where"] == "internal"),
        "ortholog_only_c": sum(s["length"] for s in o_gaps if s["where"] == "C"),
        "human_unmatched_segments": ";".join(
            f"{s['where']}:{s['start']}-{s['end']}"
            + (f"[{','.join(s['features'])}]" if s["features"] else "")
            for s in h_gaps
        ),
        "ortholog_only_segments": ";".join(f"{s['where']}:{s['start']}-{s['end']}" for s in o_gaps),
        "human_topology_label": h_label,
        "ortholog_topology_label": o_label,
        "human_signal_peptide_length": h_sp,
        "ortholog_signal_peptide_length": o_sp,
        "human_tm_helices": len(h_tm),
        "ortholog_tm_helices": len(o_tm),
        "human_tm_helices_unmatched": ";".join(
            f"TM{n + 1}:{a + 1}-{b}"
            for n, ((a, b), ok) in enumerate(zip(h_tm, tm_matched, strict=True))
            if not ok
        ),
        "human_n_orientation": h_orient,
        "ortholog_n_orientation": o_orient,
        "human_gpi": human_gpi,
        "ortholog_anchor_call": anchor_call,
        "human_pfam": ";".join(
            f"{d['name']}({d['pfam']}):{d['start']}-{d['end']}" for d in human_pfam
        ),
        "ortholog_pfam": ";".join(
            f"{d['name']}({d['pfam']}):{d['start']}-{d['end']}"
            + ("[NC]" if d["tier"] == "noise" else "")
            for d in ortholog_pfam
        ),
        "pfam_domains_lost": ";".join(lost_domains),
        "pfam_domains_lost_with_indel": ";".join(disrupted_domains),
        "pfam_domain_lost_with_indel": indel_loss,
        "pfam_domains_partially_aligned": ";".join(partial_domains),
        "pfam_domains_matched_by_clan_or_noise_tier": ";".join(weak_domains),
        "pfam_domains_gained": ";".join(gained_domains),
        "unmatched_segment_in_human_feature": bool(in_feature),
        **flags,
        "architecture_flags": ";".join(k for k in ARCHITECTURE_FLAGS if flags[k]),
        "architecture_conserved": not any(flags.values()),
    }


def recommend(row: dict) -> tuple[str, str]:
    """Mechanical rule on categorical presence/absence; no numeric cutoff.

    replace -- the ortholog form lacks a human TM helix, signal peptide or GPI anchor at
               the aligned position, or lacks a human Pfam domain whose span carries an
               indel (human residues deleted, or ortholog residues inserted, inside it).
    ask     -- nothing above is lost, but either a human Pfam domain has no same-family or
               same-clan hit at or above the noise cutoff although its span aligns without
               any indel (HMM non-detection on a diverged sequence is not evidence that the
               fold is absent), or the predicted architecture differs without a loss (an
               extra helix or signal peptide, a changed DeepTMHMM label or N-terminal side).
    keep    -- every human feature has a counterpart and label and orientation agree.

    Short indels inside a feature that is still present (``pfam_domains_partially_aligned``,
    ``unmatched_segment_in_human_feature``) are reported but do not change the call: signal
    peptides in particular vary in length between species and are cleaved.
    Reviewer judgements in JUDGEMENTS override this rule, each with its reason.
    """
    lost = [
        k
        for k in (
            "tm_helix_lost",
            "signal_peptide_lost",
            "gpi_anchor_lost",
            "pfam_domain_lost_with_indel",
        )
        if row[k]
    ]
    if lost:
        return "replace", "human feature absent: " + ",".join(lost)
    if row["pfam_domain_lost"]:
        return "ask", "Pfam domain undetected in an indel-free aligned span: " + row[
            "pfam_domains_lost"
        ]
    if row["architecture_flags"]:
        return "ask", "architecture differs without a lost feature: " + row["architecture_flags"]
    return "keep", "human architecture reproduced"


# --- seeded draw replay ---------------------------------------------------------------


def replay_draw(cache: Path) -> tuple[list[dict], dict, dict]:
    """Rerun the archived draw on its exact inputs; capture the candidates it shuffled."""
    import draw_eligible_ortholog_panel as dep

    captured = {}
    original = dep.draw

    def recording(rows, fixed, policy):
        captured.update(rows=rows, fixed=dict(fixed), policy=policy)
        return original(rows, fixed, policy)

    work = Path(tempfile.mkdtemp(prefix="ortholog-draw-replay-"))
    scored = work / "scored"
    scored.mkdir()
    shutil.copy(fetch("draw_picks", cache), scored / "ortholog_provisional_picks.tsv")
    shutil.copy(fetch("draw_score_summary", cache), scored / "ortholog_scoring_summary.json")
    shutil.copy(fetch("draw_proteins", cache), scored / "ortholog_picks.fasta")
    args = argparse.Namespace(
        old_panel=OLD_PANEL,
        old_summary=OLD_SUMMARY,
        cds=fetch("human_cds", cache),
        constructs=fetch("draw_constructs", cache),
        manifest=fetch("draw_manifest", cache),
        gpi=fetch("draw_gpi", cache),
        topology=fetch("draw_topology", cache),
        eligible=fetch("draw_eligible", cache),
        scored=scored,
        out=work / "panel",
    )
    dep.draw = recording
    try:
        dep.run(args)
    finally:
        dep.draw = original
    replayed = read_tsv(args.out / "ortholog_panel_cds_eligible.tsv")
    if replayed != read_tsv(PANEL):
        raise ValueError("Replayed draw does not reproduce the committed review panel")
    from export_ortholog_predictions import fasta

    proteins = fasta(scored / "ortholog_picks.fasta")
    shutil.rmtree(work)
    return captured["rows"], captured["fixed"], {"policy": captured["policy"], "proteins": proteins}


def stratum_queues(rows: list[dict], fixed: dict, policy: dict) -> dict:
    """Shuffle order per stratum, identical to draw_eligible_ortholog_panel.draw."""
    import random

    groups = defaultdict(list)
    for row in rows:
        groups[(row["species"], row["band"], row["topology_class"])].append(row)
    used = dict(fixed)
    queues = {}
    for species in ("mouse", "cyno"):
        rng = random.Random(policy["seed"] + (species == "cyno"))
        for band in ("near", "diverged"):
            classes = sorted(
                k.split("|")[2] for k in policy["quotas"] if k.startswith(f"{species}|{band}|")
            )
            for cls in classes:
                key = (species, band, cls)
                target = policy["quotas"]["|".join(key)]
                candidates = sorted(
                    groups[key],
                    key=lambda r: (
                        r["human_ensembl_gene"],
                        r["ortholog_ensembl_gene"],
                        r["sequence_sha256"],
                        r["picked_acc"],
                    ),
                )
                rng.shuffle(candidates)
                picked, queue = 0, []
                for position, row in enumerate(candidates, 1):
                    h = row["sequence_sha256"]
                    if h == row["human_reference_sha256"] or h in used:
                        continue
                    state = "drawn" if picked < target else "unused"
                    if state == "drawn":
                        used[h] = row["form_id"]
                        picked += 1
                    queue.append({**row, "shuffle_position": position, "draw_state": state})
                queues[key] = queue
    return queues


# --- main -----------------------------------------------------------------------------


def run(args) -> dict:
    cache = args.cache
    panel = read_tsv(PANEL)
    selected = json.loads(fetch("selected_cds", cache).read_text())
    human_cds = json.loads(fetch("human_cds", cache).read_text())
    gpi_genes = {
        r["gene_symbol"] for r in read_tsv(fetch("draw_gpi", cache)) if r["gpi_anchored"] == "1"
    }
    snapshot = pd.read_parquet(fetch("topology_review", cache))

    rows, fixed, extra = replay_draw(cache)
    queues = stratum_queues(rows, fixed, extra["policy"])
    proteins = extra["proteins"]

    def human_protein(gene: str, expected: str) -> str:
        seq = seqtools.translate(human_cds[f"lib:{gene}"]["cds"]).rstrip("*")
        if digest(seq) != expected:
            raise ValueError(f"Human library protein checksum mismatch: {gene}")
        return seq

    members = []
    for r in panel:
        enc = selected[r["form_id"]]
        seq = seqtools.translate(enc["cds"]).rstrip("*")
        if digest(seq) != r["sequence_sha256"] or digest(enc["cds"]) != r["selected_cds_sha256"]:
            raise ValueError(f"Selected CDS does not encode panel protein: {r['form_id']}")
        members.append(
            {
                **r,
                "role": "panel",
                "seq": seq,
                "human": human_protein(r["human_symbol"], r["human_reference_sha256"]),
            }
        )

    aligner = make_aligner()
    clans = load_clans(fetch("pfam_clans", cache))
    judgements = {r["form_id"]: r for r in read_tsv(JUDGEMENTS)} if JUDGEMENTS.exists() else {}
    for form_id, judged in judgements.items():
        if tuple(judged) != JUDGEMENT_COLUMNS or not all(judged.values()):
            raise ValueError(f"Judgement for {form_id} must fill {JUDGEMENT_COLUMNS}")
    unknown = set(judgements) - {m["form_id"] for m in members}
    if unknown:
        raise ValueError(f"Judgements for forms outside the panel: {sorted(unknown)}")

    def evaluate(items: list[dict]) -> list[dict]:
        human_seqs = {f"human|{m['human_symbol']}": m["human"] for m in items}
        out = []
        topo_h = topology_rows(human_seqs, "human", snapshot)
        topo_o = {}
        for species in ("mouse", "cyno"):
            sub = {m["form_id"]: m["seq"] for m in items if m["species"] == species}
            if sub:
                topo_o |= topology_rows(sub, species, snapshot)
        pfam = pfam_scan(
            {**human_seqs, **{m["form_id"]: m["seq"] for m in items}}, args.pfam_hmm, clans
        )
        for m in items:
            hk = f"human|{m['human_symbol']}"
            evidence = compare(
                m["human"],
                m["seq"],
                topo_h[hk],
                topo_o[m["form_id"]],
                pfam[hk],
                pfam[m["form_id"]],
                m["human_symbol"] in gpi_genes,
                m["anchor_call"],
                aligner,
            )
            if (m["topology_class"] == "GPI") != evidence["human_gpi"]:
                raise ValueError(f"Human GPI class disagrees with GPI set: {m['form_id']}")
            evidence["human_topology_version"] = topo_h[hk]["topology_version"]
            evidence["ortholog_topology_version"] = topo_o[m["form_id"]]["topology_version"]
            rule, why = recommend(evidence)
            judged = judgements.get(m["form_id"])
            if judged and judged["rule_recommendation"] != rule:
                raise ValueError(
                    f"Judgement for {m['form_id']} was written against a different rule result"
                )
            out.append(
                {
                    **m,
                    **evidence,
                    "rule_recommendation": rule,
                    "rule_basis": why,
                    "recommendation": judged["recommendation"] if judged else rule,
                    "recommendation_basis": judged["reason"] if judged else why,
                    "decided_by": judged["decided_by"] if judged else RULE_DECIDED_BY,
                    "implemented_by": (judged["implemented_by"] if judged else RULE_IMPLEMENTED_BY),
                    "pfam_hits": pfam,
                }
            )
        return out

    reviewed = evaluate(members)

    # Alternates: next unused candidates, in seeded shuffle order, of each stratum holding
    # a "replace" or "ask", reviewed in batches until enough of them keep the architecture.
    need = defaultdict(int)
    for m in reviewed:
        if m["recommendation"] in {"replace", "ask"}:
            need[(m["species"], m["band"], m["topology_class"])] += 1
    panel_hashes = {m["sequence_sha256"] for m in members}
    unused = {
        key: [
            q
            for q in queues[key]
            if q["draw_state"] == "unused" and q["sequence_sha256"] not in panel_hashes
        ]
        for key in need
    }
    alternates, offset = [], dict.fromkeys(need, 0)
    while True:
        batch = []
        for key, n in sorted(need.items()):
            keeps = sum(
                a["recommendation"] == "keep"
                for a in alternates
                if (a["species"], a["band"], a["topology_class"]) == key
            )
            if keeps >= n or offset[key] >= len(unused[key]):
                continue
            size = max(ALTERNATES_PER_STRATUM, 2 * n)
            for rank, q in enumerate(
                unused[key][offset[key] : offset[key] + size], offset[key] + 1
            ):
                seq = proteins[f"{q['species']}|{q['human_symbol']}|{q['picked_acc']}"]
                if digest(seq) != q["sequence_sha256"]:
                    raise ValueError("Alternate sequence checksum mismatch")
                batch.append(
                    {
                        **q,
                        "role": "alternate",
                        "queue_rank": rank,
                        "seq": seq,
                        "human": human_protein(q["human_symbol"], q["human_reference_sha256"]),
                    }
                )
            offset[key] += size
        if not batch:
            break
        alternates += evaluate(batch)

    # Assign replacements in queue order: first architecture-conserving ("keep") alternates.
    proposals = {}
    for key in sorted(need):
        flagged = sorted(
            (
                m
                for m in reviewed
                if m["recommendation"] in {"replace", "ask"}
                and (m["species"], m["band"], m["topology_class"]) == key
            ),
            key=lambda m: (m["recommendation"] != "replace", m["form_id"]),
        )
        pool = [
            a
            for a in alternates
            if (a["species"], a["band"], a["topology_class"]) == key
            and a["recommendation"] == "keep"
        ]
        pool.sort(key=lambda a: a["queue_rank"])
        for m, alt in zip(flagged, pool, strict=False):  # pool may be longer
            proposals[m["form_id"]] = alt["form_id"]
    for m in reviewed:
        m["proposed_replacement"] = proposals.get(m["form_id"], "")
        m["named_case"] = (m["species"], m["human_symbol"]) in NAMED_CASES
    if sum(m["named_case"] for m in reviewed) != len(NAMED_CASES):
        raise ValueError("A named review case is missing from the panel")

    columns = [
        "species",
        "band",
        "topology_class",
        "human_symbol",
        "ortholog_symbol",
        "form_id",
        "picked_acc",
        "picked_score",
        "sequence_sha256",
        "human_reference_sha256",
        "human_length",
        "ortholog_length",
        "length_ratio",
        "human_aligned_fraction",
        "aligned_identity",
        "human_unmatched_n",
        "human_unmatched_internal",
        "human_unmatched_c",
        "ortholog_only_n",
        "ortholog_only_internal",
        "ortholog_only_c",
        "human_unmatched_segments",
        "ortholog_only_segments",
        "human_topology_label",
        "ortholog_topology_label",
        "human_signal_peptide_length",
        "ortholog_signal_peptide_length",
        "human_tm_helices",
        "ortholog_tm_helices",
        "human_tm_helices_unmatched",
        "human_n_orientation",
        "ortholog_n_orientation",
        "human_gpi",
        "ortholog_anchor_call",
        "human_pfam",
        "ortholog_pfam",
        "pfam_domains_lost",
        "pfam_domains_lost_with_indel",
        "pfam_domain_lost_with_indel",
        "pfam_domains_matched_by_clan_or_noise_tier",
        "pfam_domains_partially_aligned",
        "pfam_domains_gained",
        "unmatched_segment_in_human_feature",
        *ARCHITECTURE_FLAGS,
        "architecture_flags",
        "architecture_conserved",
        "rule_recommendation",
        "rule_basis",
        "recommendation",
        "recommendation_basis",
        "decided_by",
        "implemented_by",
        "proposed_replacement",
        "named_case",
        "human_topology_version",
        "ortholog_topology_version",
    ]
    alt_columns = [
        "queue_rank",
        "shuffle_position",
        *[c for c in columns if c not in {"proposed_replacement", "named_case"}],
    ]
    reviewed.sort(key=lambda m: (m["species"], m["band"], m["topology_class"], m["human_symbol"]))
    alternates.sort(key=lambda m: (m["species"], m["band"], m["topology_class"], m["queue_rank"]))
    write(args.out_dir / "ortholog_correspondence_review.tsv", reviewed, columns)
    write(args.out_dir / "ortholog_correspondence_alternates.tsv", alternates, alt_columns)

    pfam_rows = []
    for item in [*reviewed, *alternates]:
        for k, label in (
            (f"human|{item['human_symbol']}", "human"),
            (item["form_id"], item["species"]),
        ):
            for d in item["pfam_hits"].get(k, []):
                pfam_rows.append({"key": k, "organism": label, "role": item["role"], **d})
    pfam_cols = [
        "key",
        "organism",
        "role",
        "pfam",
        "pfam_version",
        "name",
        "clan",
        "start",
        "end",
        "score",
        "tier",
    ]
    pfam_df = pd.DataFrame(pfam_rows, columns=pfam_cols).drop_duplicates(
        subset=["key", "pfam", "start", "end"]
    )
    pfam_df.sort_values(["key", "start", "pfam"]).to_csv(
        args.out_dir / "ortholog_correspondence_pfam_hits.tsv", sep="\t", index=False
    )

    counts = defaultdict(int)
    for m in reviewed:
        counts[m["recommendation"]] += 1
    summary = {
        "rule_version": RULE_VERSION,
        "source_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True
        ).stdout.strip(),
        "panel": str(PANEL.relative_to(HERE.parents[2])),
        "panel_sha256": sha256_file(PANEL),
        "inputs": {
            k: {"source": f"{b}://{key}", "sha256": h} for k, (b, key, h, _) in INPUTS.items()
        },
        "pfam_release": "Pfam 39.0 (EBI FTP releases/Pfam39.0); human domains at GA, ortholog counterparts same family or clan at >= NC; clan overlaps resolved as pfam_scan.pl",
        "topology_tool": topology.TOOL_V1,
        "alignment": "Bio.Align.PairwiseAligner global, BLOSUM62, open -11, extend -1, end gaps free",
        "members": len(reviewed),
        "by_species": {s: sum(m["species"] == s for m in reviewed) for s in ("mouse", "cyno")},
        "architecture_conserved": sum(m["architecture_conserved"] for m in reviewed),
        "flag_counts": {k: sum(bool(m[k]) for m in reviewed) for k in ARCHITECTURE_FLAGS},
        "unmatched_segment_in_human_feature": sum(
            m["unmatched_segment_in_human_feature"] for m in reviewed
        ),
        "recommendations": dict(sorted(counts.items())),
        "replacements_proposed": len(proposals),
        "replacements_without_conserving_alternate": sorted(
            m["form_id"]
            for m in reviewed
            if m["recommendation"] in {"replace", "ask"} and not m["proposed_replacement"]
        ),
        "alternates_reviewed": len(alternates),
        "draw_replay_reproduces_panel": True,
    }
    plot(reviewed, args.out_dir)
    (args.out_dir / "ortholog_correspondence_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def plot(reviewed: list[dict], out_dir: Path) -> None:
    """One dot per panel member: aligned fraction of the human protein vs length ratio."""
    import matplotlib.pyplot as plt

    from deliverome_analysis.plotting import get_color, save_figure, setup_plotting_style

    setup_plotting_style()
    colors = {
        "keep": get_color("neutral"),
        "ask": get_color("accent"),
        "replace": get_color("primary"),
    }
    markers = {"mouse": "o", "cyno": "s"}
    fig, ax = plt.subplots(figsize=(11, 8))
    for rec in ("keep", "ask", "replace"):
        for species in ("mouse", "cyno"):
            sub = [m for m in reviewed if m["recommendation"] == rec and m["species"] == species]
            if not sub:
                continue
            ax.scatter(
                [m["length_ratio"] for m in sub],
                [m["human_aligned_fraction"] for m in sub],
                c=colors[rec],
                marker=markers[species],
                s=70 if rec != "keep" else 40,
                alpha=0.9 if rec != "keep" else 0.5,
                edgecolors="white",
                linewidths=0.5,
                label=f"{rec}, {species} (n={len(sub)})",
            )
    # Label only the three cases named in the review brief; the other flagged forms sit in
    # the main cluster, which is the point: length and coverage do not find them.
    for m in reviewed:
        if m["named_case"]:
            ax.annotate(
                f"{m['species']} {m['human_symbol']}",
                (m["length_ratio"], m["human_aligned_fraction"]),
                xytext=(6, -4),
                textcoords="offset points",
                fontsize=14,
            )
    ax.axvline(1.0, color=get_color("line"), linewidth=1, zorder=0)
    ax.set_xlabel("Ortholog length / human library protein length")
    ax.set_ylabel("Fraction of human residues aligned")
    ax.legend(fontsize=12, loc="lower right", frameon=False)
    save_figure(fig, "ortholog_correspondence_review", output_dir=out_dir / "figures")
    plt.close(fig)


def write(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path.home() / ".cache/deliverome-analysis/s4e9-ortholog-correspondence",
    )
    parser.add_argument(
        "--pfam-hmm",
        type=Path,
        default=None,
        help="Local decompressed Pfam-A.hmm (checked against the pinned hash)",
    )
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args()
    if args.pfam_hmm is None:
        args.pfam_hmm = fetch("pfam", args.cache)
    elif sha256_file(args.pfam_hmm) != INPUTS["pfam"][3]:
        raise ValueError("Supplied Pfam-A.hmm does not match the pinned release checksum")
    print(json.dumps(run(args), indent=2))


if __name__ == "__main__":
    main()
