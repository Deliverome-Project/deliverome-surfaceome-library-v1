"""Assemble the queue of tagged L1 proteins to refold, tiered by decision value.

Becca, 2026-11-11: "for proteins that have a clear position at the terminus refold those to
double check and for n term proteins run signal reprediction and disorder etc as we had
planned to select the site then also refold those" ($230 of Modal for this folder).
Becca, 2026-11-12, delegating: "use your best judgment to complete a new and improved
evidence based l1 completely by the morning" -- so the refold verdicts are decisions, and
when the night cannot fold everything, this folder chooses what to fold by the value of
the decision it informs. Implemented by Claude Code session l1-tag-refold.

**Inputs are each rung's own final export**, after PR #283 renumbered every member, keyed
by the tagged protein's SHA-256 (or design_id), and snapshotted with their SHA-256 into
``inputs/`` beside the published models:

====================  =========================================================  ======
rung                  input (owner session, PR)                                  tier
====================  =========================================================  ======
nterm                 ``nterm_arms_decided_2026-11-12`` (nterm-reviewer-package,  1 / 2 / 3
                      #254): order 1-591 = every Igk arm + every native-site arm
                      inside the predicted fold (a refold can overturn the call)
cterm                 ``l1_cterm_tagged_proteins`` (l1-cterm-build, #284)         1 / 2 / 3
type3                 ``constructs_l1_type3_proteins`` (l1-nterm-type3, #282)     1
multipass             ``res1-fold-inputs`` (l1-nterm-res1-multipass, #280),       2 / 3
                      mature chains (Igk leader removed)
loop                  ``l1_internal_refold_queue`` (ranker-pu-segment, #245)      1 / 2 / 3
native                members' native proteins with no model of exactly that      tier of
                      sequence (``tag_refold_native_models.tsv``)                  its construct
control, sp_control   untagged calibration folds (benchmark)                      0
====================  =========================================================  ======

Tier 1: the refold can change a decision -- N-terminal Igk / inside-fold arms; C-terminal
arms whose terminus is structured (native terminal pLDDT >= 85), the M13 family, KLRF1, the
TNF ligands, ITM2B and the three new priorities l1-cterm-build named; every type III design;
both loop sites of the 100 members whose rank-2 loop is close to rank 1.
Tier 2: a stratified sample of each rung's routine builds (100 each, every k-th by SHA-256;
for multipass, all 24 leader-free arms plus 76 Igk arms) to measure its failure rate. Tier 3: the rest.

One fold per tagged SHA-256; ``tag_refold_sources.tsv`` maps every rung's id onto it. A fold
already made under an earlier queue keeps its fold_id; new ones are ``t.<sha16>`` (tagged)
or ``n.<sha16>`` (native).

    uv run --locked python library/l1/orf-design/tag_refold/build_refold_inputs.py
"""

from __future__ import annotations

import csv
import difflib
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EXP))

import build_gpcr_rescue_arms as gpcr_arms  # noqa: E402

from deliverome_analysis import library_design as ld  # noqa: E402
from deliverome_analysis import plasmid_qc as qc  # noqa: E402
from deliverome_analysis import seqtools  # noqa: E402
from deliverome_analysis.provenance import read_pinned_text  # noqa: E402

OUT = EXP / "tag_refold_queue.tsv"  # tables live at the experiment root (.gitignore)
SOURCES_OUT = EXP / "tag_refold_sources.tsv"
NATIVE_INDEX = EXP / "tag_refold_native_models.tsv"
SNAPSHOT = Path.home() / ".cache" / "deliverome-analysis" / "tag-refold-inputs"

#: The p003 51-aa tandem, full length everywhere (Becca, 2026-11-10 and 2026-11-11 "keep full").
CASSETTE = qc.P003_CASSETTE_PEPTIDE
SPYTAG003 = qc.TAG_PEPTIDES["SpyTag003"]
ALFA = qc.TAG_PEPTIDES["P-ALFA-P"].strip("P")
assert len(CASSETTE) == 51 and SPYTAG003 in CASSETTE and ALFA in CASSETTE

#: Rung inputs, by file name. Each rung session wrote its export to the shared notes folder
#: or to S3; ``--inputs-dir`` ingests a folder holding them, uploads each file once to a
#: write-once key named by its SHA-256 and records it in ``tag_refold_inputs.tsv``. A
#: rebuild without ``--inputs-dir`` reads exactly those snapshots back.
INPUT_FILES = {
    "nterm_tsv": "nterm-refold-input-2026-11-12/nterm_arms_decided_2026-11-12.tsv",
    "nterm_fasta": "nterm-refold-input-2026-11-12/nterm_arms_decided_2026-11-12.protein.fasta",
    "cterm_tsv": "l1-cterm-refold-input-2026-11-12c/l1_cterm_tagged_proteins.tsv",
    "cterm_fasta": "l1-cterm-refold-input-2026-11-12c/l1_cterm_tagged_proteins.fasta",
    "nterm_req_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-12.tsv",
    "nterm_req_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-12.fasta",
    "nterm_req_b_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-12b.tsv",
    "nterm_req_b_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-12b.fasta",
    "spacer_tsv": "res1-spacer-comparison-2026-11-12.tsv",
    "nterm_req_c_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-13a.tsv",
    "nterm_req_c_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-13a.fasta",
    "nterm_req_d_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-13b.tsv",
    "nterm_req_d_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-13b.fasta",
    "snorkel_new_tsv": "l1-snorkel-refold-input-2026-11-13/new_snorkels_to_refold.tsv",
    "nterm_req_e_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-13c.tsv",
    "nterm_req_e_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-13c.fasta",
    "snorkel_new_fasta": "l1-snorkel-refold-input-2026-11-13/l1_cterm_tagged_proteins.fasta",
    "multipass_tsv": "res1-fold-inputs-2026-11-12.tsv",
    "nterm_req_f_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-14.tsv",
    "nterm_req_f_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-14.fasta",
    "nterm_req_g_tsv": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-15_weak_sp.tsv",
    "nterm_req_g_fasta": "nterm-refold-input-2026-11-12/nterm_refold_requests_2026-11-15_weak_sp.fasta",
}
#: Rung inputs that live on S3 already (snapshotted the same way).
INPUT_KEYS = {
    "type3_tsv": "sprint4/s4e9-orf-design-fg-panel/constructs_l1_type3_proteins.tsv",
    "loop_tsv": "sprint5/s5e13-tag-site-benchmarking/l1_internal_refold_queue.tsv",
    "loop_add_tsv": "sprint5/s5e13-tag-site-benchmarking/l1_internal_refold_queue_addendum.tsv",
    "loop_alt_tsv": "sprint5/s5e13-tag-site-benchmarking/l1_refold_alternatives_requested.tsv",
    "order_fasta": "sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-13-one-arm-renamed/l1_order_inserts_recoded.fasta",
    "multipass_mature_fasta": "sprint4/s4e9-orf-design-fg-panel/multipass_arms_mature.fasta",
    "rescue_mature_fasta": "sprint4/s4e9-orf-design-fg-panel/rescue_arms_mature.fasta",
    "rescue_round_fasta": "sprint5/s5e13-tag-site-benchmarking/rescue_round/2026-11-14/rescue_round_tagged.fasta",
    "rescue_round_tsv": "sprint5/s5e13-tag-site-benchmarking/rescue_round/2026-11-14/rescue_round_sites.tsv",
    "rescue_round_b_fasta": "sprint5/s5e13-tag-site-benchmarking/rescue_round/2026-11-14b/rescue_round_tagged.fasta",
    "rescue_round_b_tsv": "sprint5/s5e13-tag-site-benchmarking/rescue_round/2026-11-14b/rescue_round_sites.tsv",
    "order_rule_fasta": "sprint5/s5e13-tag-site-benchmarking/order_rule/2026-11-15/order_rule_fold_requests.fasta",
    "order_rule_tsv": "sprint5/s5e13-tag-site-benchmarking/order_rule/2026-11-15/order_rule_fold_requests.tsv",
    "arm_rev_fasta": "sprint5/s5e13-tag-site-benchmarking/arm_revisions/2026-11-15/arm_revisions_new_constructs.fasta",
    "arm_rev_tsv": "sprint5/s5e13-tag-site-benchmarking/arm_revisions/2026-11-15/arm_revisions_new_constructs.tsv",
    "agent_screen_fasta": "sprint5/s5e13-tag-site-benchmarking/agent_site_screen_to_fold.fasta",
    "agent_screen_tsv": "sprint5/s5e13-tag-site-benchmarking/agent_site_screen_to_fold.tsv",
    "annotation_rule_fasta": "sprint5/s5e13-tag-site-benchmarking/annotation_rule_to_fold.fasta",
    "annotation_rule_tsv": "sprint5/s5e13-tag-site-benchmarking/annotation_rule_to_fold.tsv",
    "annotation_rule_b_fasta": "sprint5/s5e13-tag-site-benchmarking/annotation_rule_to_fold_2026-11-15b.fasta",
    "annotation_rule_b_tsv": "sprint5/s5e13-tag-site-benchmarking/annotation_rule_to_fold_2026-11-15b.tsv",
}
INPUTS_TABLE = EXP / "tag_refold_inputs.tsv"
NATIVES_OUT = EXP / "tag_refold_natives.tsv"
TAGGED_OUT = EXP / "tag_refold_tagged.tsv"
INPUT_PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/tagged/2026-11-12/inputs"

#: l1-reviewer-package's PR #283 head when the queue was built: its manifest names the
#: members, and the CDS records below give their proteins (loops reference the
#: ortholog-policy record; the rest the current one).
RP_COMMIT = "5b97ddfc557de0a1c78737b100daf16a1f734202"
#: main when the queue was first built (before #283 merged on 2026-11-12): its manifest is
#: the other source of native proteins, read from git so a later merge cannot change it.
MAIN_COMMIT = "5b133d0b"
CDS_KEYS = (
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-ortholog-start-trims/candidate_cds.json",
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-12-reference-fixes/candidate_cds.json",
    "sprint4/s4e9-orf-design-fg-panel/library-candidates/2026-11-10-ortholog-cds-policy/candidate_cds.json",
    ld.BASE_CDS_KEY,
)

#: Tier-1 C-terminal genes and designs (l1-cterm-build's priority list, 2026-11-11/12).
CTERM_T1_GENES = {
    "MME",
    "MMEL1",
    "ECE1",
    "PHEX",
    "KEL",
    "KLRF1",
    "ITM2B",
    "FASLG",
    "CD40LG",
    "CD70",
    "LTA",
    "LTB",
    "TNF",
}
CTERM_T1_DESIGNS = {"l1-p3137", "l1-p3003", "l1-p3112"}
NTERM_T1_MAX_ORDER = 591
SAMPLE = 100

DELEGATED = "l1-tag-refold (delegated by Becca 2026-11-12)"

#: l1-cterm-build's snorkel arms for the 26 members whose designed internal sites all
#: failed refolding (ranker abfc7f70): native + 55-aa PDGFRB JM/TM module + cassette at the
#: C-terminus. Read from that session's commit, local when queued (2026-11-12), via git.
SNORKEL_COMMIT = "3ef7f910aff5510c452632b378a2985992436cb1"
CODON = seqtools.GENETIC_CODE
translate_insert = qc.translate_insert


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fasta(path: Path) -> list[tuple[str, str]]:
    return [(header, seq.upper()) for header, seq in seqtools.read_fasta(path)]


def blocks(native: str, tagged: str) -> list[tuple[int, int, int]]:
    """(tagged_start, native_start, length) runs of identical residues, 0-based.

    Only runs of at least 6 residues count, so a cassette residue that happens to equal a
    native one (a G, an S) is never mapped onto the native protein.
    """
    sm = difflib.SequenceMatcher(None, native, tagged, autojunk=False)
    return [(b.b, b.a, b.size) for b in sm.get_matching_blocks() if b.size >= 6]


def cassette_span(tagged: str) -> tuple[int, int]:
    i = tagged.find(CASSETTE)
    if i < 0 or tagged.find(CASSETTE, i + 1) >= 0:
        raise ValueError("cassette absent or present twice")
    return i, i + len(CASSETTE)


def leader_end(b: list[tuple[int, int, int]], c0: int) -> int | str:
    """End (tagged coordinates) of a non-native N-terminal leader of 8+ residues, if any."""
    if c0 < 8 or any(t0 < c0 for t0, _, _ in b):
        return ""
    return c0


def stratified(items: list, n: int, key) -> tuple[list, list]:
    """Every k-th item by ``key`` (a SHA-256): a deterministic, unbiased sample of n."""
    items = sorted(items, key=key)
    if len(items) <= n:
        return items, []
    step = len(items) / n
    idx = {int(i * step) for i in range(n)}
    return [x for i, x in enumerate(items) if i in idx], [
        x for i, x in enumerate(items) if i not in idx
    ]


class Construct(dict):
    pass


def construct(rung: str, source_id: str, design_id: str, gene: str, native: str, tagged: str,
              *, source: str, tier: int, order: float, reason: str, extra: dict | None = None,
              published_sha: str = "") -> Construct:  # fmt: skip
    if published_sha and sha(tagged) != published_sha:
        raise SystemExit(f"{rung} {source_id}: tagged protein does not match its published SHA-256")
    b = blocks(native, tagged)
    c0, c1 = cassette_span(tagged)
    before = [n0 + n for t0, n0, n in b if t0 + n <= c0]
    return Construct(
        rung=rung, source_id=source_id, design_id=design_id, gene_symbol=gene, tier=tier,
        order=order, native_aa=len(native), tagged_aa=len(tagged),
        native_mapped_aa=sum(n for _, _, n in b), cassette_start=c0, cassette_end=c1,
        insert_after_native_residue=max(before) if before else 0,
        leader_end=leader_end(b, c0), native_sha256=sha(native), tagged_sha256=sha(tagged),
        blocks=json.dumps(b, separators=(",", ":")), native_seq=native, tagged_seq=tagged,
        source=source, reason=reason, **(extra or {}),
    )  # fmt: skip


#: Request files: a rung asks for folds by writing one, and may later overwrite or empty it
#: once served (ranker-pu-segment emptied l1_refold_alternatives_requested.tsv on
#: 2026-11-12, and a re-ingest silently dropped its 52 sites from the verdicts). For these
#: every snapshot ever ingested is kept and their rows are unioned; an empty file is never
#: recorded. Base exports (the rung tables) are replaced by the latest snapshot.
REQUEST_ROLES = {
    "loop_add_tsv", "loop_alt_tsv", "nterm_req_tsv", "nterm_req_fasta",
    "nterm_req_b_tsv", "nterm_req_b_fasta", "spacer_tsv",
    "nterm_req_c_tsv", "nterm_req_c_fasta", "nterm_req_d_tsv", "nterm_req_d_fasta",
    "snorkel_new_tsv", "snorkel_new_fasta", "multipass_mature_fasta", "rescue_mature_fasta",
    "nterm_req_e_tsv", "nterm_req_e_fasta", "nterm_req_f_tsv", "nterm_req_f_fasta",
    "nterm_req_g_tsv", "nterm_req_g_fasta",
    "rescue_round_fasta", "rescue_round_tsv", "rescue_round_b_fasta", "rescue_round_b_tsv",
    "order_rule_fasta", "order_rule_tsv", "arm_rev_fasta", "arm_rev_tsv",
    "agent_screen_fasta", "agent_screen_tsv", "annotation_rule_fasta", "annotation_rule_tsv",
    "annotation_rule_b_fasta", "annotation_rule_b_tsv",
}  # fmt: skip


def resolve_inputs(inputs_dir: Path | None) -> dict[str, list[Path]]:
    """Local paths of every rung input (a list per role), from a folder or from S3.

    With ``inputs_dir`` each file is snapshotted to a write-once S3 key and recorded in
    ``tag_refold_inputs.tsv``; without it, the recorded snapshots are read back.
    """
    from deliverome_analysis.s3 import processed

    SNAPSHOT.mkdir(parents=True, exist_ok=True)
    rows = tsv(INPUTS_TABLE) if INPUTS_TABLE.exists() else []
    if inputs_dir is not None:
        found: dict[str, Path] = {}
        for role, name in INPUT_FILES.items():
            if not (inputs_dir / name).is_file():
                raise SystemExit(f"{name}: not found under {inputs_dir}")
            found[role] = inputs_dir / name
        for role, key in INPUT_KEYS.items():
            found[role] = Path(processed.download(key))
        for role, path in sorted(found.items()):
            data = path.read_bytes()
            digest = hashlib.sha256(data).hexdigest()
            name = Path(INPUT_FILES.get(role) or INPUT_KEYS[role]).name
            key = f"{INPUT_PREFIX}/{digest[:12]}-{name}"
            if role in REQUEST_ROLES:
                if len(data.strip().splitlines()) < 2 or any(
                    r["role"] == role and r["sha256"] == digest for r in rows
                ):
                    continue  # empty (served) or already recorded: keep what we have
            else:
                rows = [r for r in rows if r["role"] != role]
            if not processed.exists(key):
                processed.upload(path, key)
            rows.append({"role": role, "file": name, "sha256": digest, "key": key})
        rows.sort(key=lambda r: (r["role"], r["key"]))
        with INPUTS_TABLE.open("w", newline="") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["role", "file", "sha256", "key"],
                delimiter="\t",
                lineterminator="\n",
            )
            w.writeheader()
            w.writerows(rows)
    out: dict[str, list[Path]] = {}
    for r in rows:
        out.setdefault(r["role"], []).append(Path(processed.download(r["key"])))
    return out


def rows_of(paths: list[Path]) -> list[dict[str, str]]:
    """Rows of every snapshot of a request file, de-duplicated (later snapshots win)."""
    seen: dict[str, dict[str, str]] = {}
    for p in paths:
        for r in tsv(p):
            seen[json.dumps(r, sort_keys=True)] = r
    return list(seen.values())


def fasta_of(paths: list[Path]) -> dict[str, str]:
    out: dict[str, str] = {}
    for p in paths:
        out.update({h.split()[0]: sq for h, sq in fasta(p)})
    return out


def native_lookup() -> tuple[dict[str, str], dict[str, str]]:
    """(design_id -> protein under the current CDS record, sha256 -> protein, any record)."""

    from deliverome_analysis import library_design as ld

    manifest = SNAPSHOT / f"candidate_manifest_tranches.{RP_COMMIT[:8]}.tsv"
    SNAPSHOT.mkdir(parents=True, exist_ok=True)
    manifest.write_text(
        read_pinned_text(
            RP_COMMIT,
            "library/l1/orf-design/candidate_manifest_tranches.tsv",
            repo=REPO,
        )
    )
    main_manifest = SNAPSHOT / f"candidate_manifest_tranches.{MAIN_COMMIT}.tsv"
    main_manifest.write_text(
        read_pinned_text(
            MAIN_COMMIT,
            "library/l1/orf-design/candidate_manifest_tranches.tsv",
            repo=REPO,
        )
    )
    did = {r["candidate_id"]: r["design_id"] for r in tsv(manifest)}
    by_design: dict[str, str] = {}
    by_sha: dict[str, str] = {}
    saved = ld.CDS_KEY
    try:
        for i, key in enumerate(CDS_KEYS):
            ld.CDS_KEY = key
            for man in (manifest, main_manifest):
                for cid, v in ld.designed_proteins(man, library=ld.LIBRARY_L1).items():
                    by_sha.setdefault(sha(v["sequence"]), v["sequence"])
                    if i == 0 and man == manifest:
                        by_design[did[cid]] = v["sequence"]
    finally:
        ld.CDS_KEY = saved
    return by_design, by_sha


def build(inp_all: dict[str, list[Path]]) -> list[Construct]:
    inp = {role: paths[-1] for role, paths in inp_all.items()}
    by_design, by_sha = native_lookup()
    # design_id for rungs that key by member id (type III, loops): the #283 design whose
    # designed protein is this native, by sequence (member ids were renumbered twice).
    design_of: dict[str, str] = {}
    for d, sq in by_design.items():
        h = sha(sq)
        design_of[h] = f"{design_of[h]};{d}" if h in design_of else d
    out: list[Construct] = []

    # --- N after SP: 1,307 decided arms ------------------------------------------------
    nt_tsv, nt_fa = inp["nterm_tsv"], inp["nterm_fasta"]
    prot = {h.split()[0]: s for h, s in fasta(nt_fa)}
    nterm = []
    for r in tsv(nt_tsv):
        tagged = prot[r["protein_sha256"]]
        native = by_design[r["design_id"]]
        k = int(r["cassette_after_residue"])
        assert tagged[k + 51 :] == native[int(r["native_after_residue"]) :], r["design_id"]
        order = int(r["order"])
        nterm.append(construct(
            "nterm", r["design_id"], r["design_id"], r["gene_symbol"], native, tagged,
            source=f"{Path(INPUT_FILES["nterm_tsv"]).name} (nterm-reviewer-package, PR #254), arm {r['arm']}",
            tier=1 if order <= NTERM_T1_MAX_ORDER else 3, order=order,
            reason="Igk or inside-fold arm: the refold can overturn the call" if order <= NTERM_T1_MAX_ORDER else "before-fold arm",
            extra={"arm": r["arm"], "inside_predicted_fold": r["inside_predicted_fold"]},
            published_sha=r["protein_sha256"],
        ))  # fmt: skip
    routine = [c for c in nterm if c["tier"] == 3]
    for c in stratified(routine, SAMPLE, key=lambda c: c["tagged_sha256"])[0]:
        c["tier"], c["reason"] = 2, "before-fold arm, stratified sample of the routine builds"
    out += nterm

    # --- C-terminal: 214 -----------------------------------------------------------------
    ct_tsv, ct_fa = inp["cterm_tsv"], inp["cterm_fasta"]
    prot = {h.split("|")[0]: s for h, s in fasta(ct_fa)}
    cterm = []
    for r in tsv(ct_tsv):
        tagged = prot[r["design_id"]]
        a, b = int(r["cassette_start_aa"]) - 1, int(r["cassette_end_aa"])
        snorkel = r["placement"].startswith("snorkel")
        # A snorkel arm is native + a 55-aa PDGFRB JM/TM module + the cassette, at the C-terminus.
        native = tagged[: int(r["native_aa"])] if snorkel else tagged[:a] + tagged[b:]
        stalk = "stalk" in r["placement"]  # a failed C-term arm moved to its stalk site: decides it
        t1 = stalk or (not snorkel and ("structured_terminus" in r["review_flags"] or r["gene_symbol"] in CTERM_T1_GENES
              or r["gene_symbol"].startswith("TNFSF") or r["design_id"] in CTERM_T1_DESIGNS))  # fmt: skip
        cterm.append(construct(
            "cterm", r["design_id"], r["design_id"], r["gene_symbol"], native, tagged,
            source=f"{Path(INPUT_FILES["cterm_tsv"]).name} (l1-cterm-build, PR #284)",
            tier=1 if t1 else (2 if snorkel else 3),
            order=0 if r["design_id"] in CTERM_T1_DESIGNS else 1,
            reason=("structured terminus or named priority" if t1 else
                    "snorkel arm (TM module + cassette): folded as budget allows" if snorkel else "routine C-terminal arm"),
            extra={"review_flags": r["review_flags"], "terminal_plddt_native": r["terminal_plddt_native"]},
            published_sha=r["protein_sha256"],
        ))  # fmt: skip
    for c in stratified(
        [c for c in cterm if c["tier"] == 3],
        SAMPLE - sum(c["tier"] == 2 for c in cterm),
        key=lambda c: c["tagged_sha256"],
    )[0]:
        c["tier"], c["reason"] = 2, "routine C-terminal arm, stratified sample"
    out += cterm

    # --- type III: 44, all tier 1 ----------------------------------------------------------
    t3_tsv = inp["type3_tsv"]
    for r in tsv(Path(t3_tsv)):
        igk = r["design"] == "igk_tag"
        # Igk designs: fold the mature chain (the leader is cleaved); tag_alone: the precursor.
        tagged = r["mature"] if igk else r["precursor"]
        shift = len(r["precursor"]) - len(tagged)
        tm = r["native_tm1_span"]
        tm_t = ""
        if tm:
            s, e = (int(x) for x in tm.split("-"))
            tm_t = f"{s - shift}-{e - shift}"
        did = r.get("design_id") or design_of.get(r["native_sha256"], "")
        out.append(construct(
            "type3", did or r["precursor_sha256"], did, r["gene_symbol"], r["native"], tagged,
            source="constructs_l1_type3_proteins.tsv (l1-nterm-type3, PR #282), "
                   + ("mature chain" if igk else "precursor"),
            tier=1, order=0, reason="type III design: N-terminal cassette next to a single TM",
            extra={"tm1_span_folded": tm_t, "design": r["design"]},
        ))  # fmt: skip

    # --- residue-1 multipass: 558 mature chains, tier 2 sample / tier 3 -------------------
    mp = []
    for r in tsv(inp["multipass_tsv"]):
        mature = r["mature_sequence"]
        c = int(r["cassette_after_mature_residue"])
        native = "M" + mature[c + 51 :]
        if sha(native) != r["native_sha256"]:
            native = by_sha.get(r["native_sha256"], "")
        if not native:
            print(f"  multipass {r['construct_id']}: native not found; not queued", file=sys.stderr)
            continue
        mp.append(construct(
            "multipass", r["construct_id"], r["design_id"], r["gene_symbol"], native, mature,
            source=f"{INPUT_FILES['multipass_tsv']} (l1-nterm-res1-multipass, PR #280), mature chain",
            tier=3, order=1, reason="routine residue-1 multipass arm", extra={"design": r["design"]},
            published_sha=r["mature_sha256"],
        ))  # fmt: skip
    # Stratified by design as well as hash (l1-nterm-res1-multipass's suggestion, 2026-11-12):
    # every leader-free arm (24, too few for a hash sample to compare) plus a hash sample of
    # the Igk arms, 100 in all.
    alone = [c for c in mp if c["design"] == "tag_alone"]
    igk = [c for c in mp if c["design"] != "tag_alone"]
    for c in alone:
        c["tier"], c["reason"] = 2, "residue-1 multipass, every leader-free arm"
    for c in stratified(igk, SAMPLE - len(alone), key=lambda c: c["tagged_sha256"])[0]:
        c["tier"], c["reason"] = 2, "residue-1 multipass, stratified sample of the Igk arms"
    out += mp

    # --- internal loops: 659 -----------------------------------------------------------------
    loops = tsv(inp["loop_tsv"])
    two = {r["candidate_id"] for r in loops if r["rank"] == "2"}
    lp = []
    for r in loops:
        native = by_sha.get(r["sequence_sha256"])
        if not native:
            print(
                f"  loop {r['candidate_id']} rank {r['rank']}: native not found; not queued",
                file=sys.stderr,
            )
            continue
        k = int(r["insert_after_residue"])
        tagged = native[:k] + CASSETTE + native[k:]
        t1 = r["candidate_id"] in two
        lp.append(construct(
            "loop", f"{r['candidate_id']}:rank{r['rank']}", design_of.get(r["sequence_sha256"], ""),
            r["gene"], native, tagged,
            source="l1_internal_refold_queue.tsv (ranker-pu-segment, PR #245)",
            tier=1 if t1 else 3, order=int(r["rank"]),
            reason="rank-1 vs close rank-2: the refold can pick the site" if t1 else "rank-1 loop site",
            extra={"rank": r["rank"], "snorkel_candidate": r["snorkel_candidate"], "loop_member": r["candidate_id"]},
            published_sha=r["tagged_protein_sha256"],
        ))  # fmt: skip
    for c in stratified(
        [c for c in lp if c["tier"] == 3], SAMPLE, key=lambda c: c["tagged_sha256"]
    )[0]:
        c["tier"], c["reason"] = 2, "rank-1 loop site, stratified sample"
    out += lp

    # --- loop addendum: 19 decided sites missing from the queue, each its member's only site
    for r in rows_of(inp_all["loop_add_tsv"]):
        native = by_sha.get(r["sequence_sha256"])
        if not native:
            print(
                f"  loop addendum {r['design_id']}: native not found; not queued", file=sys.stderr
            )
            continue
        k = int(r["insert_after_residue"])
        out.append(construct(
            "loop", f"{r['member_id']}:rank{r['rank']}", r["design_id"], r["gene"], native,
            native[:k] + CASSETTE + native[k:],
            source="l1_internal_refold_queue_addendum.tsv (ranker-pu-segment, PR #245)",
            tier=1, order=0, reason="decided loop site, the member's only one",
            extra={"rank": r["rank"], "loop_member": r["member_id"]},
            published_sha=r["tagged_protein_sha256"],
        ))  # fmt: skip

    # --- loop alternatives for picks that failed refold (ranker-pu-segment, 2026-11-12):
    # each is now the only thing deciding its member, so tier 1, ahead of everything else.
    for r in rows_of(inp_all["loop_alt_tsv"]):
        native = by_sha.get(r["sequence_sha256"])
        if not native:
            print(
                f"  loop alternative {r['design_id']}: native not found; not queued",
                file=sys.stderr,
            )
            continue
        k = int(r["insert_after_residue"])
        out.append(construct(
            "loop", f"{r['design_id']}:rank{r['rank']}:after{k}", r["design_id"], r["gene"], native,
            native[:k] + CASSETTE + native[k:],
            source="l1_refold_alternatives_requested.tsv (ranker-pu-segment, PR #245)",
            tier=1, order=-1, reason=f"alternative loop site: {r['why']}",
            extra={"rank": r["rank"]}, published_sha=r["tagged_protein_sha256"],
        ))  # fmt: skip

    # --- snorkel arms (l1-cterm-build), tier 1: the design that now decides these members
    def show(rel: str) -> str:
        return read_pinned_text(SNORKEL_COMMIT, f"library/l1/orf-design/{rel}", repo=REPO)

    fa: dict[str, str] = {}
    name = None
    for line in show("constructs_l1_cterm.fasta").splitlines():
        if line.startswith(">"):
            name = line[1:].split("|")[0].split()[0]
            fa[name] = ""
        elif name:
            fa[name] += line.strip().upper()
    for r in csv.DictReader(show("constructs_l1_cterm.tsv").splitlines(), delimiter="\t"):
        if "snorkel" not in r["placement"]:
            continue
        tagged = translate_insert(
            fa[r["candidate_id"]] if r["candidate_id"] in fa else fa[r["design_id"]]
        )
        native = tagged[: int(r["native_aa"])]
        out.append(construct(
            "cterm", r["design_id"], r["design_id"], r["gene_symbol"], native, tagged,
            source=f"constructs_l1_cterm.tsv @ {SNORKEL_COMMIT[:8]} (l1-cterm-build, PR #284), snorkel",
            tier=1, order=0, reason="snorkel arm: the members' internal sites all failed refolding",
            extra={"review_flags": "snorkel"}, published_sha=r["protein_sha256"],
        ))  # fmt: skip

    # --- built constructs on the draft-L1 integration branch whose protein was never folded
    # (verify_ordered_constructs.py): queued at tier 1 so every built member has a verdict.
    unf = EXP / "tag_refold_ordered_unfolded.tsv"
    for r in tsv(unf) if unf.exists() else []:
        native = r["native_seq"] or by_sha.get(r["native_sha256"], "")
        if not native:
            print(f"  ordered {r['construct']}: native not found; not queued", file=sys.stderr)
            continue
        rung = {"internal": "loop"}.get(r["rung"], r["rung"])
        out.append(construct(
            rung, r["construct"], r["design_id"], r["gene_symbol"], native, r["tagged_seq"],
            source="constructs_l1_* @ #299 5cfd91c0 (built, never folded)",
            tier=1, order=0, reason="built construct with no fold yet (ordered-construct check)",
            extra={"review_flags": "ordered_check"},
        ))  # fmt: skip

    # --- multipass spacer comparison (l1-nterm-res1-multipass, 2026-11-12): the arm with
    # GGGGSGGGGS between the cassette and native residue 2, for the 12 members whose
    # native-region lDDT dropped. A separate rung so it never overwrites the arm's verdict.
    for r in rows_of(inp_all.get("spacer_tsv", [])):
        mature = r["mature_sequence"]
        assert sha(mature) == r["mature_sha256"], r["construct_id"]
        native = by_sha.get(r["native_sha256"]) or "M" + mature[2 + 51 + 10 :]
        assert sha(native) == r["native_sha256"], r["construct_id"]
        out.append(construct(
            "multipass_spacer", r["construct_id"], r["design_id"], r["gene_symbol"], native, mature,
            source="res1-spacer-comparison-2026-11-12.tsv (l1-nterm-res1-multipass, PR #280)",
            tier=1, order=0, reason=f"spacer comparison against {r['compare_to_construct_id']}",
            extra={"design": r["design"]},
        ))  # fmt: skip

    # --- mature-chain files (multipass and expression-rescue arms, PR #280 / #306): header
    # "<construct_id> <member> <gene> native_sha256=<hash>". Constructs the multipass input
    # does not already carry (the SLC61A1 pair) and every rescue arm are tier 1.
    known_mp = {r["construct_id"] for r in tsv(inp["multipass_tsv"])}
    for role, rung in (("multipass_mature_fasta", "multipass"), ("rescue_mature_fasta", "rescue")):
        for path in inp_all.get(role, []):
            for head, mature in fasta(path):
                cid = head.split()[0]
                if rung == "multipass" and cid in known_mp:
                    continue
                # header hash field: native_sha256= (first files) or reference_sha256= (renamed
                # 2026-11-13); either may be absent
                h = next(
                    (
                        t.split("=", 1)[1]
                        for t in head.split()
                        if t.split("=")[0].endswith("sha256")
                    ),
                    "",
                )
                c1 = cassette_span(mature)[1]
                # The header hash is checked first. The rescue file's hashes match no protein
                # (2026-11-13), so otherwise a native is accepted only if it is a known designed
                # protein: Met + the arm's tail, or (signal-peptide swap) the single designed
                # protein ending in the arm's tail. The mismatch is recorded on the row.
                tail = mature[c1:]
                note = ""
                native = by_sha.get(h, "")
                if not native:
                    known = set(by_sha.values())
                    cands = [sq for sq in known if sq.endswith(tail)]
                    native = (
                        ("M" + tail)
                        if ("M" + tail) in known
                        else (cands[0] if len(cands) == 1 else "")
                    )
                    note = f"; header native_sha256 {h[:12]} matches no protein, native = designed protein {sha(native)[:12]}"
                if not native:
                    print(f"  {cid}: native not resolvable; not queued", file=sys.stderr)
                    continue
                out.append(construct(
                    rung, cid, cid.split(":")[0], head.split()[2], native, mature,
                    source=f"{path.name.split('-', 1)[-1]} (l1-nterm-res1-multipass), mature chain",
                    tier=1, order=0,
                    reason=("new arm" if rung == "multipass" else "expression-rescue arm") + note,
                    extra={"design": cid.split(":")[-1]},
                ))  # fmt: skip

    # --- new snorkel arms (l1-cterm-build, ranker #245 d04b1cb1)
    if "snorkel_new_tsv" in inp_all:
        sfa = {h.split("|")[0].split()[0]: sq for h, sq in fasta(inp_all["snorkel_new_fasta"][-1])}
        for r in rows_of(inp_all["snorkel_new_tsv"]):
            tagged = sfa.get(r["design_id"]) or sfa[r["candidate_id"]]
            out.append(construct(
                "cterm", r["design_id"], r["design_id"], r["gene_symbol"], tagged[: int(r["native_aa"])], tagged,
                source="new_snorkels_to_refold.tsv (l1-cterm-build)", tier=1, order=0,
                reason="new snorkel arm", extra={"review_flags": "snorkel"},
            ))  # fmt: skip

    # --- held-member rescue round (ranker-pu-segment, 2026-11-14): internal ectodomain sites,
    # at most two per held member; header "<design>:after<N> ... tagged_sha256= native_sha256=".
    # Judged with the internal-loop bars, so they sit on the loop rung.
    # The revised list (2026-11-14b: 74 constructs, VTCN1 added) is read beside the first; folds
    # of the 46 it dropped are kept, and the ranker joins only on the revised file's SHA-256s.
    for path in inp_all.get("rescue_round_fasta", []) + inp_all.get("rescue_round_b_fasta", []):
        for head, tagged in fasta(path):
            f = dict(t.split("=", 1) for t in head.split()[1:] if "=" in t)
            sid = head.split()[0]
            k = int(sid.rsplit(":after", 1)[1])
            assert tagged[k : k + 51] == CASSETTE, sid
            native = tagged[:k] + tagged[k + 51 :]
            assert sha(native) == f["native_sha256"], sid
            out.append(construct(
                "loop", f"{sid}:rescue{f['rescue_rank']}", sid.split(":")[0], f["gene"], native,
                tagged, source="rescue_round_tagged.fasta (ranker-pu-segment, held-member rescue round)",
                tier=1, order=0, reason=f"held-member rescue round, rank {f['rescue_rank']}",
                published_sha=f["tagged_sha256"],
            ))  # fmt: skip

    # --- S5E13 order rule (ranker-pu-segment, approved by the coordinator 2026-11-15): extra
    # loop arms (ranks 2-3) on internal members, "<design>:after<N>", and snorkel arms for
    # all-bad members, "<design>:snorkel" = native + the 106-aa module and cassette.
    for path in inp_all.get("order_rule_fasta", []):
        for head, tagged in fasta(path):
            f = dict(t.split("=", 1) for t in head.split()[1:] if "=" in t)
            sid = head.split()[0]
            design = sid.split(":")[0]
            if sid.endswith(":snorkel"):
                native = by_sha[f["native_sha256"]]
                assert tagged.startswith(native), sid
                out.append(construct(
                    "cterm", f"{sid}:order_rule", design, f["gene"], native, tagged,
                    source="order_rule_fold_requests.fasta (ranker-pu-segment, S5E13 order rule)",
                    tier=1, order=0, reason="all-bad member: snorkel arm (order rule)",
                    extra={"review_flags": "snorkel"}, published_sha=f["tagged_sha256"],
                ))  # fmt: skip
                continue
            k = int(sid.rsplit(":after", 1)[1])
            assert tagged[k : k + 51] == CASSETTE, sid
            native = tagged[:k] + tagged[k + 51 :]
            assert sha(native) == f["native_sha256"], sid
            out.append(construct(
                "loop", f"{sid}:order_rule_rank{f['rank'].split('.')[0]}", design, f["gene"], native,
                tagged, source="order_rule_fold_requests.fasta (ranker-pu-segment, S5E13 order rule)",
                tier=1, order=0, reason=f"extra loop arm, rank {f['rank'].split('.')[0]} (order rule)",
                published_sha=f["tagged_sha256"],
            ))  # fmt: skip

    # --- arm revisions (ranker-pu-segment, approved 2026-11-15): replacement sites forced by
    # contact rulings, and arms where v6 prefers another region; "<design>:after<N>"
    # agent-reviewer site hypotheses (2026-11-15, coordinator option b) use the same format
    # and so do the annotation-rule replacements (2026-11-15), whose headers carry no native
    # SHA-256: their native must then be exactly the design's own protein
    for path in (inp_all.get("arm_rev_fasta", []) + inp_all.get("agent_screen_fasta", [])
                 + inp_all.get("annotation_rule_fasta", [])
                 + inp_all.get("annotation_rule_b_fasta", [])):  # fmt: skip
        for head, tagged in fasta(path):
            f = dict(t.split("=", 1) for t in head.split()[1:] if "=" in t)
            sid = head.split()[0]
            k = int(sid.rsplit(":after", 1)[1])
            assert tagged[k : k + 51] == CASSETTE, sid
            native = tagged[:k] + tagged[k + 51 :]
            if "native_sha256" in f:
                assert sha(native) == f["native_sha256"], sid
            else:
                assert native == by_design[sid.split(":")[0]], sid
            out.append(construct(
                "loop", f"{sid}:{f['request']}", sid.split(":")[0], f.get("gene", ""), native, tagged,
                source=("arm_revisions_new_constructs.fasta (ranker-pu-segment, arm revisions)"
                        if "arm_revisions" in path.name
                        else f"{path.name.split('-', 1)[-1]} (ranker-pu-segment, agent site screen)"
                        if "agent_site" in path.name
                        else f"{path.name.split('-', 1)[-1]} (ranker-pu-segment, annotation rule; replaces {f.get('replaces', '')})"),
                tier=1, order=0, reason=(f"arm revision ({f['request']})" if "arm_revisions" in path.name else f["request"]), published_sha=f["tagged_sha256"],
            ))  # fmt: skip

    for r in tsv(EXP / "tag_refold_extra_requests.tsv"):
        k = int(r["insert_after_residue"]) if r["insert_after_residue"] else None
        if r.get("construct_id"):
            # a built arm sent whole (residue-1 multipass, Igk mature chain): the tagged
            # protein as sent, its native the designed protein with the stated SHA-256
            native = by_sha[r["native_sha256"]]
            out.append(construct(
                r["rung"], r["construct_id"], r["design_id"], r["gene_symbol"], native,
                r["tagged_seq"], source=f"tag_refold_extra_requests.tsv ({r['requested_by']})",
                tier=1, order=0, reason=r["reason"], published_sha=r["tagged_sha256"],
                extra={"design": r["construct_id"].split(":")[-1]},
            ))  # fmt: skip
            continue
        if r.get("tagged_seq") and r.get("native_aa"):
            # a snorkel-style construct: native, then module and cassette at the C-terminus
            tagged = r["tagged_seq"]
            native = tagged[: int(r["native_aa"])]
            assert sha(native) == r["native_sha256"], r["design_id"]
        elif r.get("tagged_seq"):
            # the requester sent the tagged protein: native = it minus the cassette, hash-checked
            tagged = r["tagged_seq"]
            assert tagged[k : k + 51] == CASSETTE, r["design_id"]
            native = tagged[:k] + tagged[k + 51 :]
            assert sha(native) == r["native_sha256"], r["design_id"]
        else:
            native = by_design[r["design_id"]]
            tagged = native[:k] + CASSETTE + native[k:]
        out.append(construct(
            r["rung"], f"{r['design_id']}:after{k}" if k is not None else f"{r['design_id']}:spacer",
            r["design_id"], r["gene_symbol"], native,
            tagged,
            source=f"tag_refold_extra_requests.tsv ({r['requested_by']})",
            tier=1, order=0, reason=r["reason"], published_sha=r["tagged_sha256"],
        ))  # fmt: skip

    # --- N after SP replacement arms for arms that failed (nterm-reviewer-package): 4 sent
    # first, then 61 (56 native-site, 5 Igk swaps). A native-site arm's native is the arm
    # minus the cassette; an Igk arm's comes from its design, and must match the arm's tail.
    for tsv_key, fa_key, label in (
        ("nterm_req_tsv", "nterm_req_fasta", "nterm_refold_requests_2026-11-12.tsv"),
        ("nterm_req_b_tsv", "nterm_req_b_fasta", "nterm_refold_requests_2026-11-12b.tsv"),
        ("nterm_req_c_tsv", "nterm_req_c_fasta", "nterm_refold_requests_2026-11-13a.tsv"),
        ("nterm_req_d_tsv", "nterm_req_d_fasta", "nterm_refold_requests_2026-11-13b.tsv"),
        ("nterm_req_e_tsv", "nterm_req_e_fasta", "nterm_refold_requests_2026-11-13c.tsv"),
        ("nterm_req_f_tsv", "nterm_req_f_fasta", "nterm_refold_requests_2026-11-14.tsv"),
        ("nterm_req_g_tsv", "nterm_req_g_fasta", "nterm_refold_requests_2026-11-15_weak_sp.tsv"),
    ):
        if tsv_key not in inp_all:
            continue
        req = fasta_of(inp_all[fa_key])
        for r in rows_of(inp_all[tsv_key]):
            tagged = req[r["tagged_sha256"]]
            c0, c1 = cassette_span(tagged)
            if r["arm"] == "native_site":
                native = tagged[:c0] + tagged[c1:]
                assert c0 == int(r["insert_after_residue"]), r["design_id"]
            else:
                native = by_design[r["design_id"]]
                assert native.endswith(tagged[c1:]), r["design_id"]
            out.append(construct(
                "nterm", r["design_id"], r["design_id"], "", native, tagged,
                source=f"{label} (nterm-reviewer-package, PR #254)",
                tier=2 if r.get("priority") == "2" else 1, order=0,
                reason=r.get("reason") or r.get("why") or (f"replacement arm for failed {r['replaces'][:12]}" if r.get("replaces") else f"Igk swap of {r.get('parent_sha256', '')[:12]}"),
                extra={"arm": r["arm"]}, published_sha=r["tagged_sha256"],
            ))  # fmt: skip
    return out


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--inputs-dir", type=Path, help="ingest rung inputs from this folder (snapshots them to S3)"
    )
    args = ap.parse_args()
    inp_all = resolve_inputs(args.inputs_dir)
    inp_paths = {role: paths[-1] for role, paths in inp_all.items()}
    constructs = build(inp_all)
    prev = {r["tagged_sha256"]: r for r in tsv(OUT)} if OUT.exists() else {}
    old_natives = (
        {r["native_sha256"]: r["native_seq"] for r in tsv(NATIVES_OUT)}
        if NATIVES_OUT.exists()
        else {}
    )
    old_tagged = (
        {r["tagged_sha256"]: r["tagged_seq"] for r in tsv(TAGGED_OUT)}
        if TAGGED_OUT.exists()
        else {}
    )
    for r in prev.values():
        if not r.get("tagged_seq"):
            r["tagged_seq"] = old_tagged[r["tagged_sha256"]]
        if not r.get("native_seq"):
            r["native_seq"] = old_natives[r["native_sha256"]]
    have_model = {r["native_sha256"] for r in tsv(NATIVE_INDEX)}

    # The order of record: every insert that will be synthesized gets a live, tier-1 row, so a
    # construct can never fall out of the verdicts because the export that first carried it
    # moved on (2026-11-13: 19 snorkel arms had only superseded rows). Matched by protein:
    # the precursor, or for an Igk design the mature chain (leader cleaved) if that is what
    # was folded. The native comes from any construct with the same tagged protein, else
    # from removing the cassette.
    if "order_fasta" in inp_paths:
        known = {c["tagged_sha256"]: c for c in constructs}
        known.update({h: r for h, r in prev.items() if h not in known})
        igk = gpcr_arms.LEADER_AA
        for head, nt in fasta(inp_paths["order_fasta"]):
            parts = head.split()
            design_id, rung = parts[-2], {"internal": "loop"}.get(parts[-1], parts[-1])
            prot = translate_insert(nt)
            # Igk multipass and type III designs are folded as the mature chain (leader cleaved);
            # every other rung, including N-after-SP Igk swaps, as the precursor. Chosen by rung,
            # not by lookup: the first queue holds never-folded precursor rows for these designs.
            mature = rung in ("multipass", "type3") and prot.startswith(igk)
            tagged = prot[len(igk) :] if mature else prot
            ref = known.get(sha(tagged))
            native = (ref or {}).get("native_seq", "")
            if not native:
                c0, c1 = cassette_span(tagged)
                native = tagged[:c0] + tagged[c1:] if not tagged.startswith(igk) else ""
            if not native:
                print(f"  order {design_id}: native not derivable; not queued", file=sys.stderr)
                continue
            constructs.append(construct(
                rung, design_id, design_id, parts[0].split("|")[1] if "|" in parts[0] else "", native, tagged,
                source="l1_order_inserts_recoded.fasta (order of record, synthesis/2026-11-13-one-arm-renamed)",
                tier=1, order=0, reason="ordered construct" + (" (mature chain folded)" if tagged != prot else ""),
            ))  # fmt: skip

    # one fold per tagged SHA-256, at the best tier any source gives it
    folds: dict[str, dict] = {}
    sources = []
    for c in sorted(constructs, key=lambda c: (c["tier"], c["order"], c["tagged_sha256"])):
        sources.append({k: c.get(k, "") for k in (
            "tagged_sha256", "rung", "source_id", "design_id", "gene_symbol", "tier", "reason",
            "arm", "design", "rank", "loop_member", "review_flags", "source")})  # fmt: skip
        f = folds.get(c["tagged_sha256"])
        if f is None:
            p = prev.get(c["tagged_sha256"])
            fid = p["fold_id"] if p else f"t.{c['tagged_sha256'][:16]}"
            folds[c["tagged_sha256"]] = {**c, "fold_id": fid, "kind": "tagged"}
    # untagged calibration controls for named members (held-member rescue round, 2026-11-14):
    # the member's designed protein, folded with the tagged arms' settings and scored against
    # the same native reference, so a member whose own untagged fold misses the bar is visible
    native_of = {
        c["design_id"]: c["native_seq"] for c in folds.values() if c.get("kind") == "tagged"
    }
    for r in tsv(EXP / "tag_refold_calibration_controls.tsv"):
        seq = native_of[r["design_id"]]
        h = sha(seq)
        # an untagged fold of this exact protein already made with these settings is reused
        folds.setdefault(h, {
            **{k: "" for k in next(iter(folds.values()))}, "fold_id": prev[h]["fold_id"] if h in prev else f"{r['design_id']}.control",
            "kind": "control", "rung": "control", "tier": 0, "order": 0,
            "gene_symbol": r["gene_symbol"], "design_id": r["design_id"], "native_seq": seq,
            "tagged_seq": seq, "native_sha256": h, "tagged_sha256": h, "native_aa": len(seq),
            "tagged_aa": len(seq), "native_mapped_aa": len(seq),
            "blocks": json.dumps([(0, 0, len(seq))]),
            "source": f"tag_refold_calibration_controls.tsv ({r['requested_by']})", "reason": r["reason"],
        })  # fmt: skip
    # natives with no model of exactly their sequence: fold untagged, at the best tier needed
    for c in list(folds.values()):
        h = c["native_sha256"]
        if h in have_model or h in folds:
            continue
        p = prev.get(h)
        fid = p["fold_id"] if p else f"n.{h[:16]}"
        folds[h] = {**{k: "" for k in c}, "fold_id": fid, "kind": "native", "rung": "native",
                    "tier": c["tier"], "order": -1, "gene_symbol": c["gene_symbol"],
                    "design_id": c["design_id"], "native_seq": c["native_seq"],
                    "tagged_seq": c["native_seq"], "native_sha256": h, "tagged_sha256": h,
                    "native_aa": c["native_aa"], "tagged_aa": c["native_aa"],
                    "blocks": json.dumps([(0, 0, len(c["native_seq"]))]),
                    "source": "native of a queued construct with no model of this exact sequence",
                    "reason": "the comparison model for its tagged forms"}  # fmt: skip
    # carry over the calibration controls and anything already folded that no rung sends any more
    for h, p in prev.items():
        if h in folds:
            continue
        if p["rung"] in ("control", "sp_control"):
            folds[h] = {**p, "kind": "control", "tier": 0}
        else:
            folds[h] = {**p, "kind": "superseded", "tier": 9,
                        "reason": "folded earlier; no longer in any rung's 2026-11-12 input"}  # fmt: skip
    cols = ["fold_id", "kind", "tier", "order", "rung", "design_id", "source_id", "gene_symbol",
            "native_aa", "tagged_aa", "native_mapped_aa", "cassette_start", "cassette_end",
            "insert_after_native_residue", "leader_end", "sp_end", "tm1_span_folded",
            "native_sha256", "tagged_sha256", "blocks", "source",
            "reason", "decided_by", "implemented_by"]  # fmt: skip
    rows = []
    for f in sorted(
        folds.values(), key=lambda f: (int(f["tier"]), float(f.get("order") or 0), f["fold_id"])
    ):
        f.setdefault("decided_by", "")
        f["decided_by"] = f["decided_by"] or f"tier and site: {DELEGATED}"
        f["implemented_by"] = "Claude Code session l1-tag-refold"
        rows.append({k: f.get(k, "") for k in cols})
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)
    uniq: dict[tuple[str, str], dict] = {}
    for src in sources:
        uniq.setdefault((src["tagged_sha256"], src["source_id"]), src)
    sources = list(uniq.values())
    with SOURCES_OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sources[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sources)
    # tagged sequences once each, in their own table (refold_queue.read_queue joins them back)
    with TAGGED_OUT.open("w", newline="") as fh:
        fh.write("tagged_sha256\ttagged_seq\n")
        for h, sq in sorted({f["tagged_sha256"]: f["tagged_seq"] for f in folds.values()}.items()):
            fh.write(f"{h}\t{sq}\n")
    # native sequences once each, in their own table, so the queue stays under 5 MB
    natives = {f["native_sha256"]: f["native_seq"] for f in folds.values() if f.get("native_seq")}
    with NATIVES_OUT.open("w", newline="") as fh:
        fh.write("native_sha256\tnative_seq\n")
        for h, sq in sorted(natives.items()):
            fh.write(f"{h}\t{sq}\n")
    from collections import Counter

    print(f"wrote {OUT.name}: {len(rows)} folds; {SOURCES_OUT.name}: {len(sources)} rung rows")
    for (k, t), n in sorted(
        Counter((r["kind"] + ":" + r["rung"], r["tier"]) for r in rows).items(),
        key=lambda x: (x[0][1], x[0][0]),
    ):
        print(f"  tier {t}  {k:<24} {n}")
    for t in (1, 2, 3):
        nat = {r["native_sha256"] for r in rows if int(r["tier"]) == t}
        print(f"  tier {t}: {sum(int(r['tier']) == t for r in rows)} folds over {len(nat)} natives")


if __name__ == "__main__":
    main()
