"""Score every tagged refold against its native model and give each a pass/flag verdict.

Reads ``tag_refold_queue.tsv`` and a local folder of folds (``<fold_id>.pdb`` +
``<fold_id>.scores.json.gz``, as ``scripts/fold_tagged_modal.py --collect-only`` writes
them). The native model is found by the SHA-256 of the native protein in
``tag_refold_native_models.tsv`` (AFDB, or a ColabFold model on S3), used only if its sequence
contains the native protein exactly; otherwise this folder's own untagged refold of the
native stands in, and with neither the construct is reported as ``no_native_model`` rather
than scored against a different protein (CLAUDE.md, "Library proteoforms").

Implemented by Claude Code session l1-tag-refold for Becca's 2026-11-11 refold request.
Every threshold is a constant below, with its reason; ``verdict_reasons`` names each test
a construct failed. A flag is a proposal for review, not a decision.

    uv run --locked python library/l1/orf-design/tag_refold/score_tagged_folds.py \\
        --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12 \\
        --out library/l1/orf-design/tag_refold_metrics.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

import structure_metrics as sm  # noqa: E402
from build_refold_inputs import ALFA, CASSETTE, SPYTAG003  # noqa: E402

CACHE = Path.home() / ".cache" / "deliverome-analysis" / "native-models"
S3_CACHE = Path.home() / ".cache" / "deliverome-analysis" / "s3" / "deliverome-processed"
COLABFOLD_PREFIX = "sprint4/s4e9-orf-design-fg-panel/predicted-structures/colabfold/"

#: Fold settings, recorded on every row (scripts/fold_tagged_modal.py).
COLABFOLD_VERSION = "1.6.2 (ghcr.io/sokrypton/colabfold:1.6.2-cuda12)"
SETTINGS = "alphafold2_ptm model_1, num_recycle=3, no relax, no templates"
MSA_SOURCE = (
    "native protein's MSA from api.colabfold.com (mmseqs2_uniref_env, filtered), cached per "
    "native SHA-256; non-native residues as gap columns"
)

# ------------------------------------------------------------------ thresholds
#: Native residues the comparison trusts: pLDDT >= 70 in the native model. Below that
#: AlphaFold itself does not commit to a position, so a different position in the tagged
#: model is not evidence the tag moved anything.
CONFIDENT_PLDDT = 70.0
#: Native-region C-alpha lDDT (confident residues) at or above which the fold is reproduced.
#: Calibrated on 17 untagged natives refolded with these same cheap settings and compared
#: with AFDB (benchmark 2026-11-12): they range 0.897-1.000, so the settings alone never
#: fall below 0.89. 0.85 leaves a margin under that floor; a tagged construct below it has
#: lost more than the settings do.
LDDT_PASS = 0.85
#: The same, over confident native residues within WINDOW of the insertion point. Across
#: 1,471 30-residue windows of those 17 controls the 1st percentile is 0.70 (median 0.99),
#: so a site window below 0.70 is rarer than 1 in 100 from the settings alone.
LOCAL_LDDT_PASS = 0.70
WINDOW = 15
#: Mean relative SASA below which an epitope is called buried. 0.25 is the usual
#: buried/exposed boundary for single residues (Tien et al. 2013 scale); for a 13-16-residue
#: peptide a mean below it means most of the epitope is packed.
EPITOPE_RSA_PASS = 0.25
#: Epitope residues (SpyTag003 or ALFA) confidently packed on the native domain: heavy
#: atoms within 4.5 Å of a native residue more than 4 positions from the junction, with
#: PAE < 10 Å both ways and the tag residue itself placed (pLDDT >= 50). Linker (GGGGS) contacts are counted and reported but never
#: flagged: a flexible linker touching the protein next to its own junction is expected.
#: Five residues is about a third of either epitope, enough to occlude a binder.
CONTACT_PAE = 10.0
CONTACT_RESIDUES_FLAG = 5
#: N-terminal rungs: does the signal peptide's hydrophobic core still read as a helix?
#: Seven residues is two helical turns, the shortest run that is a helix rather than a
#: turn; an h-region is 7-15 residues. Reported as ``sp_state`` (helix / unresolved /
#: placed_not_helical) and ``sp_note``; not part of the verdict (see below).
SP_MIN_HELIX = 7
#: A helix of 10+ residues (about the length a TM-like or signal-like helix needs to read
#: as one) inside the cassette but outside ALFA, at pLDDT >= 70, fails the construct:
#: SpyTag003 and the linkers are meant to stay flexible. A native TM1 (type III) that is one
#: helix of 15+ residues (the shortest that spans a membrane) in the native model must stay so.
CASSETTE_HELIX_FAIL = 10
#: Mean pLDDT below which AlphaFold has not placed a segment at all (its own "very low"
#: band is < 50).
PLACED_PLDDT = 50.0


def longest_run(mask: np.ndarray) -> int:
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return int(best)


_ARCHIVE: dict[str, dict[str, str]] | None = None


def archived_afdb(url: str) -> Path:
    """The archived copy of an AFDB model (``archive_afdb_models.py``), SHA-256-checked.

    AFDB is never fetched live: a model the verdicts depend on must be in
    ``tag_refold_afdb_archive.tsv`` and on ``external://alphafold-db/``, or the run stops.
    """
    global _ARCHIVE
    if _ARCHIVE is None:
        with (EXP / "tag_refold_afdb_archive.tsv").open() as fh:
            _ARCHIVE = {r["afdb_url"]: r for r in csv.DictReader(fh, delimiter="\t")}
    row = _ARCHIVE.get(url)
    if row is None:
        raise SystemExit(f"AFDB model {url} is not archived; run archive_afdb_models.py first")
    key = row["archive_key"].removeprefix("external://")
    cached = sorted(S3_CACHE.parent.glob(f"*/{key}.*"))
    if cached:
        path = cached[-1]
    else:
        from deliverome_analysis.s3 import external

        path = Path(external.download(key))
    if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise SystemExit(f"{key}: archived copy does not match its recorded SHA-256")
    return path


def native_model_for(
    native_sha: str, native: str, index: dict, folds: Path | None
) -> tuple[sm.Model | None, int, str]:
    """(model, offset of native residue 0 in the model, label), looked up by SEQUENCE.

    ``index`` is ``tag_refold_native_models.tsv`` (existing AFDB / ColabFold models by the
    SHA-256 of the protein they model). Failing that, this folder's own untagged refold of
    the native (``n.<sha16>``, or ``native.<sha12>`` from the first queue). ``folds=None``
    skips that fallback, for the controls, whose own fold is the thing being tested.
    """
    p = index.get(native_sha)
    if p and p.get("local"):
        # prefetched by main(): read the local copy, no network inside a worker
        m = sm.Model(Path(p["local"]))
        off = m.seq.find(native)
        if off < 0 and m.seq.find(native[1:]) >= 1:
            off = m.seq.find(native[1:]) - 1
        if off >= 0:
            return (
                m,
                off,
                (f"AFDB {Path(p['ref']).stem}" if p["ref"].startswith("https://") else p["ref"]),
            )
    elif p:
        try:
            ref = p["ref"]
            if ref.startswith("https://"):
                path = archived_afdb(ref)
                label = f"AFDB {Path(ref).stem}"
            else:
                from deliverome_analysis.s3 import processed

                path = processed.download(ref.removeprefix("processed://"))
                label = ref
            m = sm.Model(path)
            off = m.seq.find(native)
            if off < 0 and m.seq.find(native[1:]) >= 1:
                # non-ATG start: the model carries the literal first residue
                off = m.seq.find(native[1:]) - 1
            if off >= 0:
                return m, off, label
        except Exception as e:
            print(f"  {native_sha[:12]}: native model unavailable ({e})", file=sys.stderr)
    if folds is not None:
        for own in (folds / f"n.{native_sha[:16]}.pdb", folds / f"native.{native_sha[:12]}.pdb"):
            if own.exists():
                m = sm.Model(own)
                if m.seq == native:
                    return m, 0, f"own untagged refold {own.name}"
    return None, 0, "no_native_model"


def score(row: dict, folds: Path, index: dict) -> dict:
    fid = row["fold_id"]
    out: dict[str, object] = {k: row.get(k, "") for k in (
        "fold_id", "kind", "tier", "rung", "design_id", "source_id", "gene_symbol", "native_aa",
        "tagged_aa", "insert_after_native_residue", "tagged_sha256", "native_sha256",
    )}  # fmt: skip
    if row.get("kind") == "superseded":
        out["verdict"] = "superseded"
        out["verdict_reasons"] = row.get("reason", "")
        return out
    pdb = folds / f"{fid}.pdb"
    if not pdb.exists():
        out["verdict"] = "not_folded"
        return out
    t = sm.Model(pdb)
    if t.seq != row["tagged_seq"]:
        out["verdict"] = "error"
        out["verdict_reasons"] = "model sequence differs from the queued protein"
        return out
    scores = sm.read_scores(folds / f"{fid}.scores.json.gz")
    pae = np.array(scores["pae"], dtype=float)
    out["ptm"] = round(float(scores.get("ptm", np.nan)), 3)
    out["mean_plddt"] = round(float(t.plddt.mean()), 1)

    native = row["native_seq"]
    untagged = row.get("kind") in ("control", "native") or row["rung"] in (
        "control", "sp_control", "native_missing")  # fmt: skip
    is_control = untagged
    n_model, off, label = native_model_for(
        row["native_sha256"], native, index, None if untagged else folds
    )
    out["native_model"] = label
    if n_model is None:
        if row.get("kind") == "native":
            out["verdict"] = "native_model"
            out["verdict_reasons"] = "this fold is the comparison model for its tagged forms"
            return out
        out["verdict"] = "no_native_model"
        out["verdict_reasons"] = "no model of this exact native protein to compare against"
        return out

    tmap = [(t0 + i, n0 + i) for t0, n0, n in json.loads(row["blocks"]) for i in range(n)]
    ti = np.array([a for a, _ in tmap])
    ni = np.array([b for _, b in tmap])
    mi = ni + off
    pt, pn = t.ca[ti], n_model.ca[mi]
    conf = n_model.plddt[mi] >= CONFIDENT_PLDDT
    out["native_confident_aa"] = int(conf.sum())
    out["tm_native"] = round(sm.tm_score(pt, pn), 3)
    out["rmsd_native"] = round(sm.kabsch_rmsd(pt, pn), 2)
    if conf.sum() >= 5:
        out["tm_confident"] = round(sm.tm_score(pt[conf], pn[conf]), 3)
        out["rmsd_confident"] = round(sm.kabsch_rmsd(pt[conf], pn[conf]), 2)
        per_res = sm.lddt_ca(pt, pn, conf)
        out["lddt_confident"] = round(float(np.nanmean(per_res)), 3)
    else:
        per_res = np.full(len(ti), np.nan)
    out["dplddt_native"] = round(float(t.plddt[ti].mean() - n_model.plddt[mi].mean()), 1)

    reasons: list[str] = []
    if not is_control:
        c0, c1 = int(row["cassette_start"]), int(row["cassette_end"])
        # native index on either side of the cassette
        before = [n for tt, n in tmap if tt < c0]
        k = (max(before) + 1) if before else 0
        win = (ni >= k - WINDOW) & (ni < k + WINDOW)
        out["dplddt_site"] = (
            round(float(t.plddt[ti[win]].mean() - n_model.plddt[mi[win]].mean()), 1)
            if win.any()
            else ""
        )
        wl = per_res[win & conf]
        out["local_confident_aa"] = int((win & conf).sum())
        out["lddt_site"] = round(float(np.nanmean(wl)), 3) if np.isfinite(wl).any() else ""

        cas = np.arange(c0, c1)
        out["cassette_plddt"] = round(float(t.plddt[cas].mean()), 1)
        rsa = t.rsa()
        s0 = c0 + CASSETTE.index(SPYTAG003)
        a0 = c0 + CASSETTE.index(ALFA)
        out["spytag003_rsa"] = round(float(rsa[s0 : s0 + len(SPYTAG003)].mean()), 3)
        out["alfa_rsa"] = round(float(rsa[a0 : a0 + len(ALFA)].mean()), 3)
        out["spytag003_plddt"] = round(float(t.plddt[s0 : s0 + len(SPYTAG003)].mean()), 1)
        out["alfa_plddt"] = round(float(t.plddt[a0 : a0 + len(ALFA)].mean()), 1)

        xyz, res = t.heavy_atoms()
        far = ti[(np.abs(ti - c0) > 4) & (np.abs(ti - (c1 - 1)) > 4)]
        pairs = sm.residue_contacts(xyz, res, cas, far)
        confident = {
            (i, j)
            for i, j in pairs
            if pae[i, j] < CONTACT_PAE and pae[j, i] < CONTACT_PAE and t.plddt[i] >= PLACED_PLDDT
        }
        out["cassette_contact_res"] = len({i for i, _ in pairs})
        out["cassette_confident_contact_res"] = len({i for i, _ in confident})
        out["native_res_contacted"] = len({j for _, j in pairs})
        epitope = set(range(s0, s0 + len(SPYTAG003))) | set(range(a0, a0 + len(ALFA)))
        out["epitope_confident_contact_res"] = len({i for i, _ in confident if i in epitope})
        if confident:
            near = sorted({int(dict(tmap)[j]) + 1 for _, j in confident})
            out["contacted_native_residues"] = ",".join(map(str, near[:40]))

        lddt = out.get("lddt_confident")
        if lddt != "" and lddt is not None and float(lddt) < LDDT_PASS:
            reasons.append(f"native-region lDDT {lddt} < {LDDT_PASS}")
        if out["lddt_site"] != "" and float(out["lddt_site"]) < LOCAL_LDDT_PASS:
            reasons.append(
                f"lDDT within {WINDOW} of the site {out['lddt_site']} < {LOCAL_LDDT_PASS}"
            )
        if out["spytag003_rsa"] < EPITOPE_RSA_PASS:
            reasons.append(f"SpyTag003 buried (RSA {out['spytag003_rsa']})")
        if out["alfa_rsa"] < EPITOPE_RSA_PASS:
            reasons.append(f"ALFA buried (RSA {out['alfa_rsa']})")
        if out["epitope_confident_contact_res"] >= CONTACT_RESIDUES_FLAG:
            reasons.append(
                f"{out['epitope_confident_contact_res']} epitope residues confidently packed on the native domain"
            )

        # A helix inside the cassette other than ALFA (which is a helix by design): for a
        # type III or residue-1 tag this is the "second TM-like helix" risk.
        hel = t.helix()
        non_alfa = np.array([i for i in cas if not a0 <= i < a0 + len(ALFA)])
        run = longest_run(hel[non_alfa])
        out["cassette_longest_helix_excl_alfa"] = run
        out["spytag003_longest_helix"] = longest_run(hel[s0 : s0 + len(SPYTAG003)])
        if run >= CASSETTE_HELIX_FAIL and float(np.mean(t.plddt[non_alfa][hel[non_alfa]])) >= 70:
            reasons.append(f"confident {run}-residue helix inside the cassette (not ALFA)")
        if row.get("tm1_span_folded"):
            s1, e1 = (int(x) for x in row["tm1_span_folded"].split("-"))
            tm = np.arange(s1 - 1, e1)
            out["tm1_plddt"] = round(float(t.plddt[tm].mean()), 1)
            out["tm1_longest_helix"] = longest_run(hel[tm])
            # The same TM in the native model: fail only a helix the tag LOSES. AFDB leaves
            # some TMs unmodelled (IGSF5: no helix at pLDDT ~50), and that is not the tag.
            nat_tm = [m for tt, m in zip(ti, mi, strict=True) if s1 - 1 <= tt < e1]
            if nat_tm:
                out["tm1_longest_helix_native"] = longest_run(n_model.helix()[np.array(nat_tm)])
            if out["tm1_longest_helix"] < 15 and out.get("tm1_longest_helix_native", 0) >= 15:
                reasons.append(
                    f"native TM1 helix lost (longest helical run {out['tm1_longest_helix']} vs "
                    f"{out['tm1_longest_helix_native']} in the native model)"
                )

    # Signal peptide (native, N-after-SP rungs) or leader (Igk, tagged coordinates): is it
    # still an N-terminal helix of its own? Scored on controls too, to calibrate. AlphaFold
    # renders a segment it cannot place as an extended low-pLDDT ribbon; that is "unresolved",
    # not evidence of a lost helix, so only a CONFIDENTLY placed non-helical SP is flagged.
    ht = t.helix()
    sp_t = np.array([], dtype=int)
    if row.get("leader_end"):
        sp_t = np.arange(int(row["leader_end"]))
        out["sp_kind"] = "Igk leader"
        mature_t = ti[conf]
    elif row.get("sp_end"):
        sp = int(row["sp_end"])
        sp_t = ti[ni < sp]
        sp_n = mi[ni < sp]
        mature_t = ti[(ni >= sp + 5) & conf]
        if len(sp_t):
            out["sp_kind"] = "native signal peptide"
            out["sp_longest_helix_native"] = longest_run(n_model.helix()[sp_n])
            out["sp_plddt_native"] = round(float(n_model.plddt[sp_n].mean()), 1)
    if len(sp_t):
        xyz, res = t.heavy_atoms()
        out["sp_longest_helix"] = longest_run(ht[sp_t])
        out["sp_plddt"] = round(float(t.plddt[sp_t].mean()), 1)
        out["sp_mature_contact_res"] = len(
            {i for i, _ in sm.residue_contacts(xyz, res, sp_t, mature_t)}
        )
        if out["sp_longest_helix"] >= SP_MIN_HELIX:
            out["sp_state"] = "helix"
        elif out["sp_plddt"] < PLACED_PLDDT:
            out["sp_state"] = "unresolved"
        else:
            out["sp_state"] = "placed_not_helical"
        # Reported, never a flag: in the 7 SP controls the cheap untagged refold already
        # left 2 of 7 AFDB SP helices unresolved, so a loss against AFDB is often the
        # settings, not the tag. Cleavage is SignalP's question (nterm-reviewer-package).
        native_helix = out.get("sp_longest_helix_native", 0) >= SP_MIN_HELIX
        if out["sp_state"] != "helix" and native_helix:
            out["sp_note"] = (
                f"AFDB SP helix ({out['sp_longest_helix_native']} aa) {out['sp_state']} in this fold"
            )
    out["verdict"] = "fail" if reasons else "pass"
    if is_control:
        out["verdict"] = "control"
    out["verdict_reasons"] = "; ".join(reasons)
    out["colabfold_version"] = COLABFOLD_VERSION
    out["settings"] = SETTINGS
    out["msa_source"] = MSA_SOURCE
    out["implemented_by"] = "Claude Code session l1-tag-refold"
    out["decided_by"] = (
        ""
        if is_control
        else "verdict: l1-tag-refold (delegated by Becca 2026-11-12), thresholds in score_tagged_folds.py"
    )
    return out


def _score_one(args):
    row, folds, index = args
    try:
        return score(row, folds, index)
    except Exception as e:
        return {"fold_id": row["fold_id"], "rung": row["rung"], "tagged_sha256": row["tagged_sha256"],
                "verdict": "error", "verdict_reasons": f"{type(e).__name__}: {e}"}  # fmt: skip


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--queue", default=str(HERE.parent / "tag_refold_queue.tsv"))
    ap.add_argument("--folds", required=True)
    ap.add_argument("--out", default=str(HERE.parent / "tag_refold_metrics.tsv"))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--check",
        action="store_true",
        help="re-score from the saved models into a temporary file and require it to match --out byte for byte",
    )
    args = ap.parse_args()
    committed = Path(args.out)
    if args.check:
        import tempfile

        args.out = str(Path(tempfile.mkdtemp()) / committed.name)
    folds = Path(args.folds).expanduser()
    import importlib.util

    # this folder's refold_queue, not S5E13's module of the same name
    spec = importlib.util.spec_from_file_location(
        "tag_refold_refold_queue", HERE / "refold_queue.py"
    )
    rq = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rq)
    read_queue = rq.read_queue

    rows = read_queue(Path(args.queue))
    with (EXP / "tag_refold_native_models.tsv").open() as fh:
        index = {p["native_sha256"]: p for p in csv.DictReader(fh, delimiter="\t")}
    # Fetch native models serially first (network, politely), then score in parallel.
    # A network blip must never turn into a "no_native_model" verdict (2026-11-13: 133 rows
    # did), so every native model is fetched up front with retries, and the run stops if any
    # stays unreachable.
    import time

    from deliverome_analysis.s3 import processed

    unreachable = []
    local_of: dict[str, str] = {}
    for ref in sorted({index[r["native_sha256"]]["ref"] for r in rows
                       if r["native_sha256"] in index and r.get("kind") != "superseded"
                       and (folds / f"{r['fold_id']}.pdb").exists()}):  # fmt: skip
        for attempt in range(5):
            try:
                if ref.startswith("https://"):
                    local_of[ref] = str(archived_afdb(ref))  # stops if not archived
                    ok = True
                else:
                    # a model already in the local S3 cache (any version) is used as is: the
                    # PDB is write-once, and a HEAD per model stalled whole runs on a bad network
                    key = ref.removeprefix("processed://")
                    cached = sorted(S3_CACHE.glob(key + ".*"))
                    local_of[ref] = str(cached[-1] if cached else processed.download(key))
                    ok = True
            except Exception:
                ok = False
            if ok:
                break
            time.sleep(5 * 2**attempt)
        else:
            unreachable.append(ref)
    for p in index.values():
        if p["ref"] in local_of:
            p["local"] = local_of[p["ref"]]
    if unreachable:
        raise SystemExit(
            f"{len(unreachable)} native models unreachable, e.g. {unreachable[:3]}; rerun when the network is back"
        )
    tasks = [(r, folds, index) for r in rows]
    if args.workers <= 1:  # in-process: no worker pool to wedge (2026-11-13)
        scored = [_score_one(t) for t in tasks]
    else:
        with ProcessPoolExecutor(args.workers) as ex:
            scored = list(ex.map(_score_one, tasks, chunksize=4))
    fields: list[str] = []
    for s in scored:
        for k in s:
            if k not in fields:
                fields.append(k)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(scored)
    from collections import Counter

    print(f"wrote {args.out}")
    if args.check:
        same = Path(args.out).read_bytes() == committed.read_bytes()
        if not same:
            import difflib

            diff = list(
                difflib.unified_diff(
                    committed.read_text().splitlines(),
                    Path(args.out).read_text().splitlines(),
                    lineterm="",
                    n=0,
                )
            )
            print("\n".join(diff[:20]), file=sys.stderr)
            raise SystemExit(
                f"--check FAILED: re-scored metrics differ from {committed} ({len(diff)} diff lines)"
            )
        print(f"--check passed: re-scored metrics are byte-identical to {committed}")
    for (rung, v), n in sorted(Counter((s["rung"], s.get("verdict")) for s in scored).items()):
        print(f"  {rung:<18} {v:<16} {n}")


if __name__ == "__main__":
    main()
