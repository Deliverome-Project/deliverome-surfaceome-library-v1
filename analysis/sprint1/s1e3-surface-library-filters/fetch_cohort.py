#!/usr/bin/env python3
"""Fetch the accessible-surfaceome deep-dive cohort (compact) from the public API.

Read-only, public, edge-cached endpoint. Stdlib-only (no deps) so it runs
under the system python without a venv.

Two modes, and the difference matters for any claim about how current the
analysis is:

* default (**resume**) — skips symbols already in ``cohort.jsonl``, so a re-run
  only fills gaps. Cheap, but it NEVER updates a record it already has: an
  upstream re-annotation is invisible to it.
* ``--refresh`` — re-fetches every gene, setting the previous file aside as
  ``cohort.jsonl.superseded``. Use this before any run whose output is
  presented as reflecting current D1 state.

Either way a ``cohort_provenance.json`` sidecar records when the fetch ran, in
which mode, and how many records were carried over stale, so downstream
consumers can state what they actually used instead of asserting freshness.

Keeps only the fields the displayable/isoform/ortholog filters need; drops the
big AA sequences + evidence ledger so each line stays a few KB.
"""

import argparse
import json
import os
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

BASE = "https://api.deliverome.org/surfaceome/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cohort.jsonl")
PROV = os.path.join(HERE, "cohort_provenance.json")
LOG = os.path.join(HERE, "fetch.log")
MAX_WORKERS = 8
TIMEOUT = 30


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def get_json(url, retries=6):
    delay = 1.0
    for _ in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "deliverome-analysis/surface-library"}
            )
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(delay)
                delay = min(delay * 2, 30)
                continue
            raise
        except Exception:
            time.sleep(delay)
            delay = min(delay * 2, 30)
            continue
    raise RuntimeError(f"failed after {retries}: {url}")


_ISO_KEYS = (
    "isoform_id",
    "tm_helix_count",
    "signal_peptide_length",
    "ecd_length_residues",
    "icd_length_residues",
    "n_terminal_orientation",
    "c_terminal_orientation",
    "per_residue_topology",
    "ecd_pct_identity_to_canonical",
    "full_length_pct_identity_to_canonical",
)
_ORTH_KEYS = (
    "ortholog_symbol",
    "tm_helix_count",
    "ecd_length_residues",
    "ecd_pct_identity_to_human_canonical",
    "full_length_pct_identity_to_human_canonical",
    "per_residue_topology",
    "deeptmhmm_label",
    "tm_absent_from_model",
    "n_tm_regions_absent",
    "is_canonical",
    "type",
)
_CT_KEYS = (
    "tm_helix_count",
    "signal_peptide_length",
    "ecd_length_residues",
    "icd_length_residues",
    "n_terminal_orientation",
    "c_terminal_orientation",
    "per_residue_topology",
)


def compact(r):
    g = r.get("gene", {}) or {}
    df = r.get("deterministic_features", {}) or {}
    es = r.get("executive_summary", {}) or {}
    ct = df.get("canonical_topology") or {}
    orth = df.get("orthologs") or {}
    sl = (r.get("biological_context") or {}).get("subcellular_localization") or {}
    sbd = df.get("surface_bind") or {}
    # co-receptor partners live in accessibility_risks (this block was dropped until
    # now). Carry the partner list through so the co-receptor gate reads it from the
    # cohort record directly — removes the ad-hoc coreceptor_partners.jsonl dependency.
    crr = (r.get("accessibility_risks") or {}).get("co_receptor_requirements") or {}
    # evidence_grade_summary lives in executive_summary; fold it into filters so
    # the PR137 gates can read it directly (no separate injection needed).
    filters = r.get("filters", {}) or {}
    filters["evidence_grade_summary"] = es.get("evidence_grade_summary")
    # Membership tags now served by the API (deepDiveTier() on accessible-surfaceome
    # main): carry them through so the filter reads the API's canonical/likely verdict
    # directly instead of recomputing the preset predicate locally.
    filters["deep_dive_tier"] = r.get("deep_dive_tier")
    filters["deep_dive_facet"] = r.get("deep_dive_facet")
    # count SURFACE-Bind sites whose anchor residue is EXTRACELLULAR ('O') — the
    # display-targetability signal (intracellular sites, e.g. TRAT1/OXTR, don't count).
    per_res = ct.get("per_residue_topology") or ""
    ec_sites = sum(
        1
        for s in (sbd.get("sites") or [])
        if (a := s.get("anchor_residue")) and 1 <= a <= len(per_res) and per_res[a - 1] == "O"
    )
    return {
        "gene": g.get("hgnc_symbol"),
        "uniprot": g.get("uniprot_acc"),
        "filters": filters,
        "uniprot_family": es.get("uniprot_family"),
        "hgnc_gene_groups": es.get("hgnc_gene_groups", []) or [],
        "canonical_topology": {k: ct.get(k) for k in _CT_KEYS},
        "primary_compartment": sl.get("primary_compartment"),
        "co_receptor_partners": crr.get("partners") or [],
        "surface_bind_n_sites": sbd.get("n_sites"),
        "surface_bind_n_ec_sites": ec_sites,
        "isoforms": [
            {k: t.get(k) for k in _ISO_KEYS} for t in (df.get("isoform_topologies") or [])
        ],
        "orthologs": {
            "mouse": [{k: o.get(k) for k in _ORTH_KEYS} for o in orth.get("mouse", [])],
            "cyno": [{k: o.get(k) for k in _ORTH_KEYS} for o in orth.get("cynomolgus", [])],
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--refresh",
        action="store_true",
        help=(
            "Re-fetch EVERY gene from the API instead of only the ones missing "
            "from cohort.jsonl. Use this whenever the analysis must reflect the "
            "current D1 state: the default resume mode never updates a record it "
            "already has, so re-annotations upstream are invisible to it."
        ),
    )
    args = ap.parse_args()

    symbols = [g["gene_symbol"] for g in get_json(f"{BASE}/genes")["genes"]]

    done = set()
    if args.refresh:
        # Truncate rather than append: the writer below is append-only, so
        # without this a refresh would leave the stale record in the file
        # alongside the new one and the reader would pick up whichever came
        # first.
        if os.path.exists(OUT):
            os.replace(OUT, OUT + ".superseded")
            log(f"refresh: previous cohort moved to {os.path.basename(OUT)}.superseded")
    elif os.path.exists(OUT):
        with open(OUT) as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["gene"])
                except Exception:
                    pass

    todo = [s for s in symbols if s not in done]
    mode = "refresh" if args.refresh else "resume"
    log(f"mode={mode} cohort={len(symbols)} already_done={len(done)} todo={len(todo)}")
    if not args.refresh and done:
        log(
            f"NOTE: resume mode reuses {len(done)} existing records as-is. "
            f"Upstream re-annotations will NOT be picked up — re-run with "
            f"--refresh before any analysis that claims to be current."
        )

    lock = threading.Lock()
    counters = {"ok": 0, "fail": 0}

    def work(sym):
        try:
            c = compact(get_json(f"{BASE}/genes/{sym}"))
            with lock:
                with open(OUT, "a") as fh:
                    fh.write(json.dumps(c) + "\n")
                counters["ok"] += 1
                n = counters["ok"] + counters["fail"]
                if n % 250 == 0:
                    log(f"progress ok={counters['ok']} fail={counters['fail']} ({n}/{len(todo)})")
        except Exception as e:
            with lock:
                counters["fail"] += 1
                with open(LOG, "a") as fh:
                    fh.write(f"FAIL {sym}: {e}\n")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        list(ex.map(work, todo))
    total = len(done) + counters["ok"]
    log(f"DONE ok={counters['ok']} fail={counters['fail']} total_in_file={total}")

    # Provenance sidecar. Downstream (export_fg_library.py) stamps this into
    # its output instead of asserting the cohort is current — a claim no
    # consumer could previously check. `stale_records` is the number carried
    # over untouched from a previous run; it is 0 only after a --refresh.
    prov = {
        "fetched_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "api_base": BASE,
        "mode": mode,
        "n_genes_in_catalog": len(symbols),
        "n_fetched_this_run": counters["ok"],
        "n_failed_this_run": counters["fail"],
        "stale_records": len(done),
        "total_records": total,
    }
    with open(PROV, "w") as fh:
        json.dump(prov, fh, indent=2)
    log(f"provenance written to {os.path.basename(PROV)} (stale_records={len(done)})")


if __name__ == "__main__":
    main()
