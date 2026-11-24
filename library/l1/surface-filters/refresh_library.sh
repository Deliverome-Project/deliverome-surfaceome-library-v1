#!/usr/bin/env bash
#
# Refresh the S1E3 FG surface-protein library from CURRENT accessible-surfaceome
# D1 state, then ship the updated library to the viewer's DEV preview and show it
# live. Run this AFTER an upstream data change lands (e.g. the n_papers_found
# backfill and the intracellular_rescue topology run).
#
# ── PREPARE-ONLY ────────────────────────────────────────────────────────────
# This script is wiring; it is NOT auto-run. Review it, make sure the topology
# work for the intracellular_rescue cohort has landed (the guard enforces this),
# then run it yourself:
#     bash refresh_library.sh
#
# What changed that makes a refresh necessary:
#   • n_papers_found backfill  → the LL tier gates on `n_papers_found < 100`
#     (final_library.py) and previously SKIPPED the null a1_recovery records;
#     they now carry real counts, so LL membership moves.
#   • ~148 new intracellular_rescue records enter via a refreshed cohort.jsonl.
#     They are non-surface (tm=0) and only become library-evaluatable once their
#     DeepTMHMM topology (per_residue_topology + sequence) is populated — hence
#     the guard below.
#
# Overridable knobs (env or edit):
#   VIEWER      accessible-surfaceome viewer checkout. MUST be on the fg-library
#               branch (fix/s1e3-library-export) — `main` has no library UI.
#   DEV_BRANCH  Cloudflare Pages preview branch = the dev site.
#   FORCE=1     skip the topology guard (only if you know what you're doing).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIEWER="${VIEWER:-$HERE/../../../../accessible-surfaceome/viewer}"
DEV_BRANCH="${DEV_BRANCH:-dev}"
VIEWER_JSON="$VIEWER/public/data/fg-library.json"

echo "== S1E3 library refresh =="
echo "   viewer:      $VIEWER"
echo "   dev branch:  $DEV_BRANCH"
echo

# ── 0. GUARD: refuse to run until topology is complete ───────────────────────
# Aborts if ANY published deep-dive record still lacks
# deterministic_features.canonical_topology.sequence — i.e. the intracellular_
# rescue topology run hasn't finished. Authoritative check against public D1
# (cohort.jsonl deliberately drops the AA sequence, so we can't check it there).
# Needs CLOUDFLARE_API_TOKEN / CLOUDFLARE_ACCOUNT_ID / CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID
# (same .env the surfaceome tooling uses). Fail-closed: no creds -> abort.
if [ "${FORCE:-0}" != "1" ]; then
  echo "[guard] checking topology completeness in public D1 ..."
  python3 - <<'PY'
import json, os, sys, urllib.request

tok = os.environ.get("CLOUDFLARE_API_TOKEN")
acct = os.environ.get("CLOUDFLARE_ACCOUNT_ID")
db = os.environ.get("CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID")
if not (tok and acct and db):
    sys.exit("[guard] ABORT — CLOUDFLARE_* env not set; cannot verify topology. "
             "Set creds (see accessible-surfaceome .env) or re-run with FORCE=1.")

url = f"https://api.cloudflare.com/client/v4/accounts/{acct}/d1/database/{db}/query"
# canonical_topology is ENRICHED at serve time from topology_public, so the raw
# stored annotation_json.sequence stays null even when topology IS done. Check
# topology_public coverage — what the Worker enriches and what final_library reads
# via the API (fetch_cohort) — not the un-enriched stored record.
sql = ("SELECT COUNT(*) AS n FROM surface_annotation sa WHERE NOT EXISTS ("
       "SELECT 1 FROM topology_public tp WHERE tp.uniprot_acc = sa.uniprot_acc "
       "AND tp.cohort = 'human_canonical');")
req = urllib.request.Request(
    url, data=json.dumps({"sql": sql, "params": []}).encode(),
    headers={"Authorization": f"Bearer {tok}", "Content-Type": "application/json"},
)
try:
    body = json.load(urllib.request.urlopen(req, timeout=30))
except Exception as exc:  # noqa: BLE001
    sys.exit(f"[guard] ABORT — D1 query failed: {exc}")
if not body.get("success"):
    sys.exit(f"[guard] ABORT — D1 error: {body.get('errors')}")
n = body["result"][0]["results"][0]["n"]
if n:
    sys.exit(f"[guard] ABORT — {n} deep-dive gene(s) have no topology_public "
             f"(human_canonical) row; the Worker can't enrich their canonical_topology, "
             f"so the library would exclude them. Finish the topology sweep first, or "
             f"re-run with FORCE=1.")
print("[guard] OK — every deep-dive gene is enrichable from topology_public.")
PY
else
  echo "[guard] SKIPPED (FORCE=1)"
fi
echo

# ── 1. Refresh cohort from current D1 ────────────────────────────────────────
# --refresh re-fetches EVERY gene (default resume mode never re-reads an existing
# record, so it would miss the backfill). Picks up the 148 new + n_papers_found.
echo "[1/5] fetch_cohort.py --refresh"
uv run python "$HERE/fetch_cohort.py" --refresh

# ── 2. Refresh the catalog (UniProt surface bit the LL tier reads) ───────────
echo "[2/5] fetch_catalog.py"
uv run python "$HERE/fetch_catalog.py"

# ── 3. Rebuild expression joins so new genes have TPM for the expr_ok gate ───
echo "[3/5] build_expression_join.py (+ --full)"
uv run python "$HERE/build_expression_join.py"
uv run python "$HERE/build_expression_join.py" --full

# ── 4. Re-run final_library + export fg-library.json to the viewer + S3 ───────
# export_fg_library.py re-runs final_library.py (writes library_membership.tsv +
# tier_counts.json), uploads the source-of-truth to S3, and writes the viewer's
# committed fg-library.json.
echo "[4/5] export_fg_library.py --viewer-path $VIEWER_JSON"
uv run python "$HERE/export_fg_library.py" --viewer-path "$VIEWER_JSON"

# ── 5. Ship to DEV (git-based; NO wrangler deploy) ───────────────────────────
# The dev site is the `accessible-surfaceome` Cloudflare Pages project, which
# auto-builds from the accessible-surfaceome repo's `dev` branch on git push.
# There is deliberately no `wrangler pages deploy` here: step 4 already wrote the
# refreshed fg-library.json into the viewer checkout ($VIEWER_JSON). To ship it,
# commit that file and merge to `dev`; Pages rebuilds automatically. (An earlier
# version of this script wrangler-deployed to `--project-name surfaceome-viewer`,
# which is the PROD project — wrong target; removed.)
echo "[5/5] Ship to dev — git-based (no wrangler deploy)"
echo "  Refreshed viewer copy written to: $VIEWER_JSON"
echo "  Next, in the accessible-surfaceome repo:"
echo "    git add viewer/public/data/fg-library.json"
echo "    git commit -m 'chore(viewer): refresh fg-library'"
echo "    open a PR into 'dev' and merge."
echo "  The accessible-surfaceome Pages project then auto-builds dev at"
echo "    https://dev.accessible-surfaceome.pages.dev (and the zone custom"
echo "    domain, once configured). Verify the gene count there before prod."
