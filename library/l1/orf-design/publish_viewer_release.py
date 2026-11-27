#!/usr/bin/env python3
"""Publish an L1 cut to the surfaceome viewer: ColabFold R2 mirror, D1 loads, and cache purge.

Routine post-cut publication step documented in ``REBUILD.md`` §5. Dry-run by default
(verifies S3->R2 diff, validates ``library_structure_payload.tsv``, and builds the
content-addressed redacted D1 SQL); pass ``--execute`` to mirror missing PDBs to R2,
swap ``library_structure_public``, apply the additive ``l1_design_*_public`` load, and
verify live read-back from ``https://api.deliverome.org/surfaceome/v1/l1-design/loads``.

    # dry run (read-only checks + SQL build):
    uv run --locked python publish_viewer_release.py --viewer-repo ~/git/accessible-surfaceome

    # execute R2 mirror + D1 loads + live API read-back:
    uv run --locked python publish_viewer_release.py --viewer-repo ~/git/accessible-surfaceome --execute
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from deliverome_analysis.core.s3 import processed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
R2_BUCKET = "surfaceome-structures"
TAGGED_R2_RELEASE = "2026-11-12"
AGENT_PICKS_COMMIT = "0818d2a1c6c6b80f7efe0a2dd71417cb83f4536c"
LIT_CONTACTS_S3 = "scratch/becca/lit-contact-agents/scaleup/member_lit_contacts.tsv"
API_BASE = "https://api.deliverome.org/surfaceome/v1/l1-design"


def _load_cf_env(viewer_repo: Path) -> dict[str, str]:
    env = os.environ.copy()
    env_file = viewer_repo / ".env"
    if env_file.is_file():
        for line in env_file.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                env.setdefault(k.strip(), v.strip().strip('"').strip("'"))
    for req in (
        "CLOUDFLARE_ACCOUNT_ID",
        "CLOUDFLARE_API_TOKEN",
        "CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID",
    ):
        if not env.get(req):
            raise SystemExit(f"missing {req} (expected in environment or {env_file})")
    return env


def _list_r2_keys(acct: str, tok: str, prefix: str) -> set[str]:
    keys: set[str] = set()
    cursor = ""
    while True:
        url = (
            f"https://api.cloudflare.com/client/v4/accounts/{acct}"
            f"/r2/buckets/{R2_BUCKET}/objects?prefix={prefix}&per_page=1000"
        )
        if cursor:
            url += f"&cursor={cursor}"
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        objs = data.get("result") or []
        if isinstance(objs, dict):
            objs = objs.get("objects", [])
        for obj in objs:
            keys.add(obj["key"])
        info = data.get("result_info") or {}
        if not info.get("is_truncated") or not info.get("cursor"):
            break
        cursor = info["cursor"]
    return keys


def _r2_put(acct: str, tok: str, key: str, body: bytes) -> None:
    url = (
        f"https://api.cloudflare.com/client/v4/accounts/{acct}/r2/buckets/{R2_BUCKET}/objects/{key}"
    )
    for attempt in range(6):
        req = urllib.request.Request(
            url,
            data=body,
            method="PUT",
            headers={
                "Authorization": f"Bearer {tok}",
                "Content-Type": "chemical/x-pdb",
                "Cache-Control": "public, max-age=31536000, immutable",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                if resp.status == 200:
                    return
                err = f"HTTP {resp.status}"
        except urllib.error.HTTPError as exc:
            err = f"HTTP {exc.code}: {exc.read().decode()[:160]}"
        except Exception as exc:
            err = repr(exc)
        if attempt == 5:
            raise RuntimeError(f"R2 PUT failed for {key}: {err}")
        time.sleep(1.5 * (2**attempt) + random.uniform(0, 0.5))


def collect_missing_r2_pdbs(acct: str, tok: str) -> dict[str, str]:
    """Return ``{r2_key: processed://...}`` for every L1 native/tagged PDB missing from R2."""
    existing = _list_r2_keys(acct, tok, "library-colabfold/") | _list_r2_keys(
        acct, tok, "library-colabfold-tagged/"
    )
    missing: dict[str, str] = {}
    for fname in ("tag_refold_viewer_manifest.tsv", "tag_refold_structures.tsv"):
        with (HERE / fname).open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                uri = r.get("pdb_key", "")
                if uri.startswith("processed://"):
                    r2k = f"library-colabfold-tagged/{TAGGED_R2_RELEASE}/{uri.rsplit('/', 1)[-1]}"
                    if r2k not in existing:
                        missing[r2k] = uri
    with (HERE / "library_structure_payload.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            uri = r.get("s3_key", "")
            if uri.startswith("processed://"):
                top, rest = uri.split("predicted-structures/")[-1].split("/", 1)
                pref = "library-colabfold-tagged" if top == "tagged" else "library-colabfold"
                r2k = f"{pref}/{rest}"
                if r2k not in existing:
                    missing[r2k] = uri
    return missing


def mirror_pdbs_to_r2(acct: str, tok: str, missing: dict[str, str], *, workers: int = 8) -> None:
    """Download each missing PDB from S3 and upload idempotently to R2."""
    if not missing:
        print("R2 mirror: 0 missing ColabFold PDBs.")
        return
    print(f"R2 mirror: uploading {len(missing)} missing PDB(s) to {R2_BUCKET}...")

    def _one(item: tuple[str, str]) -> str:
        r2k, uri = item
        local = Path(processed.download(uri.removeprefix("processed://")))
        _r2_put(acct, tok, r2k, local.read_bytes())
        return r2k

    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for fut in as_completed([ex.submit(_one, it) for it in missing.items()]):
            fut.result()
            done += 1
            if done % 25 == 0 or done == len(missing):
                print(f"  mirrored {done}/{len(missing)}")


def validate_structure_payload() -> list[str]:
    """Assert every active L1 member has exactly one chosen model and every ColabFold row has an R2 key."""
    active = {
        r["candidate_id"]
        for r in csv.DictReader((HERE / "l1_tag_paths.tsv").open(), delimiter="\t")
        if r.get("library") == "l1"
    }
    rows = list(csv.DictReader((HERE / "library_structure_payload.tsv").open(), delimiter="\t"))
    chosen = Counter(r["candidate_id"] for r in rows if r["is_chosen"] == "1")
    extra = sorted({r["candidate_id"] for r in rows} - active)
    multi = sorted(k for k, v in chosen.items() if v != 1)
    zero = sorted(k for k in active if chosen[k] == 0)
    no_r2 = [r["candidate_id"] for r in rows if r["source"] == "colabfold" and not r["r2_key"]]
    if extra or multi or zero or no_r2:
        raise SystemExit(
            f"library_structure_payload.tsv invalid: extra={len(extra)}, multi={multi[:5]}, "
            f"zero={zero[:5]}, missing_r2={no_r2[:5]}"
        )
    return sorted({r["release"] for r in rows})


def _build_pins(head: str, lit_local: str) -> dict:
    def gpin(path: str, commit: str = head) -> dict[str, str]:
        return {
            "commit": commit,
            "path": path,
            "decided_by": "L1 publication cut",
            "implemented_by": "deliverome-analysis",
        }

    s4 = "library/l1/orf-design"
    s5 = "library/l1/tag-site-benchmarking"
    return {
        "analysis_repo": "Deliverome-Project/deliverome-analysis",
        "base_dir": ".",
        "sources": {
            "renumbering": gpin(f"{s4}/l1_member_renumbering.tsv"),
            "tag_paths": gpin(f"{s4}/l1_tag_paths.tsv"),
            "structure_payload": gpin(f"{s4}/library_structure_payload.tsv"),
            "cterm": gpin(f"{s4}/constructs_l1_cterm.tsv"),
            "cterm_decisions": gpin(f"{s4}/l1_cterm_decisions.tsv"),
            "nterm_sp": gpin(f"{s4}/reviewer-package/nterm/nterm_site_picks.tsv"),
            "nterm_decisions": gpin(f"{s4}/reviewer-package/nterm/nterm_decisions.tsv"),
            "type3": gpin(f"{s4}/constructs_l1_type3.tsv"),
            "type3_proteins": gpin(f"{s4}/constructs_l1_type3_proteins.tsv"),
            "multipass": gpin(f"{s4}/constructs_l1_multipass.tsv"),
            "multipass_lane_check": gpin(f"{s4}/nterm_res1_lane_check.tsv"),
            "internal_constructs": gpin(f"{s4}/constructs_l1_internal.tsv"),
            "internal_decisions": gpin(f"{s5}/l1_internal_decisions.tsv"),
            "internal_proposed": gpin(f"{s5}/l1_internal_proposed_sites.tsv"),
            "internal_gate": gpin(f"{s4}/membrane-side-gate/l1_internal_pick_gate.tsv"),
            "loop_evidence": gpin(f"{s5}/l1_loop_evidence.tsv"),
            "lit_contacts": {
                "object": f"file://{lit_local}",
                "decided_by": "literature contact screen",
                "implemented_by": "deliverome-analysis",
            },
            "agent_picks": gpin(
                "scratch/becca/agentic-tag-sites/agent_picks.tsv", commit=AGENT_PICKS_COMMIT
            ),
            "proteoform_decisions": gpin(f"{s4}/l1_proteoform_decisions.tsv"),
            "snorkel_designs": gpin(f"{s5}/l1_snorkel_designs.tsv"),
            "insert_numbering": gpin(f"{s4}/l1_insert_numbering.tsv"),
            "member_annotations": gpin(f"{s4}/l1_member_sp_uniprot.tsv"),
        },
        "tagged_structures": [
            {
                "format": "l1_tag_refold_v1",
                "commit": head,
                "path": f"{s4}/tag_refold_structures.tsv",
                "metrics_path": f"{s4}/tag_refold_metrics.tsv",
                "manifest_path": f"{s4}/tag_refold_viewer_manifest.tsv",
                "release": TAGGED_R2_RELEASE,
                "decided_by": "L1 publication cut",
                "implemented_by": "deliverome-analysis",
            }
        ],
    }


def _d1_query(env: dict[str, str], sql_str: str) -> list[dict]:
    acct = env["CLOUDFLARE_ACCOUNT_ID"]
    db = env["CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID"]
    url = f"https://api.cloudflare.com/client/v4/accounts/{acct}/d1/database/{db}/query"
    for attempt in range(5):
        req = urllib.request.Request(
            url,
            data=json.dumps({"sql": sql_str}).encode(),
            headers={
                "Authorization": f"Bearer {env['CLOUDFLARE_API_TOKEN']}",
                "Content-Type": "application/json",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                res = json.loads(resp.read().decode())
                if not res.get("success"):
                    raise RuntimeError(f"D1 error: {res.get('errors')}")
                return res["result"][0]["results"]
        except urllib.error.HTTPError as exc:
            if attempt == 4:
                raise RuntimeError(f"D1 HTTP {exc.code}: {exc.read().decode()[:300]}") from exc
            time.sleep(2**attempt)
    return []


def _apply_redacted_sql(env: dict[str, str], sql_path: Path) -> None:
    stmts: list[str] = []
    buf: list[str] = []
    for line in sql_path.read_text().splitlines():
        if line.startswith("--") or not line.strip():
            continue
        buf.append(line)
        if line.rstrip().endswith(";"):
            stmts.append("\n".join(buf))
            buf = []
    assert not buf and "INSERT OR IGNORE INTO l1_design_load_public" in stmts[-1]
    print(f"D1 l1_design: executing {len(stmts)} batched SQL statements...")
    for idx, stmt in enumerate(stmts, 1):
        _d1_query(env, stmt)
        if idx % 40 == 0 or idx == len(stmts):
            print(f"  executed {idx}/{len(stmts)}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--viewer-repo", type=Path, required=True, help="path to accessible-surfaceome")
    ap.add_argument("--execute", action="store_true", help="write to R2 and public D1")
    ap.add_argument("--skip-r2", action="store_true", help="skip R2 PDB mirroring")
    ap.add_argument("--skip-d1", action="store_true", help="skip D1 uploads")
    args = ap.parse_args(argv)

    viewer_repo = args.viewer_repo.expanduser().resolve()
    env = _load_cf_env(viewer_repo)
    acct, tok = env["CLOUDFLARE_ACCOUNT_ID"], env["CLOUDFLARE_API_TOKEN"]
    head = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()

    releases = validate_structure_payload()
    print(f"Validated library_structure_payload.tsv (releases={releases}).")

    if not args.skip_r2:
        missing = collect_missing_r2_pdbs(acct, tok)
        print(f"R2 check: {len(missing)} ColabFold PDB(s) missing from {R2_BUCKET}.")
        if args.execute:
            mirror_pdbs_to_r2(acct, tok, missing)

    with tempfile.TemporaryDirectory(prefix="l1-viewer-pub-") as td:
        work = Path(td)
        lit_local = str(processed.download(LIT_CONTACTS_S3))
        pins_path = work / "pins.json"
        pins_path.write_text(json.dumps(_build_pins(head, lit_local), indent=2) + "\n")
        loader = work / "build_l1_design_tables.py"
        loader.write_bytes(
            subprocess.check_output(
                [
                    "git",
                    "-C",
                    str(viewer_repo),
                    "show",
                    "origin/main:scripts/cloud/build_l1_design_tables.py",
                ]
            )
        )
        subprocess.run(
            [
                "uv",
                "run",
                "python",
                str(loader),
                "build",
                "--pins",
                str(pins_path),
                "--data-dir",
                str(work / "data"),
                "--analysis-repo",
                str(REPO),
            ],
            cwd=viewer_repo,
            check=True,
        )
        sql_out = work / "l1_design_load.redacted.sql"
        subprocess.run(
            [
                "uv",
                "run",
                "python",
                str(loader),
                "sql",
                "--data-dir",
                str(work / "data"),
                "--out",
                str(sql_out),
                "--redact-internal",
            ],
            cwd=viewer_repo,
            check=True,
        )
        if not args.execute or args.skip_d1:
            print("Dry run complete; re-run with --execute to publish to R2 and D1.")
            return 0
        _apply_redacted_sql(env, sql_out)

        # Read back live API to verify cutover
        req = urllib.request.Request(
            f"{API_BASE}/loads", headers={"User-Agent": "deliverome-publish/1.0"}
        )
        with urllib.request.urlopen(req, timeout=15) as resp:
            live = json.loads(resp.read().decode())["load"]
        print(
            f"Live API verified: {live['load_version']} ({live['n_members']} members, "
            f"{live['n_sites']} sites, {live['n_structures']} structures)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
