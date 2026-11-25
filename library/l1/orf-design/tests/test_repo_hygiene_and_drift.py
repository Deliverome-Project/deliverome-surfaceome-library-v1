"""Offline single-source-of-truth and anti-drift guardrails for L1.

Runs 100% offline with zero network or S3 calls in both the source repository
(``deliverome-analysis``) and an exported clean-path review repository
(``deliverome-surfaceome-library-v1`` with ``EXPORT_MANIFEST.json``).

Guards five recurring drift root causes:
1. Single source of truth for the 150-member isoform panel and architecture table.
2. Active ``read_pinned_*`` commits and paths stay synchronized with
   ``export_review_repo.PINNED_GIT_SNAPSHOTS`` and ``historical/git-snapshots/``.
3. Zero broken relative markdown links in reviewer-facing documentation.
4. Reviewer READMEs match live ``bundle_manifest.tsv``, ``data_locations.tsv``,
   and ``library_pipeline.STAGES`` counts.
5. No fail-open regressions in exported core modules, and ``bundle.py audit`` passes.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest
from botocore.exceptions import ClientError

S4E9 = Path(__file__).resolve().parents[1]
RP = S4E9 / "reviewer-package"
REPO = S4E9.parents[2]
S5E13 = REPO / "library/l1/tag-site-benchmarking"
S1E3 = REPO / "library/l1/surface-filters"

if str(S4E9) not in sys.path:
    sys.path.insert(0, str(S4E9))
if str(RP) not in sys.path:
    sys.path.insert(0, str(RP))

import export_review_repo as err  # noqa: E402
import library_pipeline  # noqa: E402
import pinned_refs  # noqa: E402

SHA_RE = re.compile(r"^[0-9a-f]{7,40}$")
MD_LINK_RE = re.compile(r'!?\[[^\]]*\]\(([^)\s]+)(?:\s+"[^"]*")?\)')
STAGE_ROW_RE = re.compile(r"^\|\s*\d+\s*\|\s*`([a-z0-9_]+)`\s*\|", re.M)
NON_COMMIT_PIN_NAMES = {"START_TRIMS_COMMIT"}
EXPORTED_REPO_ONLY_ROOT_FILES = {"LICENSE", "LICENSING.md", "CITATION.cff"}
WRAPPER_FN_NAMES = {"git_show", "show", "at", "_at_main"}


def _is_exported_review_repo() -> bool:
    return (REPO / err.MANIFEST).is_file()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _isoform_key(row: dict[str, str]) -> str:
    uniprot = row["uniprot"].strip()
    if uniprot:
        return f"iso:{row['gene']}:{uniprot}"
    first_tx = row["ensembl_transcripts"].replace(",", ";").split(";")[0].strip()
    return f"iso:{row['gene']}:{first_tx}"


def _norm_rel(path: str) -> str:
    return err._remap_path(path) if _is_exported_review_repo() else path


def _strip_markdown_code(text: str) -> str:
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    return re.sub(r"`[^`\n]+`", "", text)


def _eval_str(node: ast.AST, env: dict[str, str]) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name) and node.id in env:
        return env[node.id]
    if isinstance(node, ast.JoinedStr):
        parts: list[str] = []
        for item in node.values:
            if isinstance(item, ast.Constant) and isinstance(item.value, str):
                parts.append(item.value)
            elif isinstance(item, ast.FormattedValue):
                resolved = _eval_str(item.value, env)
                if resolved is None:
                    return None
                parts.append(resolved)
            else:
                return None
        return "".join(parts)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _eval_str(node.left, env)
        right = _eval_str(node.right, env)
        if left is not None and right is not None:
            return left + right
    return None


def _module_str_env(tree: ast.Module) -> dict[str, str]:
    env: dict[str, str] = {}
    for node in tree.body:
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
        ):
            val = _eval_str(node.value, env)
            if val is not None:
                env[node.targets[0].id] = val
        elif (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.value is not None
        ):
            val = _eval_str(node.value, env)
            if val is not None:
                env[node.target.id] = val
    return env


def _collect_local_str_env(
    fn: ast.FunctionDef | ast.AsyncFunctionDef, base: dict[str, str]
) -> dict[str, str]:
    env = dict(base)
    for node in ast.walk(fn):
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
        ):
            val = _eval_str(node.value, env)
            if val is not None:
                env[node.targets[0].id] = val
    return env


def _wrapper_specs(
    tree: ast.Module, env: dict[str, str]
) -> dict[str, tuple[str | None, str | None, int | None, int | None]]:
    """Map local read_pinned_* wrappers -> (default_commit, path_prefix, commit_idx, path_idx)."""
    specs: dict[str, tuple[str | None, str | None, int | None, int | None]] = {}
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef) or node.name not in WRAPPER_FN_NAMES:
            continue
        params = [a.arg for a in node.args.args]
        defaults: dict[str, str] = {}
        offset = len(params) - len(node.args.defaults)
        for p_name, d_node in zip(params[offset:], node.args.defaults, strict=True):
            val = _eval_str(d_node, env)
            if val is not None:
                defaults[p_name] = val
        for sub in ast.walk(node):
            if not isinstance(sub, ast.Call):
                continue
            fn = sub.func
            call_name = (
                fn.id
                if isinstance(fn, ast.Name)
                else (fn.attr if isinstance(fn, ast.Attribute) else "")
            )
            if call_name not in ("read_pinned_text", "read_pinned_bytes") or len(sub.args) < 2:
                continue
            c_arg, p_arg = sub.args[0], sub.args[1]
            c_idx = (
                params.index(c_arg.id)
                if isinstance(c_arg, ast.Name) and c_arg.id in params
                else None
            )
            default_c = defaults.get(c_arg.id) if isinstance(c_arg, ast.Name) else None
            if default_c is None:
                default_c = _eval_str(c_arg, env)
            p_idx = None
            prefix: str | None = ""
            if isinstance(p_arg, ast.Name) and p_arg.id in params:
                p_idx = params.index(p_arg.id)
            elif isinstance(p_arg, ast.JoinedStr) and p_arg.values:
                last = p_arg.values[-1]
                if (
                    isinstance(last, ast.FormattedValue)
                    and isinstance(last.value, ast.Name)
                    and last.value.id in params
                ):
                    p_idx = params.index(last.value.id)
                    prefix = _eval_str(ast.JoinedStr(values=p_arg.values[:-1]), env)
            specs[node.name] = (default_c, prefix, c_idx, p_idx)
            break
    return specs


def _extract_active_pins(path: Path, rel_to_s4e9: str) -> tuple[set[str], set[tuple[str, str]]]:
    """Return (active_commits, resolved_pairs) referenced by ``read_pinned_*`` in ``path``."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if "read_pinned_" not in text:
        return set(), set()
    tree = ast.parse(text, filename=str(path))
    env = _module_str_env(tree)
    hist_names = {n for r, n in pinned_refs.HISTORICAL if r == rel_to_s4e9}
    wrappers = _wrapper_specs(tree, env)

    active_commits: set[str] = set()
    pairs: set[tuple[str, str]] = set()

    def _visit_call(call: ast.Call, local_env: dict[str, str]) -> None:
        fn = call.func
        name = (
            fn.id
            if isinstance(fn, ast.Name)
            else (fn.attr if isinstance(fn, ast.Attribute) else "")
        )
        if name in ("read_pinned_text", "read_pinned_bytes") and len(call.args) >= 2:
            c_arg, p_arg = call.args[0], call.args[1]
            if isinstance(c_arg, ast.Name) and c_arg.id in hist_names:
                return
            c = _eval_str(c_arg, local_env)
            p = _eval_str(p_arg, local_env)
            if c and SHA_RE.match(c):
                active_commits.add(c)
                if p:
                    pairs.add((c, p))
        elif name in wrappers:
            default_c, prefix, c_idx, p_idx = wrappers[name]
            c_arg = call.args[c_idx] if c_idx is not None and c_idx < len(call.args) else None
            if isinstance(c_arg, ast.Name) and c_arg.id in hist_names:
                return
            c = _eval_str(c_arg, local_env) if c_arg is not None else default_c
            p_arg = call.args[p_idx] if p_idx is not None and p_idx < len(call.args) else None
            p_tail = _eval_str(p_arg, local_env) if p_arg is not None else None
            if c and SHA_RE.match(c):
                active_commits.add(c)
                if prefix is not None and p_tail is not None:
                    pairs.add((c, f"{prefix}{p_tail}"))
        for kw in call.keywords:
            if kw.arg and ("commit" in kw.arg or kw.arg in ("ref", "rev")):
                val = _eval_str(kw.value, local_env)
                if val and SHA_RE.match(val):
                    active_commits.add(val)

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            local_env = _collect_local_str_env(node, env)
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call):
                    _visit_call(sub, local_env)
        else:
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call):
                    _visit_call(sub, env)

    referenced = {
        n.id
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        for n in ast.walk(node)
        if isinstance(n, ast.Name)
    }
    for key, val in env.items():
        if key in hist_names or key in NON_COMMIT_PIN_NAMES or "SHA256" in key:
            continue
        if not SHA_RE.match(val):
            continue
        if (key.endswith(("_REF", "_COMMIT", "_SHA")) or key == "BEFORE") and key in referenced:
            active_commits.add(val)

    return active_commits, pairs


def _exported_python_files() -> list[tuple[str, Path]]:
    """Return ``(repo_rel_posix, abs_path)`` for every exported Python script."""
    skip_names = {"provenance.py", "export_review_repo.py", "pinned_refs.py"}
    if _is_exported_review_repo():
        manifest = json.loads((REPO / err.MANIFEST).read_text(encoding="utf-8"))["files"]
        rels = sorted(p for p in manifest if p.endswith(".py"))
    else:
        collected = set(err.collect("HEAD", include_historical=False))
        for folder in (S4E9, S1E3):
            for p in folder.rglob("*.py"):
                collected.add(p.relative_to(REPO).as_posix())
        rels = sorted(p for p in collected if p.endswith(".py"))

    out: list[tuple[str, Path]] = []
    for rel in rels:
        parts = Path(rel).parts
        if (
            "tests" in parts
            or "historical" in parts
            or "review-repo" in parts
            or "data" in parts
            or parts[-1] in skip_names
        ):
            continue
        abs_path = REPO / rel
        if abs_path.is_file():
            out.append((rel, abs_path))
    return out


# ---------------------------------------------------------------------------
# 1. Single source of truth for isoform panel & architecture classes
# ---------------------------------------------------------------------------


def test_isoform_panel_matches_active_manifest_1_to_1() -> None:
    """isoform_panel.tsv (150 rows) must match active L1 isoforms in candidate_manifest_tranches.tsv."""
    iso_rows = _read_tsv(S4E9 / "isoform_panel.tsv")
    manifest_rows = _read_tsv(S4E9 / "candidate_manifest_tranches.tsv")
    active_iso = [
        r
        for r in manifest_rows
        if r["set"] == "isoform" and r["library"] == "l1" and not r["retired_reason"].strip()
    ]

    assert len(iso_rows) == 150, f"expected 150 rows in isoform_panel.tsv, got {len(iso_rows)}"
    assert len(active_iso) == 150, (
        f"expected 150 active L1 isoforms in candidate_manifest_tranches.tsv, got {len(active_iso)}"
    )

    iso_by_key = {_isoform_key(r): r for r in iso_rows}
    manifest_by_key = {r["key"]: r for r in active_iso}
    assert len(iso_by_key) == 150, "duplicate keys in isoform_panel.tsv"
    assert len(manifest_by_key) == 150, "duplicate keys in active L1 isoform manifest"

    only_panel = sorted(set(iso_by_key) - set(manifest_by_key))
    only_manifest = sorted(set(manifest_by_key) - set(iso_by_key))
    assert not only_panel and not only_manifest, (
        f"isoform_panel.tsv vs candidate_manifest_tranches.tsv drift: "
        f"only_in_panel={only_panel[:10]}, only_in_manifest={only_manifest[:10]}"
    )

    tx_mismatches: list[tuple[str, str, set[str]]] = []
    for key, m_row in sorted(manifest_by_key.items()):
        p_row = iso_by_key[key]
        panel_txs = {
            t.strip()
            for t in p_row["ensembl_transcripts"].replace(",", ";").split(";")
            if t.strip()
        }
        if m_row["transcript"].strip() not in panel_txs:
            tx_mismatches.append((key, m_row["transcript"], panel_txs))
    assert not tx_mismatches, (
        f"manifest transcripts absent from isoform_panel.tsv: {tx_mismatches[:5]}"
    )


def test_isoform_architecture_classes_matches_panel_and_gate() -> None:
    """isoform_architecture_classes.tsv must sum to 150 panel forms and 1,749 passing gate forms."""
    iso_rows = _read_tsv(S4E9 / "isoform_panel.tsv")
    iac_rows = _read_tsv(S4E9 / "isoform_architecture_classes.tsv")
    gate_rows = _read_tsv(S4E9 / "isoform_gate_2026_11_23.tsv")

    panel_sum = sum(int(r["panel_forms"]) for r in iac_rows)
    assert panel_sum == len(iso_rows) == 150, (
        f"panel_forms sum ({panel_sum}) != len(isoform_panel.tsv) ({len(iso_rows)})"
    )

    panel_by_class = Counter(r["topology_class"] for r in iso_rows)
    iac_by_class = {r["class"]: int(r["panel_forms"]) for r in iac_rows if int(r["panel_forms"])}
    assert iac_by_class == dict(panel_by_class), (
        f"per-class panel_forms drift: table={iac_by_class}, panel={dict(panel_by_class)}"
    )

    panel_tiers = Counter(r["panel_tier"] for r in iso_rows)
    assert sum(int(r["panel_forms_observed"]) for r in iac_rows) == panel_tiers["observed"]
    assert (
        sum(int(r["panel_forms_prediction_only"]) for r in iac_rows)
        == panel_tiers["prediction_only"]
    )

    drawable_sum = sum(int(r["drawable_forms_gate_of_record"]) for r in iac_rows)
    gate_pass = sum(1 for r in gate_rows if r["verdict"] == "pass")
    assert drawable_sum == gate_pass == 1749, (
        f"drawable_forms_gate_of_record sum ({drawable_sum}) != gate pass rows ({gate_pass})"
    )


# ---------------------------------------------------------------------------
# 2. Zero stale PINNED_GIT_SNAPSHOTS vs active read_pinned_* calls
# ---------------------------------------------------------------------------


def test_pinned_git_snapshots_cover_all_active_read_pinned_calls() -> None:
    """Every active read_pinned_* commit and path is covered by PINNED_GIT_SNAPSHOTS
    and (in an exported repo) materialized in historical/git-snapshots/SHA256SUMS.json."""
    s4e9_prefix = S4E9.relative_to(REPO).as_posix() + "/"
    active_commits: set[str] = set()
    active_pairs: set[tuple[str, str]] = set()

    for rel, abs_path in _exported_python_files():
        rel_to_s4 = rel[len(s4e9_prefix) :] if rel.startswith(s4e9_prefix) else rel
        commits, pairs = _extract_active_pins(abs_path, rel_to_s4)
        active_commits |= commits
        active_pairs |= {(c, _norm_rel(p)) for c, p in pairs}

    assert active_commits, (
        "expected to discover active read_pinned_* commits across exported scripts"
    )

    declared_pairs = {(c, _norm_rel(p)) for c, p in err.PINNED_GIT_SNAPSHOTS}
    declared_commits = {c for c, _ in declared_pairs}

    missing_commits = sorted(active_commits - declared_commits)
    assert not missing_commits, (
        f"active read_pinned_* commits missing from export_review_repo.PINNED_GIT_SNAPSHOTS: "
        f"{missing_commits}"
    )

    missing_pairs = sorted(active_pairs - declared_pairs)
    assert not missing_pairs, (
        f"active read_pinned_* (commit, path) pairs missing from PINNED_GIT_SNAPSHOTS: "
        f"{missing_pairs}"
    )

    if _is_exported_review_repo():
        snap_root = REPO / "historical" / "git-snapshots"
        manifest_path = snap_root / "SHA256SUMS.json"
        assert manifest_path.is_file(), f"missing {manifest_path.relative_to(REPO)}"
        manifest: dict[str, str] = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest_commits = {k.split("/", 1)[0] for k in manifest}
        missing_exp_commits = sorted(active_commits - manifest_commits)
        assert not missing_exp_commits, (
            f"active read_pinned_* commits missing from SHA256SUMS.json: {missing_exp_commits}"
        )

        missing_snaps: list[str] = []
        bad_hashes: list[str] = []
        for commit, rel in sorted(declared_pairs | active_pairs):
            key = f"{commit}/{rel}"
            target = snap_root / key
            if key not in manifest or not target.is_file():
                missing_snaps.append(key)
            elif hashlib.sha256(target.read_bytes()).hexdigest() != manifest[key]:
                bad_hashes.append(key)
        assert not missing_snaps, f"unmaterialized pinned git snapshots: {missing_snaps[:10]}"
        assert not bad_hashes, f"checksum mismatch in pinned git snapshots: {bad_hashes[:10]}"


# ---------------------------------------------------------------------------
# 3. Zero broken relative markdown links
# ---------------------------------------------------------------------------


def test_markdown_docs_have_no_broken_relative_links() -> None:
    """Every relative markdown link in the reviewer docs must resolve to an existing path."""
    targets: list[tuple[Path, Path, set[str]]] = [
        (RP / "BRIEF.md", RP, set()),
        (RP / "README.md", RP, set()),
        (RP / "COVER_NOTE.md", RP, set()),
        (RP / "DESIGN.md", RP, set()),
        (S4E9 / "REBUILD.md", S4E9, set()),
    ]
    if _is_exported_review_repo():
        targets.append((REPO / "README.md", REPO, set()))
        if (REPO / "LICENSING.md").is_file():
            targets.append((REPO / "LICENSING.md", REPO, set()))
    else:
        # In the source repo, review-repo/README.md is the template for the exported root README.
        targets.append((RP / "review-repo" / "README.md", REPO, EXPORTED_REPO_ONLY_ROOT_FILES))

    broken: list[str] = []
    for doc_path, base_dir, allowed_missing in targets:
        assert doc_path.is_file(), f"expected documentation file missing: {doc_path}"
        cleaned = _strip_markdown_code(doc_path.read_text(encoding="utf-8"))
        for match in MD_LINK_RE.finditer(cleaned):
            raw = match.group(1).strip().strip("<>")
            if raw.startswith(("http://", "https://", "mailto:", "#")):
                continue
            rel_target = raw.split("#", 1)[0].split("?", 1)[0]
            if not rel_target or rel_target in allowed_missing:
                continue
            resolved = (base_dir / rel_target).resolve()
            if not resolved.exists():
                broken.append(f"{doc_path.relative_to(REPO)}: {raw}")

    assert not broken, "broken relative markdown links:\n  " + "\n  ".join(broken[:20])


# ---------------------------------------------------------------------------
# 4. README vs bundle_manifest.tsv / data_locations.tsv / library_pipeline.py
# ---------------------------------------------------------------------------


def test_readmes_match_manifest_locations_and_pipeline_counts() -> None:
    """Reviewer READMEs must match exact manifest, data-location, and pipeline stage counts."""
    manifest_rows = _read_tsv(RP / "bundle_manifest.tsv")
    location_rows = _read_tsv(RP / "data_locations.tsv")

    total = len(manifest_rows)
    replay_yes = sum(1 for r in manifest_rows if r["read_by_replay"] == "yes")
    replay_no = sum(1 for r in manifest_rows if r["read_by_replay"] == "no")
    homes = Counter(r["home"] for r in location_rows)
    stage_ids = [s.id for s in library_pipeline.STAGES]
    n_stages = len(stage_ids)

    assert replay_yes + replay_no == total == len(location_rows)

    root_readme = (
        (REPO / "README.md") if _is_exported_review_repo() else (RP / "review-repo" / "README.md")
    )
    readmes = [RP / "README.md", root_readme]

    stale_patterns = (
        re.compile(r"\b3,504\b"),
        re.compile(r"\b3,240\b"),
        re.compile(r"\b16\s+(?:pipeline\s+)?stages\b"),
        re.compile(r"\b23\s+(?:pipeline\s+)?stages\b"),
    )
    stage_count_re = re.compile(rf"\b{n_stages}(?:\*\*)?\s+(?:pipeline\s+)?stages\b")

    for readme in readmes:
        text = readme.read_text(encoding="utf-8")
        rel = readme.relative_to(REPO)
        for label, count in (
            ("bundle total", total),
            ("read_by_replay=yes", replay_yes),
            ("read_by_replay=no", replay_no),
            ("public_r2", homes["public_r2"]),
            ("public_source", homes["public_source"]),
            ("review_repo", homes["review_repo"]),
            ("private_release", homes["private_release"]),
        ):
            token = f"{count:,}"
            assert re.search(rf"(?<!\d){re.escape(token)}(?!\d)", text), (
                f"{rel} is missing exact {label} count {token}"
            )
        assert stage_count_re.search(text), f"{rel} does not state {n_stages} pipeline stages"
        for stale in stale_patterns:
            assert not stale.search(text), f"{rel} contains stale drift pattern {stale.pattern!r}"

    rp_text = (RP / "README.md").read_text(encoding="utf-8")
    documented_stages = STAGE_ROW_RE.findall(rp_text)
    assert documented_stages == stage_ids, (
        f"reviewer-package/README.md stage table ({len(documented_stages)}) != "
        f"library_pipeline.STAGES ({len(stage_ids)})"
    )


# ---------------------------------------------------------------------------
# 5. No fail-open regressions in core modules & bundle.py audit
# ---------------------------------------------------------------------------


def test_s3_exists_is_fail_closed_on_non_not_found_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """s3.py exists() must check 404/NoSuchKey and re-raise non-404 errors."""
    s3_path = REPO / "src" / "deliverome_analysis" / "s3.py"
    tree = ast.parse(s3_path.read_text(encoding="utf-8"), filename=str(s3_path))

    exists_fns = [
        n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "exists"
    ]
    assert exists_fns, "s3.py missing exists() method"
    exists_fn = exists_fns[0]
    handlers = [n for n in ast.walk(exists_fn) if isinstance(n, ast.ExceptHandler)]
    assert handlers, "s3.py exists() missing try/except handler"
    assert any(isinstance(n, ast.Raise) for h in handlers for n in ast.walk(h)), (
        "s3.py exists() must re-raise non-404 ClientError instead of returning False"
    )
    handler_calls = {
        (n.func.id if isinstance(n.func, ast.Name) else n.func.attr)
        for h in handlers
        for n in ast.walk(h)
        if isinstance(n, ast.Call) and isinstance(n.func, (ast.Name, ast.Attribute))
    }
    handler_consts = {
        n.value
        for h in handlers
        for n in ast.walk(h)
        if isinstance(n, ast.Constant) and isinstance(n.value, str)
    }
    assert "_is_not_found" in handler_calls or {"404", "NoSuchKey"} <= handler_consts, (
        "s3.py exists() must gate False on 404 / NoSuchKey / _is_not_found"
    )

    helper_fns = [
        n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "_is_not_found"
    ]
    assert helper_fns, "s3.py missing _is_not_found() helper"
    helper_consts = {
        n.value
        for n in ast.walk(helper_fns[0])
        if isinstance(n, ast.Constant) and isinstance(n.value, str)
    }
    assert {"404", "NoSuchKey"} <= helper_consts

    from deliverome_analysis import s3 as s3_mod

    monkeypatch.delenv("DELIVEROME_S3_MIRROR", raising=False)

    class _FakeClient:
        def __init__(self, code: str, status: int) -> None:
            self._code = code
            self._status = status

        def head_object(self, **_kwargs: object) -> dict[str, object]:
            raise ClientError(
                {
                    "Error": {"Code": self._code},
                    "ResponseMetadata": {"HTTPStatusCode": self._status},
                },
                "HeadObject",
            )

    bucket_404 = s3_mod._ReadBucket("test-bucket", lambda: _FakeClient("404", 404))
    assert bucket_404.exists("missing.tsv") is False

    bucket_nsk = s3_mod._ReadBucket("test-bucket", lambda: _FakeClient("NoSuchKey", 404))
    assert bucket_nsk.exists("missing.tsv") is False

    bucket_403 = s3_mod._ReadBucket("test-bucket", lambda: _FakeClient("403", 403))
    with pytest.raises(ClientError):
        bucket_403.exists("forbidden.tsv")


def test_membrane_side_site_verdict_rejects_non_finite_inputs() -> None:
    """membrane_side.py site_verdict() must check isfinite and reject NaN/Inf."""
    ms_path = REPO / "src" / "deliverome_analysis" / "membrane_side.py"
    tree = ast.parse(ms_path.read_text(encoding="utf-8"), filename=str(ms_path))

    fns = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "site_verdict"]
    assert fns, "membrane_side.py missing site_verdict()"
    called = {
        (n.func.id if isinstance(n.func, ast.Name) else n.func.attr)
        for n in ast.walk(fns[0])
        if isinstance(n, ast.Call) and isinstance(n.func, (ast.Name, ast.Attribute))
    }
    assert "isfinite" in called, "membrane_side.site_verdict() must check isfinite"

    from deliverome_analysis import membrane_side as ms

    th = ms.GateThresholds(
        outside=12.0,
        margin=4.0,
        core=10.0,
        min_plddt=50.0,
        min_topology_ratio=2.0,
        max_kr_excess_outside=2,
    )
    assert ms.site_verdict(float("nan"), 90.0, True, th) == "unscorable"
    assert ms.site_verdict(20.0, float("nan"), True, th) == "unscorable"
    assert ms.site_verdict(float("inf"), 90.0, True, th) == "unscorable"
    assert ms.site_verdict(float("-inf"), 90.0, True, th) == "unscorable"


def test_v8_prodomino_fail_closed_and_rescue_round_gpi_clearance() -> None:
    """v8 prodomino must fail closed on usable == 0; rescue_round must keep GPI clearance >= 12."""
    v8_path = S4E9 / "design_l1_internal_v8_prodomino.py"
    v8_tree = ast.parse(v8_path.read_text(encoding="utf-8"), filename=str(v8_path))
    mask_fns = [
        n
        for n in ast.walk(v8_tree)
        if isinstance(n, ast.FunctionDef) and n.name == "apply_membrane_mask"
    ]
    assert mask_fns, "design_l1_internal_v8_prodomino.py missing apply_membrane_mask()"
    mask_src = ast.unparse(mask_fns[0])
    assert "unusable_shas" in mask_src and "eligible.append(False)" in mask_src, (
        "apply_membrane_mask() must explicitly fail closed on unusable rows"
    )

    import design_l1_internal_v8_prodomino as v8

    pool = pd.DataFrame(
        {
            "sequence_sha256": ["sha_unusable", "sha_unsupp", "sha_ok"],
            "residue": [2, 2, 2],
            "region_index": [0, 0, 0],
        }
    )
    mm = pd.DataFrame(
        {
            "sequence_sha256": ["sha_unusable", "sha_unsupp", "sha_ok"],
            "usable": [0, 1, 1],
            "orientation_supported": [True, False, True],
            "mask": ["OOOO", "OOOO", "OOOO"],
            "junction_mask": ["EEEE", "EEEE", "EEEE"],
        }
    )
    out = v8.apply_membrane_mask(pool, mm)
    assert list(out["eligible"]) == [False, False, True], (
        f"usable==0 or orientation_supported==False must not be eligible: {list(out['eligible'])}"
    )

    rescue_path = S5E13 / "rescue_round.py"
    assert rescue_path.is_file(), f"missing {rescue_path}"
    r_tree = ast.parse(rescue_path.read_text(encoding="utf-8"), filename=str(rescue_path))
    clearance: int | None = None
    for node in r_tree.body:
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "GPI_OMEGA_UPSTREAM_CLEARANCE"
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, int)
        ):
            clearance = node.value.value
            break
    assert clearance is not None and clearance >= 12, (
        f"rescue_round.py GPI_OMEGA_UPSTREAM_CLEARANCE must be >= 12, got {clearance!r}"
    )

    uses_clearance = False
    for node in ast.walk(r_tree):
        if isinstance(node, ast.Compare) and any(isinstance(op, ast.LtE) for op in node.ops):
            cmp_src = ast.unparse(node)
            if "omega_member - GPI_OMEGA_UPSTREAM_CLEARANCE" in cmp_src:
                uses_clearance = True
                break
    assert uses_clearance, (
        "rescue_round.py must enforce r <= omega_member - GPI_OMEGA_UPSTREAM_CLEARANCE"
    )


def test_bundle_audit_exits_zero_offline() -> None:
    """bundle.py audit must exit 0 offline with zero missing, drifted, or empty prefixes."""
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src")
    env["S3_ENDPOINT_URL"] = "http://127.0.0.1:9"
    res = subprocess.run(
        [sys.executable, str(RP / "bundle.py"), "audit"],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert res.returncode == 0, f"bundle.py audit failed:\n{res.stdout}\n{res.stderr}"
    for required in (
        "pinned but not bundled: 0",
        "SHA-256 differs: 0",
        "release prefix with no bundled file: 0",
    ):
        assert required in res.stdout, f"bundle.py audit output missing {required!r}:\n{res.stdout}"


# ---------------------------------------------------------------------------
# 6. Multi-file order-of-record pin synchronization
# ---------------------------------------------------------------------------


def _str_constants_from_file(path: Path) -> dict[str, str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    return _module_str_env(tree)


def test_active_order_pins_agree_across_all_files() -> None:
    """Active order-of-record refs, FASTA keys, SHA-256s, and insert IDs must agree 100%."""
    sc_consts = _str_constants_from_file(S4E9 / "synthesis_compliance.py")
    adb_consts = _str_constants_from_file(S4E9 / "audit_domestication_barcodes.py")
    foc_consts = _str_constants_from_file(RP / "fold_order_check.py")

    input_ref = sc_consts["INPUT_REF"]
    fasta_key = sc_consts["FASTA_KEY"]
    integration_ref = adb_consts["INTEGRATION_REF"]
    recoded_key = adb_consts["RECODED_KEY"]
    recoded_sha256 = adb_consts["RECODED_SHA256"]
    recodes_ref = adb_consts["RECODES_REF"]
    foc_fasta_key = foc_consts["FASTA_KEY"]
    foc_fasta_sha256 = foc_consts["FASTA_SHA256"]
    refold_commit = foc_consts["REFOLD_COMMIT"]

    assert integration_ref == input_ref, (
        f"audit_domestication_barcodes.INTEGRATION_REF ({integration_ref}) != "
        f"synthesis_compliance.INPUT_REF ({input_ref})"
    )
    assert recoded_key == fasta_key == foc_fasta_key, (
        f"order FASTA_KEY drift: synthesis_compliance={fasta_key!r}, "
        f"audit_domestication_barcodes={recoded_key!r}, fold_order_check={foc_fasta_key!r}"
    )
    assert len(foc_fasta_sha256) >= 7 and recoded_sha256.startswith(foc_fasta_sha256), (
        f"fold_order_check.FASTA_SHA256 ({foc_fasta_sha256!r}) is not a prefix of "
        f"audit_domestication_barcodes.RECODED_SHA256 ({recoded_sha256!r})"
    )
    assert len(refold_commit) >= 7 and input_ref.startswith(refold_commit), (
        f"fold_order_check.REFOLD_COMMIT ({refold_commit!r}) is not a prefix of "
        f"synthesis_compliance.INPUT_REF ({input_ref!r})"
    )

    snapshot_commits = {commit for commit, _ in err.PINNED_GIT_SNAPSHOTS}
    for label, commit in (
        ("INPUT_REF", input_ref),
        ("RECODES_REF", recodes_ref),
        ("REFOLD_COMMIT", refold_commit),
    ):
        assert commit in snapshot_commits, (
            f"{label} ({commit!r}) missing from export_review_repo.PINNED_GIT_SNAPSHOTS"
        )

    cut_state = json.loads((S4E9 / "final_cut_state.json").read_text(encoding="utf-8"))
    assert cut_state["main"] == input_ref, (
        f"final_cut_state.json main ({cut_state['main']!r}) != INPUT_REF ({input_ref!r})"
    )
    assert cut_state["fasta_key"] == fasta_key, (
        f"final_cut_state.json fasta_key ({cut_state['fasta_key']!r}) != FASTA_KEY ({fasta_key!r})"
    )
    assert cut_state["fasta_sha256"] == recoded_sha256, (
        f"final_cut_state.json fasta_sha256 ({cut_state['fasta_sha256']!r}) != "
        f"RECODED_SHA256 ({recoded_sha256!r})"
    )
    assert cut_state["published_from"] == recodes_ref, (
        f"final_cut_state.json published_from ({cut_state['published_from']!r}) != "
        f"RECODES_REF ({recodes_ref!r})"
    )

    comp_summary = json.loads(
        (S4E9 / "synthesis_compliance_summary.json").read_text(encoding="utf-8")
    )
    assert comp_summary["recoded_fasta"] == f"processed://{fasta_key}", (
        f"synthesis_compliance_summary.json recoded_fasta ({comp_summary['recoded_fasta']!r}) != "
        f"'processed://{fasta_key}'"
    )
    assert comp_summary["recoded_fasta_sha256"] == recoded_sha256, (
        f"synthesis_compliance_summary.json recoded_fasta_sha256 "
        f"({comp_summary['recoded_fasta_sha256']!r}) != RECODED_SHA256 ({recoded_sha256!r})"
    )
    assert comp_summary["input_refs"] and all(
        ref == input_ref for ref in comp_summary["input_refs"].values()
    ), f"synthesis_compliance_summary.json input_refs drifted from INPUT_REF ({input_ref!r})"

    numbering_rows = _read_tsv(S4E9 / "l1_insert_numbering.tsv")
    compliance_rows = _read_tsv(S4E9 / "synthesis_compliance.tsv")
    assert comp_summary["inserts"] == len(numbering_rows) == len(compliance_rows), (
        f"insert count mismatch: summary={comp_summary['inserts']}, "
        f"l1_insert_numbering.tsv={len(numbering_rows)}, "
        f"synthesis_compliance.tsv={len(compliance_rows)}"
    )

    empty_insert_ids = [r["arm_id"] for r in compliance_rows if not r["insert_id"].strip()]
    assert not empty_insert_ids, (
        f"synthesis_compliance.tsv has rows with empty insert_id: {empty_insert_ids[:10]}"
    )

    numbering_by_arm = {
        r["arm_id"]: (r["insert_id"].strip(), int(r["ordered_insert_nt"])) for r in numbering_rows
    }
    compliance_by_arm = {
        r["arm_id"]: (r["insert_id"].strip(), int(r["insert_nt"])) for r in compliance_rows
    }
    assert len(numbering_by_arm) == len(numbering_rows), (
        "duplicate arm_id in l1_insert_numbering.tsv"
    )
    assert len(compliance_by_arm) == len(compliance_rows), (
        "duplicate arm_id in synthesis_compliance.tsv"
    )
    assert compliance_by_arm == numbering_by_arm, (
        "synthesis_compliance.tsv (arm_id -> insert_id, insert_nt) drifted from "
        "l1_insert_numbering.tsv"
    )
