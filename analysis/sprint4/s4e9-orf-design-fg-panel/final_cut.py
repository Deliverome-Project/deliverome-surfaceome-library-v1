"""Cut the L1 order of record from a commit on main, in three commands.

The order of record is every built arm on ``main``, recoded for synthesis, numbered and plated by
its ordered insert length. Cutting it touches a dozen scripts in a fixed order; this runs them, so
a re-cut after a change to one arm is the same three commands with a new commit and tag.

``prepare --main <sha> --tag <key-date>``
    Points every pin at ``<sha>`` (synthesis ``INPUT_REF``, the fold check's ``REFOLD_COMMIT``, the
    domestication audit's ``INTEGRATION_REF``, a new ``save_rung_builds`` entry) and the order's key
    at ``synthesis/<key-date>/``; reruns synthesis compliance on the rung builds at ``<sha>``;
    renumbers members by ordered insert length; re-keys the rung tables and their decision inputs
    to the new member ids; rebuilds the four rung arm tables and the C-terminal review, decisions
    and arms on those ids; reruns synthesis compliance so every insert carries its new member id;
    and proves that every rebuilt rung table equals ``<sha>``'s but for member ids.
    Then commit.

``publish`` (from a clean, committed tree)
    Publishes the recoded FASTA to its new key (stamped with the commit), pins its SHA-256 in the
    fold check and the audit, and saves the rung builds at ``<sha>`` to S3. Then commit.

``check``
    Runs every check on the cut: the fold check (each insert equals its fold), the consistency
    gate (every refold failure has a disposition), completeness, the order plan, the
    domestication and barcode audit, the brief numbers, the brief's prose-number test and the
    pinned-ref check. Then commit.

Nothing here uploads a rung table: S3 copies of rung and shared tables are published once, by
hand, after the coordinator's go-ahead.

Requested by the coordinator (2026-11-14: "build the cut so a small re-cut is one command");
implemented by the Claude Code session l1-reviewer-package.

    uv run --locked python final_cut.py prepare --main d6ab5093… --tag 2026-11-14-final-b
    # or, cut together with an arm change in the same PR:
    uv run --locked python final_cut.py prepare --main <branch sha> --tag … --in-pr
    git commit -am "…"
    uv run --locked python final_cut.py publish
    git commit -am "…"
    uv run --locked python final_cut.py check
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
REL = HERE.relative_to(REPO).as_posix()
STATE = HERE / "final_cut_state.json"
KEY_PREFIX = "sprint4/s4e9-orf-design-fg-panel/synthesis/"
#: Rung tables rebuilt on the new member ids, compared with main's but for member ids.
REKEYED = (
    "constructs_l1_nterm.tsv", "constructs_l1_nterm.fasta", "constructs_l1_nterm_summary.json",
    "reviewer-package/nterm/nterm_decisions.tsv",
    "constructs_l1_multipass.tsv", "constructs_l1_multipass.fasta", "multipass_arms_mature.fasta",
    "nterm_res1_multipass.tsv",
    "constructs_l1_type3.tsv", "constructs_l1_type3.fasta", "constructs_l1_type3_proteins.tsv",
    "constructs_l1_type3_proteins.fasta", "type3_decisions.tsv", "type3_nterm_calls.tsv",
    "type3_nterm_features.tsv",
    "constructs_l1_internal.tsv",
    "constructs_l1_cterm.tsv", "constructs_l1_cterm.fasta", "constructs_l1_cterm_summary.json",
    "l1_cterm_review.tsv", "l1_cterm_decisions.tsv",
)  # fmt: skip
MEMBER = re.compile(r"\bl1-m\d{4}\b")


def py(*args: str) -> None:
    print(f"$ python {' '.join(args)}", flush=True)
    subprocess.run([sys.executable, *args], cwd=HERE, check=True)


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO, check=True, capture_output=True,
                          text=True).stdout  # fmt: skip


def set_const(path: str, pattern: str, value: str) -> None:
    """Replace the first match of ``pattern`` (one group: the old value) in ``path``."""
    p = HERE / path
    text = p.read_text()
    new, n = re.subn(
        pattern, lambda m: m.group(0).replace(m.group(1), value), text, count=1, flags=re.M
    )
    if n != 1:
        raise SystemExit(f"{path}: pattern not found: {pattern}")
    p.write_text(new)


def fasta_key(tag: str) -> str:
    return f"{KEY_PREFIX}{tag}/l1_order_inserts_recoded.fasta"


def ids_only(sha: str) -> list[str]:
    """Rebuilt rung tables that differ from ``sha``'s in anything but ids (and line order).

    Every L1 id is masked on both sides (member ``l1-m``, design ``l1-p``, arm ``-aK``): from the
    per-insert cut a member is named by its design id, a dated member id survives inside paths
    to saved files, and JSON lists sorted by id reorder when ids change, which moves their
    trailing commas, so those are ignored too."""
    any_id = re.compile(r"\bl[12]-[mp]\d{4}(?:-a\d+)?\b")

    def norm(text: str) -> Counter:
        return Counter(any_id.sub("ID", ln).rstrip(",") for ln in text.splitlines())

    bad = []
    for name in REKEYED:
        if norm((HERE / name).read_text()) != norm(git("show", f"{sha}:{REL}/{name}")):
            bad.append(name)
    return bad


def prepare(sha: str, tag: str, in_pr: bool = False) -> int:
    sha = git("rev-parse", sha).strip()
    # A cut is normally pinned to main. With --in-pr it may be pinned to a commit on this
    # branch (a cut that lands together with the arm change it carries, as final-d did with
    # #352); the PR's merge commit then makes it a main ancestor, which pinned_refs checks.
    base = "HEAD" if in_pr else "origin/main"
    if subprocess.run(["git", "merge-base", "--is-ancestor", sha, base], cwd=REPO).returncode:
        raise SystemExit(f"{sha[:8]} is not on {base}")
    key = fasta_key(tag)
    set_const("synthesis_compliance.py", r'^INPUT_REF = "([0-9a-f]{40})"', sha)
    set_const("synthesis_compliance.py", r'FASTA_KEY = \(?\s*"(sprint4/[^"]+)"', key)
    set_const("reviewer-package/fold_order_check.py", r'^REFOLD_COMMIT = "([0-9a-f]+)"', sha[:8])
    set_const("reviewer-package/fold_order_check.py", r'FASTA_KEY = \(?\s*"(sprint4/[^"]+)"', key)
    set_const("audit_domestication_barcodes.py", r'^INTEGRATION_REF = "([0-9a-f]{40})"', sha)
    set_const("audit_domestication_barcodes.py", r'RECODED_KEY = \(?\s*"(sprint4/[^"]+)"', key)
    where = "cut in its PR" if in_pr else "main"
    rb = HERE / "save_rung_builds.py"
    text = rb.read_text()
    if sha not in text:
        rb.write_text(text.replace(
            "COMMITS = {\n",
            f'COMMITS = {{\n    "{sha}": "order of record {tag} ({where})",\n', 1))  # fmt: skip
    STATE.write_text(json.dumps({"main": sha, "tag": tag, "fasta_key": key}, indent=2) + "\n")

    py("synthesis_compliance.py")  # recodes and ordered lengths, on the builds at sha
    py("renumber_library_members.py", "--write")  # ids and plates by ordered length
    py("rekey_rung_tables.py", "--from", sha)
    for builder in ("build_l1_nterm_arms.py", "build_l1_multipass_arms.py",
                    "build_l1_type3_arms.py", "build_l1_internal_arms.py"):  # fmt: skip
        py(builder)
    py("review_l1_cterm_members.py")
    py("l1_cterm_decisions.py")
    py("build_l1_cterm_arms.py")
    py("synthesis_compliance.py")  # names every insert by its new member id
    py("synthesis_compliance.py", "--check")
    bad = ids_only(sha)
    if bad:
        raise SystemExit(f"rebuilt tables differ from {sha[:8]} beyond member ids: {bad}")
    print(f"prepared the cut from {sha[:8]} into {key}; every rebuilt rung table equals main's "
          "but for member ids. Commit, then run `final_cut.py publish`.")  # fmt: skip
    return 0


def publish() -> int:
    if git("status", "--porcelain", "--untracked-files=no").strip():
        raise SystemExit("commit first: the published FASTA is stamped with the commit")
    state = json.loads(STATE.read_text())
    py("save_rung_builds.py", "save")  # skips bundles already saved, so a rerun is a no-op
    if git("status", "--porcelain", "--untracked-files=no").strip():
        # The save recorded a new bundle in rung_builds_saved.json. The FASTA's stamp must name
        # a clean commit (final-d was first stamped dirty this way), so stop for a commit.
        raise SystemExit(
            "saved the rung builds: commit rung_builds_saved.json, then run publish again"
        )
    py("synthesis_compliance.py", "--publish")
    fasta_sha = json.loads((HERE / "synthesis_compliance_summary.json").read_text())[
        "recoded_fasta_sha256"
    ]
    head = git("rev-parse", "HEAD").strip()
    set_const(
        "reviewer-package/fold_order_check.py", r'FASTA_SHA256 = \(?\s*"([0-9a-f]+)"', fasta_sha[:8]
    )
    set_const("audit_domestication_barcodes.py", r'^RECODED_SHA256 = "([0-9a-f]{64})"', fasta_sha)
    set_const("audit_domestication_barcodes.py", r'^RECODES_REF = "([0-9a-f]{40})"', head)
    state.update({"fasta_sha256": fasta_sha, "published_from": head})
    STATE.write_text(json.dumps(state, indent=2) + "\n")
    print(
        f"published {state['fasta_key']} ({fasta_sha[:8]}). Commit, then run `final_cut.py check`."
    )
    return 0


def check() -> int:
    py("reviewer-package/fold_order_check.py")
    py("l1_consistency.py", "--check")
    py("l1_completeness.py", "--check")
    py("reviewer-package/order_plan.py")
    py("audit_domestication_barcodes.py")
    py("save_rung_builds.py", "check")
    py("reviewer-package/brief_numbers.py")
    py("pinned_refs.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "--no-header", "-p", "no:cacheprovider",
                    "tests/test_brief_prose_numbers.py", "tests/test_member_id_policy.py"],
                   cwd=HERE, check=True)  # fmt: skip
    print("the cut checks out. Commit, update the brief's prose from brief_numbers.json "
          "(the prose test names any number that moved), then re-render the brief.")  # fmt: skip
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--main", required=True, help="a commit on origin/main holding every rung build")
    p.add_argument("--tag", required=True, help="the dated key folder, e.g. 2026-11-14-final-b")
    p.add_argument(
        "--in-pr",
        action="store_true",
        help="allow a commit on this branch that its PR will merge into main",
    )
    sub.add_parser("publish")
    sub.add_parser("check")
    args = ap.parse_args(argv)
    if args.cmd == "prepare":
        return prepare(args.main, args.tag, args.in_pr)
    return publish() if args.cmd == "publish" else check()


if __name__ == "__main__":
    raise SystemExit(main())
