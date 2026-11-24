"""Run one ``library_pipeline.py verify`` stage, in this repo or in an exported review repo.

Two stages (``isoform_universe``, ``ortholog_repick``) rerun scripts as they were at the
commits that produced their records, which ``reproduce_isoform_universe.code_at`` reads with
``git archive``. An exported review repo has no such history, so ``export_review_repo.py``
ships each of those archives as ``historical/<commit>.tar``, made by the same command. When
that folder exists, this reads them from there; nothing else about the stage changes.

    python run_stage.py reference
"""

from __future__ import annotations

import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).parent
S4E9 = HERE.parent
REPO = S4E9.parents[2]
HISTORICAL = REPO / "historical"


def historical_code_at(commit: str, dest: Path) -> Path:
    import reproduce_isoform_universe as riu

    if not dest.exists():
        dest.mkdir(parents=True)
        with tarfile.open(HISTORICAL / f"{commit}.tar") as tar:
            tar.extractall(dest, filter="data")
    return dest / riu.AP


def main(stage: str) -> int:
    sys.path.insert(0, str(S4E9))
    if HISTORICAL.is_dir():
        import reproduce_isoform_universe as riu

        riu.code_at = historical_code_at
    import library_pipeline

    return library_pipeline.main(["verify", "--stage", stage])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
