"""The internal-rung arms rebuild byte for byte from the ranker's decisions.

The rebuild reads the domesticated ORFs from S3 (``library_design.CONSTRUCTS_KEY``), so it skips
without credentials, like the repo's other live tests. The committed arms are checked offline in
CI by ``test_l1_checks.py`` (translation, layout, fold, verdict, disposition).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from deliverome_analysis.core import s3  # noqa: E402

pytestmark = pytest.mark.skipif(
    not s3.credentials_available(), reason="rebuild reads the domesticated ORFs from S3"
)


def test_internal_arms_rebuild_byte_for_byte():
    import build_l1_internal_arms

    assert build_l1_internal_arms.main(["--check"]) == 0
