"""CODE_MAP.md lists every s4e9 script (regenerate with code_map.py)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import code_map


def test_code_map_lists_every_script():
    missing, extra = code_map.check_listed()
    assert not missing and not extra, (
        f"run code_map.py: not listed {missing}; listed but gone {extra}"
    )
