"""No s4e9 script reads a live external source without a saved copy on record.

``external_inputs.scan`` is offline: every script that makes a network call must be listed in
``external_inputs.tsv``, a ``verify`` stage must read a saved copy, and every committed copy
must match its SHA-256. The S3 copies are checked by the pipeline stage ``external_inputs``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import external_inputs


def test_every_network_call_has_a_saved_copy_on_record():
    problems = external_inputs.scan()
    assert not problems, "\n".join(problems)


def test_the_scan_sees_a_live_fetch(tmp_path, monkeypatch):
    """The detector itself works: an unlisted script with urlopen( is reported."""
    script = external_inputs.HERE / "zz_probe_live_fetch.py"
    monkeypatch.setattr(external_inputs, "fetching_scripts", lambda: ["zz_probe_live_fetch.py"])
    assert any("zz_probe_live_fetch.py" in p for p in external_inputs.scan())
    assert external_inputs.NETWORK.search("urllib.request.urlopen(url)")
    assert not external_inputs.NETWORK.search("see https://rest.uniprot.org in the docstring")
    assert not script.exists()
