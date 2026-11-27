"""The pinned inputs of the isoform and ortholog panel draws, staged from S3.

``sample_isoform_panel.py`` and ``sample_ortholog_panel.py`` were first run against files on
one laptop: a UniProt cache, an accessible-surfaceome checkout and intermediate tables
nobody committed. Every input they read is now one dated snapshot on S3, checksummed in
``panel_inputs.json``, and ``--snapshot`` on either sampler stages it and reads nothing
else. With an empty rejection list the two samplers reproduce the published 150-isoform
and 200-ortholog panels byte for byte (checked 2026-11-01).

    uv run --locked python panel_inputs.py --publish DIR   # once: upload a snapshot from DIR
    uv run --locked python panel_inputs.py --stage DIR     # fetch and verify into DIR
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).parent
MANIFEST = HERE / "panel_inputs.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage(dest: str | Path) -> Path:
    """Download every pinned input into ``dest`` and verify its checksum."""
    from deliverome_analysis.core.s3 import processed

    spec = json.loads(MANIFEST.read_text())
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    for name, entry in spec["files"].items():
        path = dest / name
        if not path.exists() or _sha256(path) != entry["sha256"]:
            processed.download(spec["prefix"] + name, path)
        if _sha256(path) != entry["sha256"]:
            raise SystemExit(f"{name}: checksum does not match panel_inputs.json")
    return dest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", type=Path)
    ap.add_argument("--publish", type=Path, help="upload DIR's files under the manifest prefix")
    args = ap.parse_args(argv)
    spec = json.loads(MANIFEST.read_text())
    if args.publish:
        from deliverome_analysis.core.s3 import processed

        for name, entry in spec["files"].items():
            path = args.publish / name
            if _sha256(path) != entry["sha256"]:
                raise SystemExit(f"{name}: local file does not match panel_inputs.json")
            processed.upload(str(path), spec["prefix"] + name)
            print(f"published processed://{spec['prefix']}{name}")
    if args.stage:
        print(f"staged {len(spec['files'])} verified inputs in {stage(args.stage)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
