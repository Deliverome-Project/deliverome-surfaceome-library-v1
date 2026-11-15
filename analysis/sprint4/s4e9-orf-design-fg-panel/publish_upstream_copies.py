"""Copy the upstream reference inputs that were saved in processed:// to external://.

CLAUDE.md puts reference inputs in the append-only external bucket under ``{source}/{version}/``;
the overwrite-capable processed bucket can be changed under a pin. The pipeline's saved inputs
included 13 copies of third-party data held in processed:// (GENCODE Ribo-seq ORFs, PeptideAtlas
query responses, MyGene.info responses, UniProt REST results). This copies each, byte for byte,
to a new external key, never overwriting one, checks the copy's SHA-256 after upload, and
records the mapping in ``upstream_copies.tsv``. The processed objects stay where they are.

Decided by the coordinator (2026-11-14, after the independent audit), under Becca's standing
instruction on uploads (new append-only keys only); implemented by the Claude Code session
l1-reviewer-package.

    uv run --locked python publish_upstream_copies.py publish
    uv run --locked python publish_upstream_copies.py check
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "upstream_copies.tsv"
P = "sprint4/s4e9-orf-design-fg-panel/"
PI = P + "library-design/pipeline-inputs/"
RIBO = [
    f"Ribo-seq_ORFs.{s}.{e}" for s in ("all", "comprehensive", "primary") for e in ("bed", "faa")
]

#: processed key -> (external key, source description)
COPIES: dict[str, tuple[str, str]] = {
    **{
        f"{PI}2026-11-01/{n}": (
            f"gencode-riboseq-orfs/fetched-2026-11-01/{n}",
            "GENCODE Ribo-seq ORFs (phase 1 all; phase 2 primary, comprehensive), "
            "https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/, fetched 2026-11-01",
        )
        for n in RIBO
    },
    f"{PI}2026-11-01/peptideatlas_build607_pools.json": (
        "peptideatlas/build-607/queries-2026-11-01/peptideatlas_build607_pools.json",
        "PeptideAtlas GetPeptides responses, atlas_build_id 607, queried 2026-11-01",
    ),
    f"{PI}2026-11-02/s1e3/mygene_response.json": (
        "mygene/2026-11-02/s1e3/mygene_response.json",
        "MyGene.info symbol resolution (plain expression join), queried 2026-11-02",
    ),
    f"{PI}2026-11-02/s1e3/mygene_response_full.json": (
        "mygene/2026-11-02/s1e3/mygene_response_full.json",
        "MyGene.info symbol resolution (--full expression join), queried 2026-11-02",
    ),
    f"{P}library-design/cterm-review/2026-11-12-second-pass/uniprot_entries.json.gz": (
        "uniprot/rest-2026-11-12/cterm-review-second-pass/uniprot_entries.json.gz",
        "UniProt REST reviewed entries for the C-terminal members, fetched 2026-11-12",
    ),
    f"{P}library-design/cterm-review/2026-11-13-slc68a1/uniprot_entries.json.gz": (
        "uniprot/rest-2026-11-13/cterm-review-slc68a1/uniprot_entries.json.gz",
        "the 2026-11-12 entries plus SLC68A1, fetched 2026-11-13",
    ),
    f"{P}library-design/panel-inputs/2026-11-01/uniprot_human_gpi.tsv": (
        "uniprot/2026_03/rest-2026-11-01/uniprot_human_gpi.tsv",
        "UniProt REST release 2026_03: keyword:KW-0336 AND organism_id:9606 AND reviewed:true",
    ),
    f"{P}reviewer-package/nterm/uniprot-propeptides/2026_04/uniprot_features.tsv": (
        "uniprot/2026_04/nterm-propeptides/uniprot_features.tsv",
        "UniProt 2026_04 PROPEP, SIGNAL and CHAIN features, fetched 2026-11-13",
    ),
}
COLS = ("processed_key", "external_key", "sha256", "bytes", "source")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def publish() -> int:
    from deliverome_analysis.s3 import external, processed

    rows = []
    for src, (dst, what) in COPIES.items():
        local = Path(processed.download(src))
        want = sha256(local)
        if external.exists(dst):
            if sha256(Path(external.download(dst))) != want:
                raise SystemExit(f"external://{dst} exists with other content; not overwritten")
        else:
            external.upload(local, dst)
            if sha256(Path(external.download(dst))) != want:
                raise SystemExit(f"external://{dst}: uploaded copy differs from the source")
        rows.append({"processed_key": f"processed://{src}", "external_key": f"external://{dst}",
                     "sha256": want, "bytes": local.stat().st_size, "source": what})  # fmt: skip
        print(f"external://{dst}")
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return 0


def check() -> str:
    from deliverome_analysis.s3 import external

    with OUT.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    bad = [r["external_key"] for r in rows
           if sha256(Path(external.download(r["external_key"].removeprefix("external://"))))
           != r["sha256"]]  # fmt: skip
    if bad:
        raise SystemExit(f"external copies changed: {bad}")
    return f"{len(rows)} upstream reference copies in external:// unchanged"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("publish", "check"))
    args = ap.parse_args(argv)
    if args.command == "publish":
        return publish()
    print(check())
    return 0


if __name__ == "__main__":
    sys.exit(main())
