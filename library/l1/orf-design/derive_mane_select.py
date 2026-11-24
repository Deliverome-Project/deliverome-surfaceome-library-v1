"""MANE Select per gene, read from the pinned Ensembl 115 GTF's ``tag "MANE_Select"``.

The reference ladder (``select_reference.py``) took MANE from the APPRIS e115v50 file's MANE
column, which omits 775 of the GTF's 19,276 MANE Select transcripts; 66 library genes then
looked "APPRIS only" though Ensembl tags a MANE Select for them (proteoform review,
2026-11-22). This writes the GTF's table, which the ladder reads as a pinned input
(``mane_select_ensembl115.tsv`` in ``library_pipeline_inputs.json``). Where the APPRIS column
names a MANE transcript the two agree for every gene; ``--check`` asserts that.

Decided by Becca Carlson, 2026-11-22 ("Yes, in this build", to rerunning the ladder with GTF
MANE); implemented by Claude Code session l1-review-export-scope (subagent).

    uv run --locked python derive_mane_select.py            # write and compare
    uv run --locked python derive_mane_select.py --publish  # upload once and pin
"""

from __future__ import annotations

import argparse
import gzip
import json
import tempfile
from pathlib import Path

import library_pipeline as lp

NAME = "mane_select_ensembl115.tsv"
KEY = "ensembl/release-115/derived/mane-select-from-gtf/" + NAME


def derive(gtf: Path) -> str:
    rows: dict[str, tuple[str, str, str]] = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "transcript" or 'tag "MANE_Select"' not in f[8]:
                continue

            def attr(k: str, a: str = f[8]) -> str:
                return a.split(f'{k} "')[1].split('"')[0] if f'{k} "' in a else ""

            gid = attr("gene_id")
            assert gid not in rows, f"{gid}: two MANE Select transcripts"
            rows[gid] = (attr("gene_name"), attr("transcript_id"), attr("transcript_biotype"))
    lines = ["gene_id\tgene_name\ttranscript_id\ttranscript_biotype"]
    lines += ["\t".join((g, *rows[g])) for g in sorted(rows)]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--publish", action="store_true")
    a = ap.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="mane-select-"))
    got = lp.stage_inputs(
        ["Homo_sapiens.GRCh38.115.gtf.gz", "appris_homo_sapiens_e115v50.txt"], work
    )
    text = derive(got["Homo_sapiens.GRCh38.115.gtf.gz"])
    table = {ln.split("\t")[0]: ln.split("\t")[2] for ln in text.splitlines()[1:]}
    appris = {}
    for line in got["appris_homo_sapiens_e115v50.txt"].open():
        f = line.rstrip("\n").split("\t")
        if len(f) > 5 and "MANE_Select" in f[5]:
            appris[f[1].split(".")[0]] = f[2].split(".")[0]
    bad = {g: (t, table.get(g)) for g, t in appris.items() if table.get(g) != t}
    assert not bad, f"APPRIS and GTF MANE disagree: {list(bad.items())[:5]}"
    print(f"{len(table)} MANE Select genes; APPRIS names {len(appris)}, all agreeing")
    out = work / NAME
    out.write_text(text)
    print(f"sha256 {lp.sha256(out)}")
    if a.publish:
        from deliverome_analysis.s3 import external

        spec = json.loads(lp.INPUTS.read_text())
        if NAME in spec["files"]:
            assert spec["files"][NAME]["sha256"] == lp.sha256(out), "pinned table differs"
            print("already pinned")
            return 0
        if external.exists(KEY):
            raise SystemExit(f"external://{KEY} exists")
        external.upload(out, KEY)
        spec["files"][NAME] = {
            "source": f"external:{KEY}",
            "sha256": lp.sha256(out),
            "origin": "derive_mane_select.py from the pinned Homo_sapiens.GRCh38.115.gtf.gz "
            '(transcripts tagged "MANE_Select"), 2026-11-22',
        }
        lp.INPUTS.write_text(json.dumps(spec, indent=1, ensure_ascii=False) + "\n")
        print(f"published external://{KEY} and pinned it in {lp.INPUTS.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
