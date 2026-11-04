"""Normalise the provenance comments on the p0016-p0021 records.

Three defects, all introduced by the annotation fixes earlier in s4e9, none of which
touch a single base:

1. **Over-wide lines.** The notes name their script by repo-relative path. With the
   12-column GenBank indent, ``analysis/sprint4/s4e9-orf-design-fg-panel/\
fix_target_array_annotation.py`` is 72 characters, so the line runs to 84 and breaks
   the 80-column convention. Biopython warns and writes it anyway, so no text is lost,
   but a stricter parser on the vendor's side may reject the file. Naming the script by
   experiment folder rather than full repo path brings every word inside the limit.
2. **A duplicate note.** ``fix_target_array_annotation.py`` ran twice -- once to trim the
   flanks, once to annotate the four sites -- and appended the same paragraph each time.
3. **Two contradictory notes.** The handle retrofit also ran twice on 2026-10-06. The
   first pass wrote "The 2-nt GG pad is unchanged"; the second changed it to CC. Read in
   isolation the first is simply wrong about the record it sits on, so it is dropped and
   the surviving note states the end state.

Comments are rebuilt from paragraphs, not patched in place, so re-running is a no-op.

    uv run --locked python tidy_construct_notes.py --dry-run
    uv run --locked python tidy_construct_notes.py
"""

from __future__ import annotations

import argparse
import tempfile
import warnings
from pathlib import Path

from Bio import BiopythonWarning, SeqIO

from deliverome_analysis import seqtools
from deliverome_analysis.s3 import processed

DESIGNS = ("p0016", "p0017", "p0018", "p0019", "p0020", "p0021")
KEY = "constructs/designed/{}.gb"
LONG_PREFIX = "analysis/sprint4/s4e9-orf-design-fg-panel/"
SHORT_PREFIX = "s4e9-orf-design-fg-panel/"
# the superseded first-pass retrofit note, identified by the claim that is now false
STALE_MARKER = "The 2-nt GG pad is unchanged."


def paragraphs(comment: str) -> list[str]:
    """Split a wrapped comment into its notes. Every note ends with 'See <script>.py'."""
    flat = " ".join(comment.split())
    out, start = [], 0
    for i in range(len(flat)):
        if flat.startswith(".py", i) and (i + 3 == len(flat) or flat[i + 3] == " "):
            out.append(flat[start : i + 3].strip())
            start = i + 3
    tail = flat[start:].strip()
    if tail:
        out.append(tail)
    return [p for p in out if p]


def tidy(comment: str) -> tuple[str, list[str]]:
    """Return the rebuilt comment and a note of what was dropped or rewritten."""
    changes: list[str] = []
    kept: list[str] = []
    for para in paragraphs(comment):
        if STALE_MARKER in para:
            changes.append("dropped the superseded 'GG pad is unchanged' note")
            continue
        if LONG_PREFIX in para:
            para = para.replace(LONG_PREFIX, SHORT_PREFIX)
            changes.append("shortened a script path")
        if para in kept:
            changes.append("dropped a duplicate note")
            continue
        kept.append(para)
    return "\n".join(kept), changes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    staged = []
    for name in DESIGNS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            rec = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        before = str(rec.annotations.get("comment", ""))
        after, changes = tidy(before)
        if not changes:
            print(f"{name}: comment already tidy")
            continue
        if long := seqtools.overlong_genbank_words(after):
            raise SystemExit(f"{name}: words still too wide for GenBank: {long}")
        if tidy(after)[1]:
            raise SystemExit(f"{name}: tidy() is not idempotent")
        n_before, n_after = len(paragraphs(before)), len(paragraphs(after))
        print(f"{name}  notes {n_before} -> {n_after}   {'; '.join(sorted(set(changes)))}")
        rec.annotations["comment"] = after
        staged.append((name, rec, str(rec.seq).upper(), len(rec.features)))

    if args.dry_run:
        print("\ndry run — S3 not written")
        return 0

    for name, rec, seq, n_feat in staged:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"{name}.gb"
            SeqIO.write(rec, path, "genbank")
            if over := [ln for ln in path.read_text().splitlines() if len(ln) > 80]:
                raise SystemExit(f"{name}: {len(over)} line(s) still over 80 columns")
            processed.upload(path, KEY.format(name), content_type="chemical/seq-na-genbank")
        back = SeqIO.read(processed.download(KEY.format(name)), "genbank")
        if str(back.seq).upper() != seq:
            raise SystemExit(f"{name}: sequence changed, which this script must never do")
        if len(back.features) != n_feat:
            raise SystemExit(f"{name}: feature count changed {n_feat} -> {len(back.features)}")
        if tidy(str(back.annotations.get("comment", "")))[1]:
            raise SystemExit(f"{name}: re-read comment is still untidy")
    print(f"\nuploaded {len(staged)} record(s); sequence and features confirmed unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
