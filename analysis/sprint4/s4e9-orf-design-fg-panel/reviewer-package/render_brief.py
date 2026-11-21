"""Render the brief and the cover note to self-contained PDF and Word files, and copy them out.

For reading the brief outside the repository (e.g. a synced Google Drive folder, where a
.md file shows only as plain text). Each document is converted to HTML with its figures
embedded, its repository-relative links rewritten to GitHub links pinned to the commit, and
a footer naming that commit; headless Chrome prints it to PDF, and python-docx writes a .docx
with the figures embedded (Google Docs cannot fetch a Markdown file's relative image paths,
but it opens a .docx with its images in place). The Markdown files are copied
alongside. Render from a commit: a dirty tree is refused, since the PDF would cite a commit
that does not contain what it shows.

    uv run --locked python render_brief.py --dest "/path/to/FG Library Design"
"""

from __future__ import annotations

import argparse
import base64
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import markdown

HERE = Path(__file__).parent
REPO = HERE.parents[3]
#: Reviewers get the review repository, never this one, so every file link points there. Its
#: paths match this repository's, but its commits do not, so links are pinned to its branch
#: rather than to the commit this render came from; that commit is named in the footer instead.
GITHUB = "https://github.com/Deliverome-Project/l1-library-review/blob"
REVIEW_REF = "main"
#: This repository, for the footer's provenance line only — not for any link.
SOURCE_REPO = "deliverome-analysis"
DOCS = {
    "DESIGN.md": "L1_design_overview",
    "BRIEF.md": "L1_design_brief",
    "COVER_NOTE.md": "L1_design_brief_cover_note",
}
CHROME_CANDIDATES = (
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "google-chrome",
    "chromium",
)
#: Manrope is the house typeface. Chrome fetches it for the PDF; the .docx names it directly,
#: so a reader with Manrope installed sees it and anyone else falls back to their sans-serif.
FONT = "Manrope"
CSS = """
@import url('https://fonts.googleapis.com/css2?family=Manrope:wght@400;600;700&display=swap');
body { font-family: 'Manrope', -apple-system, Helvetica, Arial, sans-serif; font-size: 10.5pt;
       line-height: 1.4; max-width: 7.2in; margin: 0 auto; color: #1f1718; }
h1 { font-size: 17pt; } h2 { font-size: 13.5pt; margin-top: 1.4em; } h3 { font-size: 11.5pt; }
table { border-collapse: collapse; margin: 0.6em 0; font-size: 8.5pt; }
th, td { border: 1px solid #cfc4bf; padding: 3px 6px; vertical-align: top; }
th { background: #f3ece5; }
code { font-size: 9pt; }
img { max-width: 100%; }
.stamp { color: #6f5d5a; font-size: 8.5pt; border-top: 1px solid #cfc4bf; margin-top: 2em; }
"""


def chrome() -> str:
    for c in [os.environ.get("CHROME", ""), *CHROME_CANDIDATES]:
        if c and (Path(c).exists() or shutil.which(c)):
            return c
    raise SystemExit("no Chrome/Chromium found; set CHROME to its executable")


def git(*args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True
    ).stdout.strip()


def exported_paths() -> tuple[set[str], list[str]]:
    """(traced files, whole-folder prefixes) the review repository publishes."""
    sys.path.insert(0, str(HERE))
    import export_review_repo as ex

    traced = set((HERE / "review_repo_files.txt").read_text().split())
    return traced | ex.package_paths(), list(ex.ALWAYS)


def unbundled_citations(text: str) -> list[str]:
    """S3 locations ``text`` cites that the reviewer bundle does not carry.

    A reviewer works offline from the bundle, where ``bucket://key`` sits at
    ``s3/<bucket>/<key>``; a cited key outside it cannot be opened. A folder (a key ending
    in ``/``) passes when the bundle holds anything under it. ``bundle.py collect`` bundles
    every location cited here, so a failure means the manifest predates the citation.
    """
    import bundle

    with (HERE / "bundle_manifest.tsv").open() as fh:
        bundled = {f"{r['bucket']}/{r['key']}" for r in csv.DictReader(fh, delimiter="\t")}
    missing = []
    for bucket, key in bundle.cited_locations(text):
        loc = f"{bucket}/{key}"
        ok = any(b.startswith(loc) for b in bundled) if key.endswith("/") else loc in bundled
        if not ok:
            missing.append(f"{bucket}://{key}")
    return sorted(missing)


def body_html(md_path: Path, commit: str) -> str:
    """The document as HTML, with repository-relative links pointed at the review repository.

    A link to a file the review repository does not publish would 404 for the reader it is
    written for, so the render refuses rather than emitting one.
    """
    # `fenced_code` is not optional: without it a ``` block is parsed as inline code spans and
    # its lines are joined into one, which silently ran the rerun recipe's six shell commands
    # together on a single line in every Word and PDF copy.
    # `toc` is here only for the ids it puts on headings, which is what makes an in-document
    # link to the appendix resolve in the HTML and the PDF.
    text = md_path.read_text()
    if offline := unbundled_citations(text):
        listed = "\n  ".join(offline)
        raise SystemExit(
            f"{md_path.name} cites {len(offline)} S3 location(s) the reviewer bundle does not "
            f"carry, so a reviewer working offline cannot open them:\n  {listed}\n"
            "Re-collect the bundle (bundle.py collect, then manifest; collect bundles every "
            "location these documents cite), or stop citing them."
        )
    body = markdown.markdown(text, extensions=["tables", "fenced_code", "toc"])
    traced, always = exported_paths()
    unpublished: list[str] = []

    def link(m: re.Match) -> str:
        target = m.group(1)
        if re.match(r"^(https?:|mailto:|#)", target):
            return m.group(0)
        path, _, anchor = target.partition("#")
        rel = (md_path.parent / path).resolve().relative_to(REPO).as_posix()
        if not (rel in traced or any(rel == a or rel.startswith(a + "/") for a in always)):
            unpublished.append(rel)
        return f'href="{GITHUB}/{REVIEW_REF}/{rel}' + (f"#{anchor}" if anchor else "") + '"'

    out = re.sub(r'href="([^"]+)"', link, body)
    if unpublished:
        listed = "\n  ".join(sorted(set(unpublished)))
        raise SystemExit(
            f"{md_path.name} links {len(set(unpublished))} file(s) the review repository does "
            f"not publish, so they would 404 for a reviewer:\n  {listed}\n"
            "Either add them to export_review_repo.py's export set, or name them in the text "
            "without linking them."
        )
    return out


def stamp_text(md_path: Path, commit: str) -> tuple[str, str]:
    """(sentence, url): the source commit as text, linking the review repository's own copy."""
    url = f"{GITHUB}/{REVIEW_REF}/{md_path.relative_to(REPO).as_posix()}"
    return f"Rendered {datetime.now(UTC):%Y-%m-%d} from {SOURCE_REPO} {commit[:10]}.", url


def to_html(md_path: Path, commit: str) -> str:
    def embed(m: re.Match) -> str:
        png = base64.b64encode((md_path.parent / m.group(1)).read_bytes()).decode()
        return f'src="data:image/png;base64,{png}"'

    body = re.sub(r'src="(figures/[^"]+\.png)"', embed, body_html(md_path, commit))
    sentence, url = stamp_text(md_path, commit)
    stamp = f'<p class="stamp"><a href="{url}">{sentence}</a></p>'
    return f'<html><head><meta charset="utf-8"><style>{CSS}</style></head><body>{body}{stamp}</body></html>'


def _hyperlink(paragraph, url: str, text: str, style: dict) -> None:
    """python-docx has no hyperlink API; add the w:hyperlink run by hand."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    rid = paragraph.part.relate_to(
        url,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True,
    )
    link = OxmlElement("w:hyperlink")
    link.set(qn("r:id"), rid)
    run = OxmlElement("w:r")
    props = OxmlElement("w:rPr")
    for tag, val in (("w:color", "1F4E9A"), ("w:u", "single")):
        el = OxmlElement(tag)
        el.set(qn("w:val"), val)
        props.append(el)
    for flag in ("bold", "italic"):
        if style.get(flag):
            props.append(OxmlElement("w:b" if flag == "bold" else "w:i"))
    run.append(props)
    t = OxmlElement("w:t")
    t.text = text
    t.set(qn("xml:space"), "preserve")
    run.append(t)
    link.append(run)
    paragraph._p.append(link)


def _inline(paragraph, el, base: Path, style: dict | None = None) -> None:
    """Append an HTML element's inline content (text, b/i/code, links, images) to a paragraph."""
    from docx.shared import Inches, Pt

    style = dict(style or {})
    tag = el.tag if isinstance(el.tag, str) else ""
    if tag in ("strong", "b"):
        style["bold"] = True
    elif tag in ("em", "i"):
        style["italic"] = True
    elif tag == "code":
        style["code"] = True
    if tag == "img":
        paragraph.add_run().add_picture(str(base / el.get("src")), width=Inches(6.3))
        return
    if tag == "a" and el.get("href"):
        href = el.get("href")
        text = " ".join(el.text_content().split("\n"))
        if href.startswith("#"):
            # An in-document anchor resolves in the HTML and the PDF. Word would turn it into
            # an external link to a URL that does not exist, so it reads as plain text here.
            _run(paragraph, text, style, Pt)
        else:
            _hyperlink(paragraph, href, text, style)
        return
    if tag == "br":
        paragraph.add_run().add_break()
    if el.text and tag not in ("br",):
        _run(paragraph, el.text, style, Pt)
    for child in el:
        _inline(paragraph, child, base, style)
        if child.tail:
            _run(paragraph, child.tail, style, Pt)


def _run(paragraph, text: str, style: dict, pt=None) -> None:
    from docx.shared import Pt

    if not style.get("pre"):
        # Markdown source wraps paragraphs at ~100 columns; a newline there is a soft wrap,
        # not a line break, and Word and Google Docs would otherwise show it as one.
        text = " ".join(text.split("\n"))
    run = paragraph.add_run(text)
    run.bold, run.italic = style.get("bold"), style.get("italic")
    if style.get("code"):
        run.font.name, run.font.size = "Courier New", Pt(9)


def to_docx(md_path: Path, commit: str, out: Path) -> None:
    """A Word copy with the figures embedded, which Google Docs opens with them in place."""
    import lxml.html
    from docx import Document
    from docx.oxml.ns import qn
    from docx.shared import Pt

    doc = Document()
    doc.styles["Normal"].font.name = FONT
    doc.styles["Normal"].font.size = Pt(10.5)
    # Headings and the title use theme fonts (asciiTheme etc.), which override a plain font
    # name, so Word and Google Docs showed them in Calibri. Drop the theme attributes and name
    # Manrope for every script.
    for name in ("Title", "Heading 1", "Heading 2", "Heading 3", "Heading 4"):
        style = doc.styles[name]
        style.font.name = FONT
        rfonts = style.element.get_or_add_rPr().get_or_add_rFonts()
        for attr in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
            rfonts.attrib.pop(qn(attr), None)
        for attr in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
            rfonts.set(qn(attr), FONT)
    root = lxml.html.fragment_fromstring(body_html(md_path, commit), create_parent="div")
    base = md_path.parent
    for el in root:
        tag = el.tag
        if tag in ("h1", "h2", "h3", "h4"):
            _inline(doc.add_heading(level=int(tag[1]) - 1 or 0), el, base)
        elif tag == "p":
            _inline(doc.add_paragraph(), el, base)
        elif tag in ("ul", "ol"):
            for li in el.findall("li"):
                _inline(
                    doc.add_paragraph(style="List Bullet" if tag == "ul" else "List Number"),
                    li,
                    base,
                )
        elif tag == "table":
            rows = el.findall(".//tr")
            ncol = max(len(r.findall("th") + r.findall("td")) for r in rows)
            table = doc.add_table(rows=len(rows), cols=ncol)
            table.style = "Table Grid"
            for i, r in enumerate(rows):
                for j, c in enumerate(r.findall("th") + r.findall("td")):
                    cell = table.cell(i, j).paragraphs[0]
                    _inline(cell, c, base, {"bold": c.tag == "th"})
                    for run in cell.runs:  # tables are dense; body size made them hard to read
                        run.font.size = Pt(8.5)
            doc.add_paragraph()
        elif tag == "pre":
            _run(doc.add_paragraph(), el.text_content(), {"code": True, "pre": True})
        elif tag == "hr":
            doc.add_paragraph()
        else:
            _inline(doc.add_paragraph(), el, base)
    sentence, url = stamp_text(md_path, commit)
    p = doc.add_paragraph()
    _hyperlink(p, url, sentence, {"italic": True})
    doc.save(out)


def overwrite(src: Path, dst: Path) -> None:
    """Write ``src``'s bytes into ``dst`` in place, with a fresh modification time.

    A synced Drive folder keeps a file's ID, and so every link to it, only while the same
    file is rewritten; a copy that carries the source's older modification time (``copy2``)
    or replaces the file can reach Drive as a new file, leaving readers on a stale one.
    """
    data = src.read_bytes()
    with open(dst, "r+b" if dst.exists() else "wb") as fh:
        fh.write(data)
        fh.truncate()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dest", type=Path, required=True)
    a = ap.parse_args(argv)
    if git("status", "--porcelain", "--", str(HERE)):
        raise SystemExit("reviewer-package has uncommitted changes; commit before rendering")
    commit = git("rev-parse", "HEAD")
    dest = a.dest.expanduser()
    if not dest.is_dir():
        raise SystemExit(f"{dest} is not a folder")
    exe, work = chrome(), Path(tempfile.mkdtemp(prefix="render-brief-"))
    for name, stem in DOCS.items():
        html = work / f"{stem}.html"
        html.write_text(to_html(HERE / name, commit))
        pdf = work / f"{stem}.pdf"
        subprocess.run(
            [
                exe,
                "--headless=new",
                "--disable-gpu",
                "--no-pdf-header-footer",
                f"--print-to-pdf={pdf}",
                html.as_uri(),
            ],
            check=True,
            capture_output=True,
        )
        if not pdf.is_file() or pdf.stat().st_size == 0:
            raise SystemExit(f"Chrome produced no PDF for {name}")
        docx = work / f"{stem}.docx"
        to_docx(HERE / name, commit, docx)
        for f in (pdf, docx):
            overwrite(f, dest / f.name)
        overwrite(HERE / name, dest / f"{stem}.md")
        print(f"wrote {stem}.pdf, .docx and .md to {dest} ({commit[:10]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
