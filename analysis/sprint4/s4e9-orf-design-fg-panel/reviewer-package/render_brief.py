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
import os
import re
import shutil
import subprocess
import tempfile
from datetime import UTC, datetime
from pathlib import Path

import markdown

HERE = Path(__file__).parent
REPO = HERE.parents[3]
GITHUB = "https://github.com/Deliverome-Project/deliverome-analysis/blob"
DOCS = {"BRIEF.md": "L1_design_brief", "COVER_NOTE.md": "L1_design_brief_cover_note"}
CHROME_CANDIDATES = (
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "google-chrome",
    "chromium",
)
CSS = """
body { font-family: -apple-system, Helvetica, Arial, sans-serif; font-size: 10.5pt;
       line-height: 1.4; max-width: 7.2in; margin: 0 auto; color: #1f1718; }
h1 { font-size: 17pt; } h2 { font-size: 13.5pt; margin-top: 1.4em; } h3 { font-size: 11.5pt; }
table { border-collapse: collapse; margin: 0.6em 0; font-size: 9.5pt; }
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


def body_html(md_path: Path, commit: str) -> str:
    """The document as HTML, with repository-relative links pinned to ``commit``."""
    body = markdown.markdown(md_path.read_text(), extensions=["tables"])

    def link(m: re.Match) -> str:
        target = m.group(1)
        if re.match(r"^(https?:|mailto:|#)", target):
            return m.group(0)
        path, _, anchor = target.partition("#")
        rel = (md_path.parent / path).resolve().relative_to(REPO).as_posix()
        return f'href="{GITHUB}/{commit}/{rel}' + (f"#{anchor}" if anchor else "") + '"'

    return re.sub(r'href="([^"]+)"', link, body)


def stamp_text(md_path: Path, commit: str) -> tuple[str, str]:
    """(sentence, url) naming the commit a rendering was made from."""
    url = f"{GITHUB}/{commit}/{md_path.relative_to(REPO).as_posix()}"
    return f"Rendered {datetime.now(UTC):%Y-%m-%d} from deliverome-analysis {commit[:10]}.", url


def to_html(md_path: Path, commit: str) -> str:
    def embed(m: re.Match) -> str:
        png = base64.b64encode((md_path.parent / m.group(1)).read_bytes()).decode()
        return f'src="data:image/png;base64,{png}"'

    body = re.sub(r'src="(figures/[^"]+\.png)"', embed, body_html(md_path, commit))
    stamp = (
        f'<p class="stamp">Rendered {datetime.now(UTC):%Y-%m-%d} from '
        f'<a href="{GITHUB}/{commit}/{md_path.relative_to(REPO).as_posix()}">'
        f"deliverome-analysis {commit[:10]}</a>.</p>"
    )
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
        _hyperlink(paragraph, el.get("href"), " ".join(el.text_content().split("\n")), style)
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
    from docx.shared import Pt

    doc = Document()
    doc.styles["Normal"].font.name = "Arial"
    doc.styles["Normal"].font.size = Pt(10.5)
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
            shutil.copy2(f, dest / f.name)
        shutil.copy2(HERE / name, dest / f"{stem}.md")
        print(f"wrote {stem}.pdf, .docx and .md to {dest} ({commit[:10]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
