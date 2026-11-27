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
import functools
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
DEV_VIEWER = "https://surfaceome.deliverome.org"
#: Technical abbreviations that collide with human gene symbols (e.g. CPM = counts per million,
#: GPI = glycosylphosphatidylinositol anchor) and must never be auto-linked as proteins.
NON_PROTEIN_TOKENS = frozenset(
    {
        "ALFA",
        "APPRIS",
        "AUC",
        "CDS",
        "CPM",
        "CSV",
        "DNA",
        "EFS",
        "EUR",
        "FACS",
        "GPI",
        "GTF",
        "HEK",
        "HPA",
        "IRES",
        "ISS",
        "LTR",
        "MANE",
        "MSA",
        "NGS",
        "NMD",
        "ORF",
        "PAE",
        "PCR",
        "PDB",
        "PDF",
        "RBP",
        "RNA",
        "RSA",
        "SBS",
        "SHA",
        "TCEP",
        "TPM",
        "TSV",
        "URL",
        "WPRE",
    }
)
#: Reviewers get the review repository, never this one, so every file link points there. Its
#: paths match this repository's, but its commits do not, so links are pinned to its branch
#: rather than to the commit this render came from; that commit is named in the footer instead.
GITHUB = "https://github.com/Deliverome-Project/deliverome-surfaceome-library-v1/blob"
REVIEW_REF = "main"
#: This repository, for the footer's provenance line only — not for any link.
SOURCE_REPO = "deliverome-analysis"
DOCS = {
    "DESIGN.md": "L1_design_overview",
    "BRIEF.md": "L1_design_brief",
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
@page {
  margin: 0.72in 0.72in 0.82in 0.72in;
  @bottom-center {
    content: "Page " counter(page) " of " counter(pages);
    font-family: 'Manrope', -apple-system, Helvetica, Arial, sans-serif;
    font-size: 8.5pt;
    color: #6f5d5a;
  }
}
body { font-family: 'Manrope', -apple-system, Helvetica, Arial, sans-serif; font-size: 10.5pt;
       line-height: 1.4; max-width: 7.2in; margin: 0 auto; color: #1f1718; }
h1 { font-size: 17pt; } h2 { font-size: 13.5pt; margin-top: 1.4em; } h3 { font-size: 11.5pt; }
table { border-collapse: collapse; margin: 0.6em 0; font-size: 8.5pt; }
th, td { border: 1px solid #cfc4bf; padding: 3px 6px; vertical-align: top; }
th { background: #f3ece5; }
code { font-size: 9pt; }
img { max-width: 100%; }
p.caption { font-size: 8.5pt; color: #4a3f3d; margin-top: 0.3em; margin-bottom: 1.0em; }
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


@functools.lru_cache(maxsize=1)
def _protein_viewer_links() -> tuple[re.Pattern[str], dict[str, str], set[str]]:
    """Build (regex, symbol->url, active_design_ids) for dev.surfaceome.deliverome.org links."""
    s4e9 = HERE.parent
    urls: dict[str, str] = {}
    active_dids: set[str] = set()
    manifest = s4e9 / "candidate_manifest_tranches.tsv"
    if manifest.exists():
        with manifest.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("library", "").lower() == "l1" and not r.get("retired_reason"):
                    did = r["design_id"].strip()
                    active_dids.add(did)
                    if r.get("set") == "library_reference":
                        gene = r["key"].split(":", 1)[1] if ":" in r["key"] else r["key"]
                        if gene and gene not in NON_PROTEIN_TOKENS:
                            urls[gene] = f"{DEV_VIEWER}/v1-review/?m={did}"
    if "TFRC" in urls:
        urls["TfR1"] = urls["TFRC"]
    for path, col in (
        (HERE / "gate_lists/01_surface_evidence.tsv", "gene"),
        (s4e9 / "s3e8_gap_candidates.tsv", "gene_symbol"),
    ):
        if not path.exists():
            continue
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            key = col if reader.fieldnames and col in reader.fieldnames else None
            if not key and reader.fieldnames:
                key = reader.fieldnames[0]
            for r in reader:
                gene = (r.get(key or "") or "").strip()
                if gene and len(gene) >= 3 and gene not in urls and gene not in NON_PROTEIN_TOKENS:
                    urls[gene] = f"{DEV_VIEWER}/{gene}/"
    parts = [r"l1-p\d{4}(?:-a\d+)?"] + [re.escape(g) for g in sorted(urls, key=len, reverse=True)]
    pat = re.compile(r"(?<![A-Za-z0-9_/-])(" + "|".join(parts) + r")(?![A-Za-z0-9_-]|\.[A-Za-z])")
    return pat, urls, active_dids


def _split_protein_links(text: str, pat: re.Pattern[str], urls: dict[str, str], dids: set[str]):
    """Split text into alternating plain strings and (label, href) link tuples."""
    out: list[str | tuple[str, str]] = []
    pos = 0
    for m in pat.finditer(text):
        tok = m.group(1)
        href = urls.get(tok)
        if not href and tok.startswith("l1-p"):
            base_did = tok.split("-a", 1)[0]
            if base_did in dids:
                href = f"{DEV_VIEWER}/v1-review/?m={base_did}"
        if not href:
            continue
        if m.start() > pos:
            out.append(text[pos : m.start()])
        out.append((tok, href))
        pos = m.end()
    if pos == 0:
        return None
    if pos < len(text):
        out.append(text[pos:])
    return out


def _link_proteins_in_html(body: str) -> str:
    """Turn mentioned protein names and active L1 design IDs into dev-viewer links."""
    from lxml import html as lhtml

    pat, urls, dids = _protein_viewer_links()
    if not urls and not dids:
        return body
    root = lhtml.fromstring(f"<div>{body}</div>")
    skip_tags = {"a", "pre", "h1", "h2", "h3", "h4", "h5", "h6", "script", "style"}

    def _skip_code(el) -> bool:
        if el.tag != "code":
            return False
        txt = el.text_content()
        return any(ch in txt for ch in ("/", "_", "://", "=")) or bool(
            re.search(r"\.[a-z]{2,4}\b", txt)
        )

    def _walk(el) -> None:
        if not isinstance(el.tag, str) or el.tag in skip_tags or _skip_code(el):
            return
        children = list(el)
        if el.text:
            pieces = _split_protein_links(el.text, pat, urls, dids)
            if pieces:
                el.text = ""
                insert_at = 0
                for item in pieces:
                    if isinstance(item, str):
                        if insert_at == 0:
                            el.text = item
                        else:
                            el[insert_at - 1].tail = (el[insert_at - 1].tail or "") + item
                    else:
                        a = lhtml.Element("a", href=item[1])
                        a.text = item[0]
                        el.insert(insert_at, a)
                        insert_at += 1
        for child in children:
            _walk(child)
            if child.tail:
                pieces = _split_protein_links(child.tail, pat, urls, dids)
                if pieces:
                    child.tail = ""
                    idx = list(el).index(child)
                    cur = child
                    for item in pieces:
                        if isinstance(item, str):
                            cur.tail = (cur.tail or "") + item
                        else:
                            a = lhtml.Element("a", href=item[1])
                            a.text = item[0]
                            idx += 1
                            el.insert(idx, a)
                            cur = a

    _walk(root)
    return "".join(lhtml.tostring(c, encoding="unicode") for c in root)


def body_html(md_path: Path, commit: str, *, clean_paths: bool = False) -> str:
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
    import export_review_repo as ex

    unpublished: list[str] = []
    bad_r2_dirs: list[str] = []
    heading_ids = set(re.findall(r'<h[1-6][^>]*\sid="([^"]+)"', body))
    broken_anchors: list[str] = []

    def link(m: re.Match) -> str:
        target = m.group(1)
        if target.startswith("#"):
            if target[1:] not in heading_ids:
                broken_anchors.append(target)
            return m.group(0)
        if re.match(r"^(https?:|mailto:)", target):
            # Public Cloudflare R2 buckets do not serve directory listings; trailing-slash URLs 404.
            if ".r2.dev/" in target and target.endswith("/"):
                bad_r2_dirs.append(target)
            return m.group(0)
        path, _, anchor = target.partition("#")
        rel = (md_path.parent / path).resolve().relative_to(REPO).as_posix()
        if not (rel in traced or any(rel == a or rel.startswith(a + "/") for a in always)):
            unpublished.append(rel)
        dest_rel = ex._remap_path(rel) if clean_paths else rel
        return f'href="{GITHUB}/{REVIEW_REF}/{dest_rel}' + (f"#{anchor}" if anchor else "") + '"'

    out = re.sub(r'href="([^"]+)"', link, body)
    if broken_anchors:
        listed = "\n  ".join(sorted(set(broken_anchors)))
        raise SystemExit(
            f"{md_path.name} has {len(set(broken_anchors))} broken internal #anchor link(s):\n  {listed}"
        )
    if bad_r2_dirs:
        listed = "\n  ".join(sorted(set(bad_r2_dirs)))
        raise SystemExit(
            f"{md_path.name} links {len(set(bad_r2_dirs))} R2 directory URL(s) ending in '/' "
            f"(R2 public buckets return 404 on folder URLs; link a specific file or DOI instead):\n  {listed}"
        )
    if unpublished:
        listed = "\n  ".join(sorted(set(unpublished)))
        raise SystemExit(
            f"{md_path.name} links {len(set(unpublished))} file(s) the review repository does "
            f"not publish, so they would 404 for a reviewer:\n  {listed}\n"
            "Either add them to export_review_repo.py's export set, or name them in the text "
            "without linking them."
        )
    # Mark figure legends (<p><em>Figure N...</em></p>) with class="caption" for HTML/PDF styling.
    out = re.sub(r"<p>(<em>Figure\s+\d+)", r'<p class="caption">\1', out)
    # Link mentioned protein names and L1 design IDs directly to the dev surfaceome viewer.
    out = _link_proteins_in_html(out)
    return out


def stamp_text(md_path: Path, commit: str, *, clean_paths: bool = False) -> tuple[str, str]:
    """(sentence, url): the source commit as text, linking the review repository's own copy."""
    import export_review_repo as ex

    rel = md_path.relative_to(REPO).as_posix()
    dest_rel = ex._remap_path(rel) if clean_paths else rel
    url = f"{GITHUB}/{REVIEW_REF}/{dest_rel}"
    return f"Rendered {datetime.now(UTC):%Y-%m-%d} from {SOURCE_REPO} {commit[:10]}.", url


def to_html(md_path: Path, commit: str, *, clean_paths: bool = False) -> str:
    def embed(m: re.Match) -> str:
        png = base64.b64encode((md_path.parent / m.group(1)).read_bytes()).decode()
        return f'src="data:image/png;base64,{png}"'

    # Every repository-relative image is embedded, not only figures/: a figure drawn elsewhere in
    # the repository (../figures/..., the internalization schematic) otherwise prints as a broken
    # link in the PDF.
    body = re.sub(
        r'src="((?!https?:|data:)[^"]+\.png)"',
        embed,
        body_html(md_path, commit, clean_paths=clean_paths),
    )
    sentence, url = stamp_text(md_path, commit, clean_paths=clean_paths)
    stamp = f'<p class="stamp"><a href="{url}">{sentence}</a></p>'
    return f'<html><head><meta charset="utf-8"><style>{CSS}</style></head><body>{body}{stamp}</body></html>'


def _is_caption(el) -> bool:
    """True for a Markdown figure legend paragraph (``*Figure N...*``)."""
    if el.tag != "p":
        return False
    children = [c for c in el if isinstance(c.tag, str)]
    return (
        len(children) == 1
        and children[0].tag == "em"
        and not (el.text or "").strip()
        and not (children[0].tail or "").strip()
        and bool(re.match(r"^Figure\s+\d+", children[0].text_content().strip()))
    )


def _bookmark_name(html_id: str, used: set[str]) -> str:
    """Valid OOXML bookmark name (<=40 chars, ASCII letters/digits/underscore, starts with letter)."""
    slug = re.sub(r"[^A-Za-z0-9_]+", "_", html_id).strip("_")
    base = f"sec_{slug}"[:36]
    name = base
    n = 1
    while name in used:
        name = f"{base[:33]}_{n}"
        n += 1
    used.add(name)
    return name


def _add_bookmark(paragraph, bm_id: int, name: str) -> None:
    """Wrap ``paragraph``'s content in a Word bookmark so internal hyperlinks jump to it."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    start = OxmlElement("w:bookmarkStart")
    start.set(qn("w:id"), str(bm_id))
    start.set(qn("w:name"), name)
    end = OxmlElement("w:bookmarkEnd")
    end.set(qn("w:id"), str(bm_id))
    p = paragraph._p
    pPr = p.find(qn("w:pPr"))
    insert_at = list(p).index(pPr) + 1 if pPr is not None else 0
    p.insert(insert_at, start)
    p.append(end)


def _hyperlink_run(text: str, style: dict):
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    run = OxmlElement("w:r")
    props = OxmlElement("w:rPr")
    for tag, val in (("w:color", "1F4E9A"), ("w:u", "single")):
        el = OxmlElement(tag)
        el.set(qn("w:val"), val)
        props.append(el)
    for flag in ("bold", "italic"):
        if style.get(flag):
            props.append(OxmlElement("w:b" if flag == "bold" else "w:i"))
    if style.get("size_pt"):
        half_pt = str(round(style["size_pt"] * 2))
        for tag in ("w:sz", "w:szCs"):
            sz = OxmlElement(tag)
            sz.set(qn("w:val"), half_pt)
            props.append(sz)
    run.append(props)
    t = OxmlElement("w:t")
    t.text = text
    t.set(qn("xml:space"), "preserve")
    run.append(t)
    return run


def _hyperlink(paragraph, url: str, text: str, style: dict) -> None:
    """python-docx has no hyperlink API; add an external ``w:hyperlink`` run by hand."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    rid = paragraph.part.relate_to(
        url,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True,
    )
    link = OxmlElement("w:hyperlink")
    link.set(qn("r:id"), rid)
    link.append(_hyperlink_run(text, style))
    paragraph._p.append(link)


def _internal_hyperlink(paragraph, anchor: str, text: str, style: dict) -> None:
    """Internal in-document hyperlink targeting a heading bookmark (works in Word and Google Docs)."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    link = OxmlElement("w:hyperlink")
    link.set(qn("w:anchor"), anchor)
    link.append(_hyperlink_run(text, style))
    paragraph._p.append(link)


def _inline(
    paragraph, el, base: Path, style: dict | None = None, bookmarks: dict[str, str] | None = None
) -> None:
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
            target = (bookmarks or {}).get(href[1:])
            if target:
                _internal_hyperlink(paragraph, target, text, style)
            else:
                _run(paragraph, text, style, Pt)
        else:
            _hyperlink(paragraph, href, text, style)
        return
    if tag == "br":
        paragraph.add_run().add_break()
    if el.text and tag not in ("br",):
        _run(paragraph, el.text, style, Pt)
    for child in el:
        _inline(paragraph, child, base, style, bookmarks)
        if child.tail:
            _run(paragraph, child.tail, style, Pt)


def _run(paragraph, text: str, style: dict, pt=None) -> None:
    from docx.shared import Pt, RGBColor

    if not style.get("pre"):
        # Markdown source wraps paragraphs and bullet lists with indentation; collapsing any
        # whitespace run (newlines + continuation spaces) to a single space prevents extra
        # gaps inside wrapped list items in Word and Google Docs.
        text = re.sub(r"\s+", " ", text)
    run = paragraph.add_run(text)
    run.bold, run.italic = style.get("bold"), style.get("italic")
    if style.get("code"):
        run.font.name, run.font.size = "Courier New", Pt(style.get("size_pt") or 9)
    elif style.get("size_pt"):
        run.font.size = Pt(style["size_pt"])
    if style.get("color_rgb"):
        run.font.color.rgb = RGBColor(*style["color_rgb"])


def _add_page_number_footer(doc) -> None:
    """Add a centered 'Page X of Y' footer to every section of a .docx document."""
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Pt, RGBColor

    def _append_field(paragraph, field_name: str) -> None:
        run = paragraph.add_run()
        run.font.name = FONT
        run.font.size = Pt(8.5)
        run.font.color.rgb = RGBColor(0x6F, 0x5D, 0x5A)
        fld_begin = OxmlElement("w:fldChar")
        fld_begin.set(qn("w:fldCharType"), "begin")
        instr = OxmlElement("w:instrText")
        instr.set(qn("xml:space"), "preserve")
        instr.text = f" {field_name} "
        fld_sep = OxmlElement("w:fldChar")
        fld_sep.set(qn("w:fldCharType"), "separate")
        fld_end = OxmlElement("w:fldChar")
        fld_end.set(qn("w:fldCharType"), "end")
        run._r.append(fld_begin)
        run._r.append(instr)
        run._r.append(fld_sep)
        run._r.append(fld_end)

    for section in doc.sections:
        footer = section.footer
        footer.is_linked_to_previous = False
        p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        for txt, fld in (("Page ", "PAGE"), (" of ", "NUMPAGES")):
            r = p.add_run(txt)
            r.font.name = FONT
            r.font.size = Pt(8.5)
            r.font.color.rgb = RGBColor(0x6F, 0x5D, 0x5A)
            _append_field(p, fld)


def to_docx(md_path: Path, commit: str, out: Path, *, clean_paths: bool = False) -> None:
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
    root = lxml.html.fragment_fromstring(
        body_html(md_path, commit, clean_paths=clean_paths), create_parent="div"
    )
    base = md_path.parent
    used_bm: set[str] = set()
    bookmarks: dict[str, str] = {}
    for el in root:
        if el.tag in ("h1", "h2", "h3", "h4") and el.get("id"):
            bookmarks[el.get("id")] = _bookmark_name(el.get("id"), used_bm)
    bm_id = 0
    for el in root:
        tag = el.tag
        if tag in ("h1", "h2", "h3", "h4"):
            p = doc.add_heading(level=int(tag[1]) - 1 or 0)
            _inline(p, el, base, bookmarks=bookmarks)
            if el.get("id") and el.get("id") in bookmarks:
                _add_bookmark(p, bm_id, bookmarks[el.get("id")])
                bm_id += 1
        elif tag == "p":
            p_style = {"size_pt": 8.5, "color_rgb": (0x4A, 0x3F, 0x3D)} if _is_caption(el) else None
            _inline(doc.add_paragraph(), el, base, p_style, bookmarks)
        elif tag in ("ul", "ol"):
            for li in el.findall("li"):
                _inline(
                    doc.add_paragraph(style="List Bullet" if tag == "ul" else "List Number"),
                    li,
                    base,
                    bookmarks=bookmarks,
                )
        elif tag == "table":
            rows = el.findall(".//tr")
            ncol = max(len(r.findall("th") + r.findall("td")) for r in rows)
            table = doc.add_table(rows=len(rows), cols=ncol)
            table.style = "Table Grid"
            for i, r in enumerate(rows):
                for j, c in enumerate(r.findall("th") + r.findall("td")):
                    cell = table.cell(i, j).paragraphs[0]
                    _inline(cell, c, base, {"bold": c.tag == "th", "size_pt": 8.5}, bookmarks)
                    for run in cell.runs:  # tables are dense; body size made them hard to read
                        run.font.size = Pt(8.5)
            doc.add_paragraph()
        elif tag == "pre":
            _run(doc.add_paragraph(), el.text_content(), {"code": True, "pre": True})
        elif tag == "hr":
            doc.add_paragraph()
        else:
            _inline(doc.add_paragraph(), el, base, bookmarks=bookmarks)
    sentence, url = stamp_text(md_path, commit, clean_paths=clean_paths)
    p = doc.add_paragraph()
    _hyperlink(p, url, sentence, {"italic": True, "size_pt": 8.5})
    _add_page_number_footer(doc)
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
    ap.add_argument(
        "--clean-paths",
        action="store_true",
        help="rewrite repo links to library/l1/* clean paths (use once deliverome-surfaceome-library-v1 is re-exported with export_review_repo.py)",
    )
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
        html.write_text(to_html(HERE / name, commit, clean_paths=a.clean_paths))
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
        to_docx(HERE / name, commit, docx, clean_paths=a.clean_paths)
        for f in (pdf, docx):
            overwrite(f, dest / f.name)
        overwrite(HERE / name, dest / f"{stem}.md")
        print(f"wrote {stem}.pdf, .docx and .md to {dest} ({commit[:10]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
