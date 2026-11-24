"""EndoNBV2 internalization assay, re-set to sit in a document.

The drawing is the Deliverome platform figure itself -- same art, same
glyphs, same colours, same story left to right.  Nothing here redraws it.
What this script changes is the one thing that stops it going on a page:

    at 6.5 in wide, the master's own lettering lands between 4.4 and 8.4 pt.

So the master's `text` group is dropped and every label is re-set in Manrope
at a single scale, the two reagent definitions move out of the artwork into a
key strip across the top where they have room, the sheet is cropped to the
ink, and a source line is added at the foot.  The art is placed once, as a
whole, and only the type moves.

    python3 library/l1/figures/fig_internalization_assay.py

See `assets/MANIFEST.md` for where the art comes from and what is vendored.
It is a schematic, not a data figure -- nothing on it is measured -- so, like
the other schematics in scratch/, it writes SVG directly rather than going
through `save_figure()` and carries no provenance stamp.
"""

from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ASSETS = HERE / "assets"
OUT_DIR = HERE / "figures"
STEM = "fro_internalization_assay"

MASTER = "platform-internalization_20260522"

# ---------------------------------------------------------------------- ink
# Read off the master itself, so the lettering matches the art it names.
AMBER = "#F4AA28"      # SpyCatcher / SpyTag
AMBER_DK = "#B07415"   # ...dark enough to read as type
GREEN = "#61BB46"      # ALFA nanobody / ALFA tag
GREEN_DK = "#3E8A2C"
INK = "#2B2B2B"
NOTE = "#6E6E6E"
RULE = "#E2E2E2"
FONT = "Manrope, 'Helvetica Neue', Arial, sans-serif"

# ------------------------------------------------------------------- layout
# The master is 2011.58 x 970.12 and its ink runs x 89..1997, y 47..930.
# Crop to that, drop a key strip above it and a source line below.
CROP_X, CROP_Y = 56.0, 20.0
ART_W, ART_H = 1980.0, 924.0
KEY_H = 104.0            # the key strip, above the art
FOOT_H = 10.0            # whitespace under the art; no source line on the
                         # sheet -- the citation belongs in the caption

W = ART_W
H = KEY_H + ART_H + FOOT_H

# 1980 units = 6.5 in, so 1 unit = 0.2364 pt.
BIG = 38      # 9.0 pt -- the key strip
MED = 36      # 8.5 pt -- labels that carry a step
SML = 34      # 8.0 pt -- the rest
TAG = 30      # 7.1 pt -- the two inline tag names, and only those

MARGIN = 44


def art(x: float, y: float) -> tuple[float, float]:
    """A point in the master's coordinates, in sheet coordinates."""
    return x - CROP_X, y - CROP_Y + KEY_H


# -------------------------------------------------------------------- text
def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, *, size=SML, fill=INK, weight=500, anchor="start"):
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}">{esc(s)}</text>'
    )


_NARROW = set("iljt.,;:!'|()[]-/ ")
_WIDE = set("mwMW%")


def text_width(s: str, size: float, bold: bool = False) -> float:
    w = 0.0
    for ch in s:
        if ch in _NARROW:
            w += 0.30
        elif ch in _WIDE:
            w += 0.83
        elif ch.isupper() or ch.isdigit():
            w += 0.68
        else:
            w += 0.52
    return w * size * (1.07 if bold else 1.0)


def lines(x, y, rows, *, size=SML, fill=INK, weight=500, lead=None,
          anchor="start"):
    lead = lead or size * 1.22
    return "".join(
        text(x, y + i * lead, row, size=size, fill=fill, weight=weight,
             anchor=anchor)
        for i, row in enumerate(rows)
    )


# --------------------------------------------------------------- the master
class Asset:
    """One vendored SVG, namespaced so several can share a document."""

    _n = 0

    def __init__(self, name: str, *, drop: tuple[str, ...] = ()):
        Asset._n += 1
        raw = (ASSETS / f"{name}.svg").read_text(encoding="utf-8")
        vb = re.search(r'viewBox="([^"]+)"', raw).group(1).split()
        self.name = name
        self.x0, self.y0, self.w, self.h = (float(v) for v in vb)
        body = raw[raw.index(">", raw.index("<svg")) + 1:raw.rindex("</svg>")]
        for gid in drop:
            body = _drop_group(body, gid)
        tag = f"a{Asset._n}" + re.sub(r"[^a-z0-9]+", "", name.lower())[:10]
        body = re.sub(r"cls-(\d+)", rf"{tag}-cls-\1", body)
        body = re.sub(r'id="([^"]+)"', rf'id="{tag}-\1"', body)
        body = re.sub(r"url\(#([^)]+)\)", rf"url(#{tag}-\1)", body)
        self.body = body

    def place(self, x, y, *, w=None, h=None, deg=0.0, anchor="c"):
        s = (w / self.w) if w else ((h / self.h) if h else 1.0)
        fx = {"c": 0.5, "n": 0.5, "s": 0.5, "e": 1.0, "w": 0.0,
              "nw": 0.0, "ne": 1.0, "sw": 0.0, "se": 1.0}[anchor]
        fy = {"c": 0.5, "n": 0.0, "s": 1.0, "e": 0.5, "w": 0.5,
              "nw": 0.0, "ne": 0.0, "sw": 1.0, "se": 1.0}[anchor]
        ax, ay = self.x0 + fx * self.w, self.y0 + fy * self.h
        return (
            f'<g transform="translate({x:.2f},{y:.2f}) rotate({deg:.2f}) '
            f'scale({s:.5f}) translate({-ax:.2f},{-ay:.2f})">{self.body}</g>'
        )


def _drop_group(body: str, gid: str) -> str:
    """Remove <g id="gid"> ... </g>, counting nested groups."""
    m = re.search(rf'<g id="{gid}"[^>]*>', body)
    if not m:
        raise SystemExit(f"{gid!r} not found in the master; did it change?")
    i, depth = m.end(), 1
    while depth:
        nxt = re.search(r"<g\b|</g>", body[i:])
        if not nxt:
            raise SystemExit(f"unbalanced groups while dropping {gid!r}")
        i += nxt.end()
        depth += 1 if nxt.group(0) == "<g" else -1
    return body[: m.start()] + body[i:]


ART = Asset(MASTER, drop=("text",))
HEX_SPY = Asset("fluorophore-hexagons-orange_050")
HEX_ALFA = Asset("fluorophore-hexagons-green_051")


# ----------------------------------------------------------------- the key
KEY = (
    (HEX_SPY, "SpyCatcher", AMBER_DK, "normalizes for expression"),
    (HEX_ALFA, "anti-ALFA nanobody", GREEN_DK,
     "cleavable on the cell surface"),
)


def key_strip() -> str:
    """The two reagent definitions, lifted out of the art and given room.

    In the master these sit inside the drawing as three- and four-line blocks.
    At page size there is no width there for type this size, and they are the
    two things a reader has to have before the drawing means anything, so they
    come out to the top where the full sheet width is available.
    """
    o = []
    widths = []
    for _, name, _, rest in KEY:
        widths.append(58 + text_width(f"{name} — {rest}", BIG, True))
    gap = (W - 2 * MARGIN - sum(widths)) / (len(KEY) - 1)
    if gap < 48:
        raise SystemExit(f"key strip is {-gap + 48:.0f} units too wide")
    x = MARGIN
    for (swatch, name, colour, rest), w in zip(KEY, widths):
        o.append(swatch.place(x + 20, 54, h=58, anchor="c"))
        o.append(
            f'<text x="{x + 58:.1f}" y="66" font-family="{FONT}" '
            f'font-size="{BIG}" fill="{INK}">'
            f'<tspan font-weight="700" fill="{colour}">{esc(name)}</tspan>'
            f'<tspan font-weight="400"> — {esc(rest)}</tspan></text>'
        )
        x += w + gap
    o.append(f'<line x1="{MARGIN}" y1="{KEY_H - 14:.0f}" '
             f'x2="{W - MARGIN}" y2="{KEY_H - 14:.0f}" stroke="{RULE}" '
             f'stroke-width="2"/>')
    return "".join(o)


# ------------------------------------------------------- the in-art labels
def labels() -> str:
    """Every label the master carried, re-set at the page's type scale.

    Positions are the master's own anchors except where the bigger type needs
    the room: the two tag names are pushed apart, and the ones that were
    centred by eye are centred on the thing they name.
    """
    o = []

    # The two tags on the probe diagram.  The master sets these at 19 units
    # and tucks them either side of the bars; at 30 there is only one clear
    # pocket left, the white block under the orange bar (sheet x 415..510,
    # y 300..345), bounded by the probe's own stem on the left and the next
    # receptor below.  SpyTag goes there and ALFATag keeps the master's
    # position, which is already clear at this size.
    x, y = art(471, 252)
    o.append(text(x, y, "SpyTag", size=TAG, fill=AMBER_DK, weight=700))
    x, y = art(616, 252)
    o.append(text(x, y, "ALFATag", size=TAG, fill=GREEN_DK, weight=700))

    # The barcode, centred on the nucleus it sits in (sheet x 354) rather
    # than on the master's text box, which started 15 units to the right.
    x, y = art(415, 503)
    o.append(text(x, y, "barcode", size=MED, fill=INK, weight=500,
                  anchor="middle"))

    # the condition, under the first arrow
    x, y = art(714, 508)
    o.append(text(x, y, "37 °C", size=MED, fill=INK, weight=600,
                  anchor="middle"))

    # the cut, beside the scissors the master already draws
    x, y = art(1024, 74)
    o.append(lines(x, y, ("Reducing agent", "or protease"), size=MED,
                   weight=500))

    # the axis of the readout column
    x, y = art(1638, 664)
    o.append(f'<g transform="translate({x:.1f},{y:.1f}) rotate(-90)">'
             + text(0, 0, "Normalized internalization", size=MED, fill=INK,
                    weight=600) + "</g>")

    # What happens to the sorted cells.  Raised and pushed right off the
    # master's anchor: at this size the third line ran through the right-hand
    # sort arrow, and the block now sits in the gap between the readout cells
    # and the fork (sheet y 790..880, arrow at 890).
    x, y = art(1841, 706)
    o.append(lines(x, y, ("Sort and", "sequence", "barcodes"), size=SML,
                   weight=500))
    return "".join(o)


# -------------------------------------------------------------------- build
def build() -> str:
    # The master embeds its glows as raster <image> with xlink:href, so the
    # sheet has to declare that namespace or nothing in it renders at all.
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" '
         f'xmlns:xlink="http://www.w3.org/1999/xlink" '
         f'viewBox="0 0 {W:.0f} {H:.0f}" width="{W:.0f}" height="{H:.0f}">',
         f'<rect width="{W:.0f}" height="{H:.0f}" fill="#ffffff"/>']

    o.append(key_strip())
    o.append(f'<g transform="translate({-CROP_X:.1f},'
             f'{-CROP_Y + KEY_H:.1f})">{ART.body}</g>')
    o.append(labels())

    o.append("</svg>")
    return "".join(o)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    svg = build()
    svg_path = OUT_DIR / f"{STEM}.svg"
    svg_path.write_text(svg, encoding="utf-8")

    # qlmanage only ever writes square thumbnails, so rasterise a padded
    # square copy and crop the band back out with sips.
    side = 3120
    scale = side / W
    inner = (f'<g transform="scale({scale:.6f}) translate(0,'
             f'{(W - H) / 2:.2f})">')
    square = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'viewBox="0 0 {side} {side}" '
        f'width="{side}" height="{side}"><rect width="{side}" height="{side}" '
        f'fill="#ffffff"/>{inner}'
        + svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
        + "</g></svg>"
    )
    tmp = OUT_DIR / "_square.svg"
    tmp.write_text(square, encoding="utf-8")
    subprocess.run(["qlmanage", "-t", "-s", str(side), "-o", str(OUT_DIR),
                    str(tmp)], check=True, capture_output=True)
    png = OUT_DIR / f"{STEM}.png"
    png.unlink(missing_ok=True)
    (OUT_DIR / "_square.svg.png").replace(png)
    tmp.unlink()
    subprocess.run(["sips", "-c", str(round(side * H / W)), str(side),
                    str(png)], check=True, capture_output=True)
    pt = 468.0 / W
    print(f"wrote {svg_path}  ({W:.0f} x {H:.0f}, place at 6.5 in wide)")
    print(f"wrote {png}  ({side} px wide, ~{side / 6.5:.0f} dpi at 6.5 in)")
    print(f"type at 6.5 in: {BIG * pt:.1f} / {MED * pt:.1f} / "
          f"{SML * pt:.1f} pt, tag names {TAG * pt:.1f} pt")


if __name__ == "__main__":
    main()
