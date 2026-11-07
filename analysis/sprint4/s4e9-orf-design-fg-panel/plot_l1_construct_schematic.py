"""Linear schematic of the L1 library ORF construct, for the 2026-11 SAB deck.

The barcoded insert between its two Esp3I shuttle sites is the subject (schematic, not to
scale), with the three tag-cassette positions in three colours. Below it, small and grey, is
that insert in the two selection backbones, to scale from cPPT to the 3' LTR: P-0039
(EFS -> insert -> IRES2 -> PuroR, one transcript) and P-0050 (PuroR from its own SV40
promoter).

Two PNGs on one canvas, for an on-click build in the deck:

- ``l1_orf_construct.png`` -- the base figure, white background.
- ``l1_orf_construct_internal.png`` -- transparent overlay: the internal site outlined with a
  callout, and the N- and C-terminal sites faded. Every base artist is drawn into the
  overlay too, invisible, so both crop to the same bounding box and stack pixel-for-pixel;
  the script refuses to finish if their sizes differ.

Every coordinate is read from the designed GenBank maps on S3, never typed in: AXL in
P-0039 (``constructs/designed/p0024.gb``) and the same insert in P-0050
(``constructs/designed/p0053.gb``). The Esp3I overhangs are read back out of the sequence
and the tag-site counts come from ``l1_tag_paths.tsv`` (2,700 L1 members). AXL is the
example ORF; the insert architecture is the same for every member, only the ORF and the
cassette position change.

    uv run --locked python plot_l1_construct_schematic.py
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from PIL import Image

from deliverome_analysis import s3
from deliverome_analysis.plotting import (
    COLORS,
    SEQUENTIAL_PALETTES,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
MAPS = {"P-0039": "constructs/designed/p0024.gb", "P-0050": "constructs/designed/p0053.gb"}
ESP3I = "CGTCTC"
BASE, OVERLAY = "l1_orf_construct", "l1_orf_construct_internal"

INK = COLORS["dark"]
GREY = COLORS["neutral"]
C = {
    "orf": COLORS["primary"],
    "barcode": COLORS["warning"],
    "primer": SEQUENTIAL_PALETTES["lavender"][5],
    "cs1": SEQUENTIAL_PALETTES["teal"][4],
    "esp3i": INK,
    "other": COLORS["line"],
}
#: one colour per tag position, none of them used elsewhere in the insert
TAG = {"N": COLORS["accent"], "internal": COLORS["secondary"], "C": COLORS["quaternary"]}
FS = 15  # insert text; nothing in the figure goes below 12 pt


def load(key: str):
    with tempfile.TemporaryDirectory() as tmp:
        return SeqIO.read(s3.processed.download(key, Path(tmp) / "m.gb"), "genbank")


def span(rec, label: str, nth: int = 0) -> tuple[int, int]:
    """1-based inclusive (start, end) of the nth feature with this label."""
    hits = [f for f in rec.features if f.qualifiers.get("label", [""])[0] == label]
    if len(hits) <= nth:
        raise SystemExit(f"{rec.id}: no feature labelled {label!r}")
    f = hits[nth]
    return int(f.location.start) + 1, int(f.location.end)


def overhangs(rec) -> tuple[str, str]:
    """The 4-nt overhangs the two inward-facing Esp3I sites leave, read from the sequence."""
    seq = str(rec.seq).upper()
    (s5, e5), (s3_, e3) = span(rec, "Esp3I (shuttle out)", 0), span(rec, "Esp3I (shuttle out)", 1)
    up, down = seq[s5 - 1 : e5], seq[s3_ - 1 : e3]
    assert up[:6] == ESP3I and down[-6:] == "GAGACG", (up, down)
    # CGTCTC(1/5): one spacer base, then the overhang; mirrored for the reverse site
    oh5 = seq[e5 : e5 + 4]
    oh3 = seq[e3 - 7 - 4 : e3 - 7]
    return oh5, oh3


def box(ax, x0, x1, y, h, color, text="", tc="white", fs=FS, ec=INK, lw=0.8, **kw):
    ax.add_patch(Rectangle((x0, y), x1 - x0, h, facecolor=color, edgecolor=ec, lw=lw, **kw))
    if text:
        ax.text((x0 + x1) / 2, y + h / 2, text, ha="center", va="center", color=tc, fontsize=fs)


def promoter_arrow(ax, x, y):
    ax.plot([x, x], [y, y + 1.1], color=GREY, lw=1.0)
    ax.add_patch(
        FancyArrowPatch(
            (x, y + 1.1), (x + 1.8, y + 1.1), arrowstyle="-|>", mutation_scale=8, color=GREY, lw=1.0
        )
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "figures")
    args = ap.parse_args(argv)

    recs = {name: load(key) for name, key in MAPS.items()}
    a, b = recs["P-0039"], recs["P-0050"]
    oh5, oh3 = overhangs(a)
    assert overhangs(b) == (oh5, oh3)
    # the insert is identical in both: everything from the 5' Esp3I to the 3' Esp3I
    ins_a = span(a, "Esp3I (shuttle out)", 0)[0], span(a, "Esp3I (shuttle out)", 1)[1]
    ins_b = span(b, "Esp3I (shuttle out)", 0)[0], span(b, "Esp3I (shuttle out)", 1)[1]
    assert ins_a == ins_b
    assert str(a.seq)[ins_a[0] - 1 : ins_a[1]] == str(b.seq)[ins_b[0] - 1 : ins_b[1]]

    paths = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t")
    n = paths["placement"].value_counts()
    n_nterm = int(n[[k for k in n.index if k.startswith("N ")]].sum())
    n_int, n_cterm, n_all = int(n["internal"]), int(n["C-terminal"]), len(paths)
    assert n_nterm + n_int + n_cterm == n_all

    def nt(rec, label):
        s, e = span(rec, label)
        return e - s + 1

    setup_plotting_style()
    fig, ax = plt.subplots(figsize=(13, 5))
    fig.subplots_adjust(0, 0, 1, 1)  # 10 data units per inch on both axes
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.axis("off")

    # ---- the insert: the subject of the figure, schematic ----------------------------
    y, h = 25.0, 7.0
    parts = [
        # (key, width, label in box, length label, colour key)
        ("esp5", 6.0, "Esp3I", "", "esp3i"),
        ("kozak", 6.5, "Kozak", f"{nt(a, 'Kozak (AgeI-XbaI)')} nt", "other"),
        ("orf", 26.5, "", "", "orf"),
        ("stop", 5.0, "stop", "", "other"),
        ("handle", 9.0, "5' handle", f"{nt(a, "5' handle")} nt", "primer"),
        ("barcode", 8.5, "barcode", f"{nt(a, 'barcode')} nt", "barcode"),
        ("rt", 7.5, "RT site", f"{nt(a, 'gap + RT site')} nt", "primer"),
        ("cs1", 8.0, "10x CS1", f"{nt(a, '10x Capture Sequence 1')} nt", "cs1"),
        ("esp3", 6.0, "Esp3I", "", "esp3i"),
    ]
    x = 16.0
    pos = {}
    for key, w, label, length, ck in parts:
        dark = ck in ("esp3i", "orf", "barcode")
        box(ax, x, x + w, y, h, C[ck], label, tc="white" if dark else INK, lw=1.0)
        if length:
            ax.text(x + w / 2, y - 0.8, length, ha="center", va="top", fontsize=13, color=INK)
        pos[key] = (x, x + w)
        x += w
    zx0, zx1 = pos["esp5"][0], pos["esp3"][1]

    # Esp3I overhangs, read from the map
    ax.text(pos["esp5"][0], y - 3.6, f"→ {oh5}", ha="left", va="top", fontsize=13,
            family="monospace", color=INK, fontweight="bold")  # fmt: skip
    ax.text(pos["esp3"][1], y - 3.6, f"{oh3} ←", ha="right", va="top", fontsize=13,
            family="monospace", color=INK, fontweight="bold")  # fmt: skip

    # tag cassette: one site per member, three possible positions, three colours
    ox0, ox1 = pos["orf"]
    ax.text(ox0 + 0.27 * (ox1 - ox0), y - 0.8, "ORF (variable)", ha="center", va="top",
            fontsize=13, color=INK)  # fmt: skip
    sites = {
        "N": (ox0 + 3.5, "N-term", n_nterm),
        "internal": ((ox0 + ox1) / 2, "internal", n_int),
        "C": (ox1 - 3.5, "C-term", n_cterm),
    }
    tw = 2.6
    site_artists: dict[str, list] = {}
    for k, (cx, name, count) in sites.items():
        before = len(ax.patches), len(ax.texts)
        box(ax, cx - tw / 2, cx + tw / 2, y - 0.6, h + 1.2, TAG[k], lw=1.0)
        ax.text(cx, y + h + 1.2, f"{name}\n{count:,}", ha="center", va="bottom", fontsize=14,
                color=TAG[k] if k != "N" else COLORS["warning"], fontweight="bold",
                linespacing=1.15)  # fmt: skip
        site_artists[k] = ax.patches[before[0] :] + ax.texts[before[1] :]
    ax.text(ox0, y + h + 7.8, f"tag cassette: one site per member, {n_all:,} members",
            ha="left", va="bottom", fontsize=13, color=INK)  # fmt: skip

    # bracket: what moves in one Golden Gate step
    by = 47.0
    ax.plot([zx0, zx0, zx1, zx1], [by - 1, by, by, by - 1], color=INK, lw=1)
    ax.text((zx0 + zx1) / 2, by + 0.4, "barcoded ORF insert, released intact by Esp3I (BsmBI)",
            ha="center", va="bottom", fontsize=FS, color=INK)  # fmt: skip

    # ---- the two backbones: small and grey, to scale ----------------------------------
    start = span(a, "cPPT/CTS")[0]
    end_a = span(a, "3' LTR (Delta-U3)")[1]
    sx0, sx1 = 16.0, 98.0
    k = (sx1 - sx0) / (end_a - start)

    def X(bp):
        return sx0 + (bp - start) * k

    rows = [(a, 8.5, "P-0039  IRES-puro"), (b, 2.5, "P-0050  SV40-puro")]
    rh = 2.4
    face = COLORS["light"]
    for rec, ry, name in rows:
        end = span(rec, "3' LTR (Delta-U3)")[1]
        ax.plot([X(start), X(end)], [ry + rh / 2] * 2, color=GREY, lw=0.8, zorder=0)
        ax.text(0.5, ry + rh / 2, name, ha="left", va="center", fontsize=12, color=GREY)
        feats = [("cPPT/CTS", ""), ("EFS_v1", "EFS")]
        if rec is a:
            feats += [("IRES2_v1", "IRES2")]
        feats += [("PuroR", "PuroR"), ("WPRE", "WPRE"), ("3' LTR (Delta-U3)", "")]
        for label, text in feats:
            s, e = span(rec, label)
            if label == "IRES2_v1":  # annotated as three adjacent pieces; draw as one
                e = span(rec, label, 2)[1]
            if label == "PuroR" and rec is b:  # ATG annotated separately from the body
                e = span(rec, label, 1)[1]
            box(ax, X(s - 1), X(e), ry, rh, face, text, tc=GREY, fs=12, ec=GREY, lw=0.6)
        promoter_arrow(ax, X(span(rec, "EFS_v1")[0] - 1), ry + rh)
        if rec is b:
            # the enhancer and promoter are one block on the slide
            s, e = span(rec, "SV40 enhancer")[0], span(rec, "SV40 promoter")[1]
            box(ax, X(s - 1), X(e), ry, rh, face, "SV40", tc=GREY, fs=12, ec=GREY, lw=0.6)
            promoter_arrow(ax, X(s - 1), ry + rh)
        # the insert, as one block, tinted to tie it to the expansion above
        i0, i1 = ins_a
        box(ax, X(i0 - 1), X(i1), ry, rh, SEQUENTIAL_PALETTES["maroon"][5], "insert",
            tc=GREY, fs=12, ec=GREY, lw=0.6)  # fmt: skip

    # zoom connector from the insert in P-0039 to its expansion
    i0, i1 = X(ins_a[0] - 1), X(ins_a[1])
    ax.add_patch(Polygon([(i0, 8.5 + rh + 0.2), (i1, 8.5 + rh + 0.2), (zx1, 20.0), (zx0, 20.0)],
                         closed=True, facecolor=COLORS["line"], edgecolor="none", alpha=0.45,
                         zorder=-1))  # fmt: skip

    # ---- overlay: the internal site, for the on-click build ---------------------------
    base_artists = list(ax.patches) + list(ax.texts) + list(ax.lines)
    cx = sites["internal"][0]
    over: list[tuple] = []  # (artist, its alpha when shown)
    # fade the two terminal sites back toward the ORF so only the internal one stands out
    for kk in ("N", "C"):
        sx = sites[kk][0]
        over.append((Rectangle((sx - tw / 2 - 0.1, y - 0.7), tw + 0.2, h + 1.4,
                               facecolor=C["orf"], edgecolor="none", zorder=5), 0.75))  # fmt: skip
        over.append((Rectangle((sx - 5.0, y + h + 0.9), 10.0, 6.6, facecolor="white",
                               edgecolor="none", zorder=5), 0.8))  # fmt: skip
    over.append((Rectangle((cx - tw / 2 - 0.6, y - 1.4), tw + 1.2, h + 2.8, facecolor="none",
                           edgecolor=INK, lw=2.5, zorder=6), 1.0))  # fmt: skip
    over.append((FancyArrowPatch((cx + 3.0, 17.3), (cx + 0.4, y - 1.6), arrowstyle="-|>",
                                 mutation_scale=16, color=TAG["internal"], lw=2.2,
                                 zorder=7), 1.0))  # fmt: skip
    for art, _ in over:
        ax.add_patch(art)
    callout = ax.text(cx + 9.0, 15.4, f"{n_int:,} members: loop insertion", ha="center",
                      va="center", fontsize=15, fontweight="bold", color="white", zorder=7,
                      bbox={"boxstyle": "round,pad=0.45", "facecolor": TAG["internal"],
                            "edgecolor": "none"})  # fmt: skip
    over += [(callout, 1.0), (callout.get_bbox_patch(), 1.0)]

    # both PNGs draw every artist, so `bbox_inches="tight"` crops them identically; each
    # hides the other's artists with alpha 0 rather than removing them
    for art, _ in over:
        art.set_alpha(0)
    save_figure(fig, BASE, output_dir=args.out_dir)
    for art in base_artists:
        art.set_alpha(0)
    for art, alpha in over:
        art.set_alpha(alpha)
    save_figure(fig, OVERLAY, output_dir=args.out_dir, transparent=True)

    sizes = {f: Image.open(args.out_dir / f"{f}.png").size for f in (BASE, OVERLAY)}
    if len(set(sizes.values())) != 1:
        raise SystemExit(f"base and overlay differ in size: {sizes}")
    print(f"  {sizes[BASE][0]} x {sizes[BASE][1]} px each")
    print(f"  overhangs {oh5} / {oh3}; tag sites N {n_nterm}, internal {n_int}, C {n_cterm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
