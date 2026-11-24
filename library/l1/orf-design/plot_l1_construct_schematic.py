"""Linear schematic of the L1 library ORF construct (reviewer brief Figure 2a, SAB deck).

The barcoded insert between its two Esp3I shuttle sites is the subject (schematic, not to
scale), with the three tag-cassette positions in three colors and the L1 readout region
(``build_constructs.l1_readout``: stop, CC, T7 promoter, 5' handle JT832, barcode, Feldman's
padlock construct and RT site with T->A at 9, 10x CS1). Below it, to scale from cPPT to the
3' LTR, that insert in two vectors:

- the library layout, drawn first and in full color: EFS -> insert -> IRES2 -> PuroR, one
  transcript, so a cell survives puromycin only if it transcribes the ORF;
- the alternative not used, in grey: PuroR from its own SV40 promoter, no IRES (the P-0050
  layout). The TfR1 construct built this way, P-0008, gave puromycin-resistant cells with low
  TfR1 expression.

Two PNGs on one canvas, for an on-click build in the deck:

- ``l1_orf_construct.png`` -- the base figure, white background.
- ``l1_orf_construct_internal.png`` -- transparent overlay: the internal site outlined with a
  callout, and the N- and C-terminal sites faded. Every base artist is drawn into the
  overlay too, invisible, so both crop to the same bounding box and stack pixel-for-pixel;
  the script refuses to finish if their sizes differ.

Every coordinate is read from the worked-example GenBank maps on S3 that
``l1_example_vector_maps.py`` builds, never typed in: the AXL L1 insert in the library layout
(``constructs/scratch/l1-example-axl-ires-puro.gb``) and in the SV40-puro layout
(``constructs/scratch/l1-example-axl-sv40-puro.gb``). The readout blocks read from the map must
spell ``build_constructs.l1_readout`` exactly, the Esp3I overhangs are read back out of the
sequence, and the tag-site counts come from ``l1_tag_paths.tsv``. AXL is the example ORF; the
insert architecture is the same for every member, only the ORF and the cassette position change.

    uv run --locked python plot_l1_construct_schematic.py
"""

from __future__ import annotations

import argparse
import itertools
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from build_constructs import l1_readout
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from PIL import Image

from deliverome_analysis import constructs, s3
from deliverome_analysis.plotting import (
    COLORS,
    SEQUENTIAL_PALETTES,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
#: worked-example maps from l1_example_vector_maps.py; constructs.map_key finds the folder
MAPS = {"library": "l1-example-axl-ires-puro", "alternative": "l1-example-axl-sv40-puro"}
ESP3I = "CGTCTC"
BASE, OVERLAY = "l1_orf_construct", "l1_orf_construct_internal"

INK = COLORS["dark"]
GREY = COLORS["neutral"]
C = {
    "orf": COLORS["primary"],
    "barcode": COLORS["warning"],
    "primer": SEQUENTIAL_PALETTES["lavender"][5],
    "t7": SEQUENTIAL_PALETTES["amber"][5],
    "cs1": SEQUENTIAL_PALETTES["teal"][4],
    "esp3i": INK,
    "other": COLORS["line"],
}
#: one color per tag position, none of them used elsewhere in the insert
TAG = {"N": COLORS["accent"], "internal": COLORS["secondary"], "C": COLORS["quaternary"]}
#: the library-layout row, full color; the alternative row is all grey
ROW = {
    "EFS": (SEQUENTIAL_PALETTES["lavender"][4], INK),
    "insert": (COLORS["primary"], "white"),
    "IRES2": (SEQUENTIAL_PALETTES["teal"][5], INK),
    "PuroR": (COLORS["success"], "white"),
    "WPRE": (COLORS["light"], INK),
    "": (COLORS["light"], INK),
}
FS = 15  # insert text; nothing in the figure goes below 12 pt

#: the readout region after the ORF, in map order: (map label, box text, colour key)
READOUT = [
    ("stop", None, None),
    ("pad", None, None),
    ("T7 promoter (readout)", "T7", "t7"),
    ("5' handle (JT832)", "5' handle\n(JT832)", "primer"),
    ("barcode (20 N placeholder)", "barcode", "barcode"),
    ("padlock construct (Feldman, 17 nt)", "padlock", "primer"),
    ("RT site (Feldman, T9A)", "RT site\n(T9A)", "primer"),
    ("10x Capture Sequence 1", "10x CS1", "cs1"),
]


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
        ax.text((x0 + x1) / 2, y + h / 2, text, ha="center", va="center", color=tc, fontsize=fs,
                linespacing=1.0)  # fmt: skip


def promoter_arrow(ax, x, y, color=GREY):
    ax.plot([x, x], [y, y + 1.1], color=color, lw=1.0)
    ax.add_patch(
        FancyArrowPatch(
            (x, y + 1.1),
            (x + 1.8, y + 1.1),
            arrowstyle="-|>",
            mutation_scale=8,
            color=color,
            lw=1.0,
        )
    )


def nt(rec, label) -> int:
    s, e = span(rec, label)
    return e - s + 1


def check_maps(a, b) -> tuple[str, str, tuple[int, int]]:
    """Assert the two maps share one insert and the L1 readout; return its overhangs and span."""
    oh5, oh3 = overhangs(a)
    assert overhangs(b) == (oh5, oh3)
    # the insert is identical in both: everything from the 5' Esp3I to the 3' Esp3I
    ins_a = span(a, "Esp3I (shuttle out)", 0)[0], span(a, "Esp3I (shuttle out)", 1)[1]
    ins_b = span(b, "Esp3I (shuttle out)", 0)[0], span(b, "Esp3I (shuttle out)", 1)[1]
    assert ins_a == ins_b
    assert str(a.seq)[ins_a[0] - 1 : ins_a[1]] == str(b.seq)[ins_b[0] - 1 : ins_b[1]]
    assert "IRES2_v1" not in {f.qualifiers.get("label", [""])[0] for f in b.features}

    # the readout region, read from the map, must be exactly the L1 layout
    seq = str(a.seq).upper()
    spans = [span(a, lab) for lab, *_ in READOUT]
    assert all(e + 1 == s for (_, e), (s, _) in itertools.pairwise(spans))
    # the ORF keeps its own stop codon (AXL's is not necessarily l1_readout's TAA); the rest of
    # the readout must be l1_readout exactly
    stop_codon = seq[spans[0][0] - 1 : spans[0][1]]
    assert stop_codon in ("TAA", "TAG", "TGA"), stop_codon
    assert seq[spans[1][0] - 1 : spans[-1][1]] == l1_readout("N" * 20)[3:], "readout != l1_readout"
    return oh5, oh3, ins_a


def tag_counts() -> tuple[int, int, int, int]:
    """CDS counts per tag placement (N-terminal, internal, C-terminal) and in all."""
    paths = pd.read_csv(HERE / "l1_tag_paths.tsv", sep="\t")
    n = paths["placement"].value_counts()
    n_nterm = int(n[[k for k in n.index if k.startswith("N ")]].sum())
    n_int, n_cterm, n_all = int(n["internal"]), int(n["C-terminal"]), len(paths)
    assert n_nterm + n_int + n_cterm == n_all
    return n_nterm, n_int, n_cterm, n_all


def draw_vectors(ax, a, b, ins_a, zx0, zx1, y) -> None:
    """The two vectors to scale under the insert, and the zoom connector up to it."""
    start = span(a, "cPPT/CTS")[0]
    end_a = span(a, "3' LTR (Delta-U3)")[1]
    sx0, sx1 = zx0, zx1
    k = (sx1 - sx0) / (end_a - start)

    def X(bp):
        return sx0 + (bp - start) * k

    rh = 3.0
    rows = [
        (
            a,
            25.0,
            True,
            "Library layout (chosen): EFS \N{EN DASH} ORF \N{EN DASH} IRES2 \N{EN DASH} PuroR, one transcript",
        ),
        (b, 12.0, False, "Alternative, not used: PuroR from its own SV40 promoter, no IRES"),
    ]
    for rec, ry, chosen, name in rows:
        line = INK if chosen else GREY
        end = span(rec, "3' LTR (Delta-U3)")[1]
        ax.plot([X(start), X(end)], [ry + rh / 2] * 2, color=line, lw=0.8, zorder=0)
        ax.text(sx0, ry - 0.9, name, ha="left", va="top", fontsize=14 if chosen else 13,
                color=line, fontweight="bold" if chosen else "normal")  # fmt: skip
        feats = [("cPPT/CTS", ""), ("EFS_v1", "EFS")]
        if chosen:
            feats += [("IRES2_v1", "IRES2")]
        else:
            feats += [("SV40 enhancer", "SV40")]
        feats += [("PuroR", "PuroR"), ("WPRE", "WPRE"), ("3' LTR (Delta-U3)", "")]
        for label, text in feats:
            s, e = span(rec, label)
            if label == "IRES2_v1":  # annotated as three adjacent pieces; draw as one
                e = span(rec, label, 2)[1]
            if label == "SV40 enhancer":  # the enhancer and promoter are one block here
                e = span(rec, "SV40 promoter")[1]
            fc, tc = ROW[text if text in ROW else ""] if chosen else (COLORS["light"], GREY)
            box(ax, X(s - 1), X(e), ry, rh, fc, text, tc=tc, fs=13, ec=line, lw=0.6)
            if label == "SV40 enhancer":
                promoter_arrow(ax, X(s - 1), ry + rh, line)
        promoter_arrow(ax, X(span(rec, "EFS_v1")[0] - 1), ry + rh, line)
        i0, i1 = ins_a
        fc, tc = ROW["insert"] if chosen else (SEQUENTIAL_PALETTES["maroon"][5], GREY)
        box(ax, X(i0 - 1), X(i1), ry, rh, fc, "insert", tc=tc, fs=13, ec=line, lw=0.6)

    # what was seen with the alternative layout
    ax.text(sx0, 12.0 - 4.2,
            "Deliverome scientists saw TfR1 silencing in the SV40-puro TfR1 construct (P-0008):\n"
            "puromycin-resistant cells with low TfR1 expression.",
            ha="left", va="top", fontsize=13, color=INK, style="italic",
            linespacing=1.2)  # fmt: skip

    # zoom connector from the insert in the library layout to its expansion
    i0, i1 = X(ins_a[0] - 1), X(ins_a[1])
    ax.add_patch(Polygon([(i0, 25.0 + rh + 0.2), (i1, 25.0 + rh + 0.2), (zx1, y - 6.6), (zx0, y - 6.6)],
                         closed=True, facecolor=COLORS["line"], edgecolor="none", alpha=0.45,
                         zorder=-1))  # fmt: skip


def overlay_artists(ax, sites, y, h, tw, n_int) -> list[tuple]:
    """Add the internal-site overlay's artists; return them with their shown alpha."""
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
    over.append((FancyArrowPatch((cx + 3.0, y - 8.0), (cx + 0.4, y - 1.6), arrowstyle="-|>",
                                 mutation_scale=16, color=TAG["internal"], lw=2.2,
                                 zorder=7), 1.0))  # fmt: skip
    for art, _ in over:
        ax.add_patch(art)
    callout = ax.text(cx + 9.0, y - 9.8, f"{n_int:,} members: loop insertion", ha="center",
                      va="center", fontsize=15, fontweight="bold", color="white", zorder=7,
                      bbox={"boxstyle": "round,pad=0.45", "facecolor": TAG["internal"],
                            "edgecolor": "none"})  # fmt: skip
    over += [(callout, 1.0), (callout.get_bbox_patch(), 1.0)]
    return over


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "figures")
    args = ap.parse_args(argv)

    recs = {name: load(constructs.map_key(m)) for name, m in MAPS.items()}
    a, b = recs["library"], recs["alternative"]
    oh5, oh3, ins_a = check_maps(a, b)
    n_nterm, n_int, n_cterm, n_all = tag_counts()

    setup_plotting_style()
    fig, ax = plt.subplots(figsize=(13, 6.4))
    fig.subplots_adjust(0, 0, 1, 1)
    ax.set_xlim(0, 100)
    ax.set_ylim(3, 67)
    ax.axis("off")

    # ---- the insert: the subject of the figure, schematic ----------------------------
    y, h = 41.0, 7.0
    stop_nt = nt(a, "stop") + nt(a, "pad")
    parts = [
        # (key, width, label in box, length label, colour key)
        ("esp5", 5.0, "Esp3I", "", "esp3i"),
        ("kozak", 6.5, "Kozak", f"{nt(a, 'Kozak (AgeI-XbaI)')} nt", "other"),
        ("orf", 24.5, "", "", "orf"),
        ("stop", 6.5, "stop +\nCC pad", f"{stop_nt} nt", "other"),
    ]
    widths = {"T7": 5.0, "5' handle\n(JT832)": 8.5, "barcode": 8.5, "padlock": 7.5,
              "RT site\n(T9A)": 7.5, "10x CS1": 7.5}  # fmt: skip
    for lab, text, ck in READOUT[2:]:
        parts.append((lab, widths[text], text, f"{nt(a, lab)} nt", ck))
    parts.append(("esp3", 5.0, "Esp3I", "", "esp3i"))
    x = 100 - sum(p[1] for p in parts) - 1.0
    pos = {}
    for key, w, label, length, ck in parts:
        dark = ck in ("esp3i", "orf", "barcode")
        fs = 12 if key == "stop" else (13 if "\n" in label else 14)
        box(ax, x, x + w, y, h, C[ck], label, tc="white" if dark else INK, fs=fs, lw=1.0)
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

    # tag placement per CDS, three possible positions, three colours
    ox0, ox1 = pos["orf"]
    # left of center, so the overlay's arrow to the internal site does not cross it
    ax.text(ox0 + 0.3 * (ox1 - ox0), y - 0.8, "ORF (variable)", ha="center", va="top", fontsize=13,
            color=INK)  # fmt: skip
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
    ax.text(ox0, y + h + 7.8, f"tag placement per CDS, {n_all:,} CDSs (an internal CDS may have 2-3 sites)",
            ha="left", va="bottom", fontsize=13, color=INK)  # fmt: skip

    # bracket: what moves in one Golden Gate step
    by = 63.0
    ax.plot([zx0, zx0, zx1, zx1], [by - 1, by, by, by - 1], color=INK, lw=1)
    ax.text((zx0 + zx1) / 2, by + 0.4, "barcoded ORF insert, released intact by Esp3I (BsmBI)",
            ha="center", va="bottom", fontsize=FS, color=INK)  # fmt: skip

    # ---- the two vectors, to scale: the library layout first, the alternative grey ----
    draw_vectors(ax, a, b, ins_a, zx0, zx1, y)

    # ---- overlay: the internal site, for the on-click build ---------------------------
    base_artists = list(ax.patches) + list(ax.texts) + list(ax.lines)
    over = overlay_artists(ax, sites, y, h, tw, n_int)

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
    print("  readout: " + " | ".join(f"{lab} {nt(a, lab)} nt" for lab, *_ in READOUT))
    print(f"  library layout {len(a.seq):,} bp; alternative {len(b.seq):,} bp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
