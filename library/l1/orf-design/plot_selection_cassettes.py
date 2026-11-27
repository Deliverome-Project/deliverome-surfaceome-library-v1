"""Selection cassettes side by side, for the promoter/selection slide of the 2026-11 SAB deck.

Five TFRC test constructs, read from their ordered maps on S3 (``constructs/ordered/``):
P-0003, P-0011 and P-0012 carry an EFS / EF1a / PGK promoter driving the ORF, then
IRES2 -> mScarlet3-P2A-PuroR on the same mRNA; P-0009 is P-0003 with the original EMCV IRES2
instead of the Golden Gate-domesticated one; P-0008 drops the IRES and mScarlet3 and drives
PuroR from its own SV40 promoter.

"Domesticated" means the IRES2 has had its internal BbsI sites removed by single-base edits
(Fragmid pF_AA333; s1e2 README §4.9), so BbsI cuts the vector only at the assembly junctions.
The script checks this on the maps: zero BbsI sites in P-0003's IRES2, two in P-0009's.

Promoter colours match the S5E11 slides (EFS = primary, EF1a = secondary, PGK = accent), so
the ORF and the other cassette parts take colours that none of the three promoters use.

Box widths are schematic (not to scale) so every element can carry a legible label; what
is read from the maps is which elements each plasmid has and their order, both asserted.

Rows run in the slide's reveal order: PGK (the dimmest), EF1a, EFS, EFS with the original IRES2,
then SV40-puro. ``--layers`` also writes one cumulative frame per row (``selection_cassettes_step{k}``)
for an on-click build; later rows are drawn fully transparent, so every frame has the same canvas.

    uv run --locked python plot_selection_cassettes.py [--layers]
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.patches import FancyArrowPatch, Rectangle

from deliverome_analysis.core import s3
from deliverome_analysis.core.plotting import (
    COLORS,
    SEQUENTIAL_PALETTES,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent
INK = COLORS["dark"]
#: promoter -> colour, the same COLORS keys as S5E11's construct colours
#: plus S5E11's light maroon for P-0009, which its plots colour apart from P-0003
PROMOTER = {
    "EFS": COLORS["primary"],
    "EFS, original IRES": SEQUENTIAL_PALETTES["maroon"][4],
    "EF1a": COLORS["secondary"],
    "PGK": COLORS["accent"],
}
C = {
    "sv40": COLORS["neutral"],
    "orf": COLORS["light"],
    "ires": SEQUENTIAL_PALETTES["lavender"][5],
    "mscarlet": COLORS["quaternary"],
    "p2a": COLORS["line"],
    "puro": COLORS["dark"],
}
#: boxes drawn with white text
DARK = {"EFS", "EF1a", "sv40", "mscarlet", "puro"}
BBSI = "GAAGAC"  # palindrome-free; reverse complement GTCTTC
# (map key, row label, [(GenBank label, box text, colour key, width)]) in plasmid order
ORF = ("TFRC_GG_domesticated", "ORF · tag · BC", "orf", 23)
PROMOTER_W = 11
ROWS = [
    ("p012", "P-0012", "PGK_v1", "PGK"),
    ("p011", "P-0011", "EF1a_v1", "EF1a"),
    ("p003", "P-0003", "EFS_v1", "EFS"),
    ("p009", "P-0009", "EFS_v1", "EFS"),
    ("p008", "P-0008", "EFS_v1", "EFS"),
]
#: GenBank label of the IRES2 feature: domesticated ("_v1") everywhere except P-0009
IRES_LABEL = {"p009": "IRES2"}
IRES_UNIT = [
    ("IRES2_v1", "IRES2", "ires", 12),
    ("mScarlet3", "mScarlet3", "mscarlet", 18),
    ("P2A", "P2A", "p2a", 8),
    ("PuroR", "PuroR", "puro", 12),
]
SV40_UNIT = [("SV40 promoter", "SV40", "sv40", 12), ("PuroR", "PuroR", "puro", 12)]


def load(key: str):
    with tempfile.TemporaryDirectory() as tmp:
        return SeqIO.read(s3.processed.download(key, Path(tmp) / "m.gb"), "genbank")


def first_start(rec, label: str) -> int:
    starts = [
        int(f.location.start) for f in rec.features if f.qualifiers.get("label", [""])[0] == label
    ]
    if not starts:
        raise SystemExit(f"{rec.id}: no feature labelled {label!r}")
    return min(starts)


def bbsi_sites(rec, label: str) -> int:
    """BbsI sites (either strand) inside the features carrying ``label``."""
    spans = [f.location for f in rec.features if f.qualifiers.get("label", [""])[0] == label]
    seq = str(rec.seq[min(s.start for s in spans) : max(s.end for s in spans)]).upper()
    return seq.count(BBSI) + seq.count("GTCTTC")


def promoter_arrow(ax, x, y):
    ax.plot([x, x], [y, y + 1.5], color=INK, lw=1.3)
    ax.add_patch(
        FancyArrowPatch(
            (x, y + 1.5), (x + 2.6, y + 1.5), arrowstyle="-|>", mutation_scale=10, color=INK, lw=1.3
        )
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "figures")
    ap.add_argument("--layers", action="store_true", help="also write one cumulative frame per row")
    args = ap.parse_args(argv)

    setup_plotting_style()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 64)
    ax.axis("off")

    h, x0 = 3.6, 16.0
    ys = [57.0, 48.0, 39.0, 29.0, 6.0]
    rows = []  # each row's artists, for --layers
    for (key, name, prom_label, prom_text), y in zip(ROWS, ys, strict=True):
        before = set(ax.get_children())
        rec = load(f"constructs/ordered/{key}.gb")
        ires_label = IRES_LABEL.get(key, "IRES2_v1")
        unit = SV40_UNIT if key == "p008" else [(ires_label, *IRES_UNIT[0][1:]), *IRES_UNIT[1:]]
        prom_colour = "EFS, original IRES" if key == "p009" else prom_text
        parts = [(prom_label, prom_text, prom_colour, PROMOTER_W), ORF, *unit]
        # the map must carry every element, in this order, or the row is not drawn
        starts = [first_start(rec, p[0]) for p in parts]
        assert starts == sorted(starts), (key, starts)
        labels = {f.qualifiers.get("label", [""])[0] for f in rec.features}
        if key == "p008":
            assert not labels & {"IRES2_v1", "IRES2", "mScarlet3"}, key
        elif key == "p009":
            # the original IRES2: not the domesticated feature, and its two BbsI sites intact
            assert "IRES2_v1" not in labels and bbsi_sites(rec, "IRES2") == 2, key
        else:
            assert bbsi_sites(rec, "IRES2_v1") == 0, key

        ax.text(0.5, y + h / 2, name, ha="left", va="center", fontsize=13, fontweight="bold")
        x = x0
        for label, text, ck, w in parts:
            face = PROMOTER.get(ck, C.get(ck))
            original = label == "IRES2"
            ax.add_patch(Rectangle((x, y), w, h, facecolor=face, edgecolor=INK, lw=0.8))
            if original:  # dashed frame around the one element that differs from P-0003
                ax.add_patch(Rectangle((x, y - 1.0), w, h + 2.0, facecolor="none",
                                       edgecolor=INK, lw=2.0, ls="--", zorder=3))  # fmt: skip
            ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=12,
                    color="white" if ck in DARK else INK)  # fmt: skip
            if ck in PROMOTER or ck == "sv40":
                promoter_arrow(ax, x, y + h)
            if original:
                ax.text(x + w / 2, y - 1.8, "original IRES2, not domesticated",
                        ha="center", va="top", fontsize=12, fontweight="bold", color=INK)  # fmt: skip
            x += w
        # how PuroR gets translated
        if key == "p008":
            ax.text(x0, y - 1.0, "PuroR from its own promoter:\nselection independent of the ORF mRNA",
                    ha="left", va="top", fontsize=12, style="italic", color=INK)  # fmt: skip
        if key == "p009":  # under the last IRES row, so it arrives with the IRES group complete
            ax.text(x0, y - 6.5, "PuroR rides on the ORF mRNA:\ntranslated only via the IRES",
                    ha="left", va="top", fontsize=12, style="italic", color=INK)  # fmt: skip
        rows.append(set(ax.get_children()) - before)

    save_figure(fig, "selection_cassettes", output_dir=args.out_dir)
    if args.layers:
        for k in range(1, len(rows) + 1):
            for i, artists in enumerate(rows):
                for artist in artists:
                    artist.set_alpha(1.0 if i < k else 0.0)
            save_figure(fig, f"selection_cassettes_step{k}", output_dir=args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
