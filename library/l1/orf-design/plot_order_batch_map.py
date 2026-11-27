"""One page showing the block layout of every insert in the p0022-p0047 order.

For cross-checking a batch against an external registry (Benchling) before ordering:
each row is one construct drawn to scale, so a wrong cassette position, a missing
barcode or an off-by-one ORF split is visible without opening 23 GenBank files.

Blocks come from ``annotate_insert.blocks()``, which locates every block by SEARCHING
for its own sequence in the built insert and asserts the spans tile it exactly. So this
figure cannot drift from what was built -- if a construct's layout changed, the drawing
changes with it. Nothing here is a hand-typed coordinate.

    uv run --locked python plot_order_batch_map.py
    uv run --locked python plot_order_batch_map.py --out-dir ~/Downloads/...
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from annotate_insert import SOURCES, blocks
from build_constructs import read_fasta

from deliverome_analysis.core.plotting import (
    CATEGORICAL_PALETTE,
    create_figure,
    get_color,
    save_figure,
    setup_plotting_style,
)

HERE = Path(__file__).parent

#: Which barcode table holds which arm, and whether rows are gated on a column.
TABLES = (
    ("barcodes.tsv", None),
    ("barcodes_loop.tsv", None),
    ("barcodes_tfrc.tsv", "in_order_batch"),
    ("barcodes_gpcr.tsv", None),
    ("barcodes_sp_swap.tsv", None),
)

#: Block label -> palette role. Grouped by what the block IS, not by sequence order:
#: everything that encodes the receptor is one colour, the tag another, the amplicon
#: readout a third, so a glance separates "my protein" from "the machinery".
ROLE = {
    "Kozak (AgeI-XbaI)": "regulatory",
    "ORF, N-terminal": "orf",
    "ORF, remainder": "orf",
    "stop": "regulatory",
    "pad": "regulatory",
    "5' handle": "readout",
    "barcode": "barcode",
    "gap + RT site": "readout",
    "10x Capture Sequence 1": "readout",
}


#: (group, one-line reason this construct is in the batch). Keyed by plasmid, because
#: several genes appear more than once with different reasons -- SLC2A4 twice for two
#: loops, CCR5 three times across a loop and a leader pair.
GROUPS = {
    "p0022": (
        "Single-pass type I, native SP",
        "size ceiling — the deliberate test against the 5-kb ORF cap",
    ),
    "p0023": (
        "Single-pass type I, native SP",
        "slow internaliser: the low comparator the GPCR arms lack",
    ),
    "p0024": (
        "Single-pass type I, native SP",
        "multi-route — should resist full block by any single inhibitor",
    ),
    "p0027": (
        "Single-pass type I, native SP",
        "the only curated macropinocytosis primary — gives EIPA a target",
    ),
    "p0029": (
        "Single-pass type I, native SP",
        "added after the S2E5 selection; no role recorded in that table",
    ),
    "p0026": (
        "Native SP vs Igκ swap (pair)",
        "highest-novelty lead; keeps its own SP, which the cassette ABOLISHES",
    ),
    "p0048": (
        "Native SP vs Igκ swap (pair)",
        "same gene, SP 1-30 replaced by Igκ — rescues the call (p(SP) 0.9998)",
    ),
    "p0025": ("Multipass with a native SP", "raft route in a transporter, with a cAMP rate switch"),
    "p0028": ("GPI-anchored", "GPI-anchored: no cytoplasmic tail, so no sorting motif can operate"),
    "p0037": ("Single-pass type II", "ceiling control; already in hand as p003"),
    "p0030": (
        "Multipass, internal loop",
        "GLUT4: the one route separable pharmacologically (CME + cholesterol)",
    ),
    "p0031": (
        "Multipass, internal loop",
        "same gene, second loop — position contrast 252 residues away",
    ),
    "p0032": (
        "Multipass, internal loop",
        "floor control: signal here is assay background. Loop site, 5 aa from TM",
    ),
    "p0033": (
        "Multipass, internal loop",
        "floor control, same loop as p0032 — 20 aa from TM rather than 5",
    ),
    "p0034": ("GPCR, internal loop", "CCR5 ECL3; position contrast vs p0042"),
    "p0035": ("GPCR, internal loop", "AVPR2 ECL3; position contrast vs p0047"),
    "p0038": ("GPCR N-term, leader pair", "CNR1 leader-free; weakest expressor, no tether"),
    "p0041": ("GPCR N-term, leader pair", "CNR1 + Igκ; does a leader rescue it?"),
    "p0042": ("GPCR N-term, leader pair", "CCR5 leader-free; C20-C269 tether"),
    "p0043": ("GPCR N-term, leader pair", "CCR5 + Igκ; attributes a tether failure"),
    "p0044": ("GPCR N-term, leader pair", "OR51E2 leader-free; C4-C168 tether"),
    "p0045": ("GPCR N-term, leader pair", "OR51E2 + Igκ; can a leader rescue an OR?"),
    "p0046": ("GPCR N-term, leader-free only", "ACKR3; C34-C287 tether, NO Igκ twin"),
    "p0047": ("GPCR N-term, leader-free only", "AVPR2; strongest expressor, no tether"),
}
#: Draw order, coarse architecture first so like sits with like.
GROUP_ORDER = [
    "Single-pass type I, native SP",
    "Native SP vs Igκ swap (pair)",
    "Multipass with a native SP",
    "GPI-anchored",
    "Single-pass type II",
    "Multipass, internal loop",
    "GPCR, internal loop",
    "GPCR N-term, leader pair",
    "GPCR N-term, leader-free only",
]


def _rows() -> list[dict]:
    out: list[dict] = []
    for table, gate in TABLES:
        path = HERE / table
        if not path.exists():
            continue
        with path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if gate and r.get(gate) != "1":
                    continue
                out.append(r)
    return out


def _sites() -> dict[tuple[str, str], str]:
    """(gene, variant) -> insert_after_residue, for labelling the loop arms."""
    path = HERE / "constructs_loop.tsv"
    if not path.exists():
        return {}
    with path.open() as fh:
        return {
            (r["gene"], r["variant"]): r["insert_after_residue"]
            for r in csv.DictReader(fh, delimiter="\t")
        }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "figures")
    args = ap.parse_args(argv)

    inserts: dict[str, str] = {}
    for name in SOURCES:
        if (HERE / name).exists():
            inserts.update(read_fasta(HERE / name))
    sites = _sites()

    plan = []
    for r in sorted(_rows(), key=lambda x: x["plasmid"]):
        gene = r["gene"]
        variant = r.get("variant") or ""
        # Each table names its arm in a different column: `construct` for the GPCR
        # arms, `variant` for the loop arms, `arm` for TFRC, nothing for the N-term
        # eight. Try them in that order rather than assuming one schema.
        candidates = [
            r.get("construct"),
            f"{gene}_{variant}" if variant else None,
            f"{gene}_{r.get('arm')}" if r.get("arm") else None,
            gene,
        ]
        key = next((c for c in candidates if c and c in inserts), None)
        if key is None:
            raise SystemExit(
                f"{r['plasmid']}: no insert matched any of "
                f"{[c for c in candidates if c]!r}; have {sorted(inserts)[:6]}..."
            )
        site = sites.get((gene, variant))
        # `ORF, N-terminal` means three different things depending on the arm, so the
        # arm type has to be carried through: for an igk arm that block is the Igk
        # leader + Asp + the receptor's Met1 (66 nt), for a panel arm it is the native
        # signal peptide (87 nt on HEG1), and for a leader-free GPCR it is Met1 alone
        # (3 nt). Drawing all three the same colour is what hid the igk arms.
        if key.endswith("_igk"):
            kind, where = "igk", "N-term · Igκ leader"
        elif key.endswith("_plain"):
            kind, where = "plain", "N-term · leader-free"
        elif site:
            kind, where = "loop", f"loop, after {site}"
        elif gene == "TFRC":
            kind, where = "cterm", "C-term"
        else:
            kind, where = "nterm_sp", "N-term · native SP"
        grp, why = GROUPS[r["plasmid"]]
        plan.append((r["plasmid"], gene, where, inserts[key], kind, grp, why))

    # like with like: architecture group first, plasmid order within it
    plan.sort(key=lambda row: (GROUP_ORDER.index(row[5]), row[0]))

    setup_plotting_style(font_scale=1.0)
    fig, ax = create_figure(figsize=(17, 0.46 * (len(plan) + 8) + 2.0))
    colors = {
        "orf": get_color("primary"),
        "cassette": get_color("secondary"),
        "barcode": get_color("accent"),
        "readout": get_color("neutral"),
        "regulatory": get_color("light"),
        "igk": get_color("quaternary"),
        # CATEGORICAL_PALETTE[4]: a darker shade of the ORF red, because a native
        # signal peptide IS part of the native protein. get_color("accent_blue")
        # returns #3D6B60 -- the same value as "secondary" -- so it drew as cassette.
        "sp": CATEGORICAL_PALETTE[4],
    }

    # Lay out bottom-up with a gap between groups, so a group's title has clear air
    # instead of landing on the row above it.
    rows_rev = list(reversed(plan))
    span = max(len(row[3]) for row in plan)
    ypos, y, prev = {}, 0.0, None
    for i, row in enumerate(rows_rev):
        if prev is not None and row[5] != prev:
            y += 1.1  # the gap the group title sits in
        ypos[i] = y
        y += 1.0
        prev = row[5]
    top = y

    for grp in {r[5] for r in rows_rev}:
        members = [i for i, r2 in enumerate(rows_rev) if r2[5] == grp]
        lo, hi = ypos[min(members)], ypos[max(members)]
        ax.axhspan(lo - 0.5, hi + 0.5, color="#F3ECE5", zorder=0)
        ax.text(
            -600,
            hi + 0.72,
            grp.upper(),
            ha="left",
            va="center",
            fontsize=8,
            fontweight="bold",
            color="#6F5D5A",
        )

    for i, (plasmid, gene, where, insert, kind, _grp, why) in enumerate(rows_rev):
        y = ypos[i]
        for label, start, end, _btype in blocks(insert):
            if "cassette" in label:  # every cassette variant carries it in the label
                role = "cassette"
            elif label == "ORF, N-terminal" and kind == "igk":
                role = "igk"
            elif label == "ORF, N-terminal" and kind == "nterm_sp":
                role = "sp"
            else:
                role = ROLE.get(label, "regulatory")
            ax.barh(
                y,
                end - start,
                left=start,
                height=0.62,
                color=colors[role],
                edgecolor="white",
                linewidth=0.4,
            )
        ax.text(
            -60,
            y,
            f"{plasmid}  {gene}",
            ha="right",
            va="center",
            fontsize=8,
            fontfamily="monospace",
        )
        ax.text(
            span + 220,
            y,
            where,
            ha="left",
            va="center",
            fontsize=7.5,
            alpha=0.95 if kind == "igk" else 0.7,
            fontweight="bold" if kind == "igk" else "normal",
        )
        ax.text(span + 1500, y, why, ha="left", va="center", fontsize=7.5, color="#6F5D5A")

    ax.set_yticks([])
    ax.set_ylim(-0.9, top + 0.4)
    # index the insert explicitly: `*_x, i` took the LAST tuple element, which became
    # `kind` when the arm type was added, so the axis maxed at ~688 nt and
    # bbox_inches="tight" then stretched the canvas out to reach the row labels.
    ax.set_xlim(-620, span + 4300)
    ax.set_xlabel("position in the synthesised insert (nt)")
    # ticks only where constructs actually are; the right half is the text column
    ax.set_xticks([x for x in range(0, span + 1, 1000)])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    handles = [
        plt_patch(colors["orf"], "ORF (designed L1 CDS, domesticated)"),
        plt_patch(colors["cassette"], "tag cassette (51 aa)"),
        plt_patch(colors["igk"], "Igκ leader + Asp"),
        plt_patch(colors["sp"], "native signal peptide"),
        plt_patch(colors["barcode"], "barcode"),
        plt_patch(colors["readout"], "amplicon handles + 10x CS1"),
        plt_patch(colors["regulatory"], "Kozak / stop / pad"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    save_figure(fig, "s4e9_order_batch_map", output_dir=str(args.out_dir))
    print(f"{len(plan)} constructs drawn -> {args.out_dir}/s4e9_order_batch_map.png")
    return 0


def plt_patch(color: str, label: str):
    from matplotlib.patches import Patch

    return Patch(facecolor=color, label=label)


if __name__ == "__main__":
    raise SystemExit(main())
