"""The two figures for the brief's "Tag after the signal peptide" section.

    uv run --locked python signal_peptide_figures.py

1. ``figures/signal_peptide_rule_tree.png``: how each CDS in the N-after-SP decision table
   reaches its outcome. Every count is measured here from ``nterm/nterm_decisions.tsv``: each
   row is assigned to one rule by ``rung_numbers.NTERM_RULES`` (the same classifier that
   writes ``brief_numbers.json``), the counts are checked against ``brief_numbers.json``
   (``rungs.n_after_sp``), and every tree leaf is checked to sum back to the table. Where a
   CDS with no N-terminal construct went next (an internal site, or held) is read from
   ``l1_tag_paths.tsv``.
2. ``figures/signal_peptide_example.png``: three example CDSs over their first 60 residues.
   Top track: SignalP 6.0's cleavage-site probability for the untagged native protein and for
   the tagged protein at each candidate site (``nterm/nterm_members.tsv``,
   ``nterm/nterm_positions.tsv``; the Igκ forms from ``nterm/signalp_run_verdicts.tsv``).
   Bottom track: the three disorder tracks the interval reads (metapredict, NetSurfP-3.0 and
   AlphaFold 1 - pLDDT, from the frozen predictor snapshot that ``nterm/nterm_placement.py``
   reads, decoded by its ``evidence``), the cleavage site, the structure onset, the admissible
   interval and the chosen site. The interval is recomputed with
   ``nterm/admissible_interval.py`` and must equal the committed ``nterm_members.tsv`` row.

**A gap, stated rather than hidden.** SignalP 6.0 per-residue probabilities were never saved:
every pinned run used ``--format none`` (``nterm/signalp6_modal.py``), whose ``output.json``
leaves ``Plot_txt`` empty, and the re-prediction sweep stored per-protein calls only. So the
top track shows one cleavage-site probability per protein (native) or per candidate site
(tagged), not a per-residue curve.

Requested by the coordinator for the L1 design brief (2026-11-21); implemented by a Claude Code
session (branch ``claude/l1-signal-peptide-figures``).
"""

from __future__ import annotations

import csv
import json
import sys
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

from deliverome_analysis.plotting import get_color, save_figure, setup_plotting_style

HERE = Path(__file__).resolve().parent
S4E9 = HERE.parent
NTERM = HERE / "nterm"
FIGURES = HERE / "figures"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(NTERM))
sys.path.insert(0, str(S4E9))

import admissible_interval as ai  # noqa: E402
import rung_numbers as rn  # noqa: E402

# Construct colours: light fill, palette edge (deliverome_analysis.plotting families).
OUTCOME_STYLE = {
    "native": ("#CCE8E4", get_color("secondary"), "native signal peptide"),
    "igk": ("#FAECD4", get_color("warning"), "Igκ leader replaces it"),
    "none": ("#EDE7E3", get_color("neutral"), "no N-terminal construct"),
    "moved": ("#E4E0F8", get_color("info"), "moved to another placement"),
}


def _tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------------------------
# Figure 1: the rule tree
# ---------------------------------------------------------------------------------------------


def rule_counts() -> tuple[dict[str, list[dict]], list[dict]]:
    """Rows per NTERM_RULES class, each row in exactly one class; checked against the brief."""
    rows = _tsv(NTERM / "nterm_decisions.tsv")
    by_rule: dict[str, list[dict]] = {name: [] for name in rn.NTERM_RULES}
    for r in rows:
        hits = [name for name, test in rn.NTERM_RULES.items() if test(r)]
        assert len(hits) == 1, (r["design_id"], hits)
        by_rule[hits[0]].append(r)
    brief = json.loads((HERE / "brief_numbers.json").read_text())["rungs"]["n_after_sp"]
    assert {k: len(v) for k, v in by_rule.items()} == brief["by_rule"], "rule counts drifted"
    assert len(rows) == brief["decisions"], "nterm_decisions.tsv and brief_numbers.json disagree"
    assert dict(Counter(r["arm"] for r in rows)) == brief["arm"]
    assert brief["built"] == sum(r["arm"] in ("native_site", "igk_swap") for r in rows)
    return by_rule, rows


def leaf(label: str, outcome: str, rows: list[dict]) -> dict:
    return {"label": label, "outcome": outcome, "rows": rows, "children": []}


def node(label: str, children: list[dict], width: float = 2.3) -> dict:
    return {"label": label, "outcome": None, "children": children, "width": width}


def count(n: dict) -> int:
    return len(n["rows"]) if not n["children"] else sum(count(c) for c in n["children"])


def _v1_only_node(get) -> dict:
    """The branch for a signal peptide only DeepTMHMM v1 calls."""
    return node(
        "Only DeepTMHMM calls a signal peptide; SignalP 6 does not. Does a second source "
        "(UniProt, DeepTMHMM v2) confirm a cleaved one?",
        [
            leaf(
                "Yes: the Igκ leader replaces the called peptide",
                "igk",
                get("igk: only DeepTMHMM v1 calls the signal peptide; a second source confirms it"),
            ),
            leaf(
                "Yes, but it places the cleavage elsewhere: the leader is re-anchored there",
                "igk",
                get(
                    "igk: only DeepTMHMM v1 calls it; a second source places the cleavage elsewhere"
                )
                + get("igk: only DeepTMHMM v1 calls it; a signal anchor, its rung move deferred"),
            ),
            leaf(
                "No second source confirms it",
                "none",
                get(
                    "no arm: only DeepTMHMM v1 calls a signal peptide and no second source confirms it"
                ),
            ),
            leaf(
                "Displayed leader-free when N-terminally tagged in an earlier screen: "
                "tagged at residue 1, multipass",
                "moved",
                get("routed to the residue-1 multipass rung")
                + get(
                    "routed to another rung: a signal anchor or first TM helix, not a signal peptide"
                ),
            ),
        ],
        width=2.6,
    )


def build_tree(by_rule: dict[str, list[dict]]) -> dict:
    def get(name: str) -> list[dict]:
        return by_rule[name]

    refold_none = get("no arm: refold fails and no other SignalP-passing form")
    refold_none_native = [
        r for r in refold_none if "failed igk_swap" not in r["rejected_alternative"]
    ]
    refold_none_igk = [r for r in refold_none if "failed igk_swap" in r["rejected_alternative"]]
    assert len(refold_none_native) + len(refold_none_igk) == len(refold_none)

    v1_only = _v1_only_node(get)
    yes = node(
        "Yes: SignalP 6 keeps the native cleavage at one or more candidate sites",
        [
            leaf(
                "A cleavage-keeping site lies before the folded domain",
                "native",
                get(
                    "native: re-prediction-gated site before the predicted fold (rule 1, Becca 2026-11-11)"
                ),
            ),
            leaf(
                "Cleavage-keeping sites lie only inside the fold: the most disordered junction",
                "native",
                get("native: only preserved sites lie inside the predicted fold (rule 2)"),
            ),
            leaf(
                "The tag moves the cleavage site, but it stays upstream of the cassette",
                "native",
                get("native: cleavage shifts but stays upstream of the cassette (rule 3)"),
            ),
            leaf(
                "Ortholog CDS: same site rule, on the ortholog SignalP predictions",
                "native",
                get("native: ortholog SignalP release admitted (rule 5)"),
            ),
            leaf(
                "First SignalP prediction came from a later run; a site keeps cleavage",
                "native",
                get("native: round-2 SignalP on members with no prior call"),
            ),
            node(
                "The site lies in a UniProt pro-peptide, removed in the mature protein: "
                "moved to the mature N-terminus",
                [
                    leaf(
                        "SignalP keeps the cleavage at the mature N-terminus: built there",
                        "native",
                        get("native: moved past a UniProt N-terminal pro-peptide")
                        + get("igk: moved past a UniProt N-terminal pro-peptide"),
                    ),
                    leaf(
                        "No form at the mature N-terminus passes",
                        "none",
                        get(
                            "no arm: the mature N-terminus after a pro-peptide has no passing form"
                        ),
                    ),
                ],
            ),
            leaf(
                "But the native signal peptide is weak (cleavage-site probability < 0.5): "
                "Igκ leader",
                "igk",
                get(
                    "igk: weak predicted signal peptide (cleavage-site probability < 0.5; Becca 2026-11-14)"
                ),
            ),
            node(
                "The tagged model fails the refold check at the chosen site",
                [
                    leaf(
                        "The next cleavage-keeping site passes",
                        "native",
                        get("native: next SignalP-passing site after a failed refold"),
                    ),
                    leaf(
                        "Every site fails, the fold is intact: built and flagged",
                        "native",
                        get(
                            "native: refold fails at every SignalP-passing form; built and flagged (disposition rule)"
                        ),
                    ),
                    leaf(
                        "Every site fails: built at the first pick, flagged (flag, never hold)",
                        "native",
                        get(rn.FLAGGED),
                    ),
                    leaf("No other cleavage-keeping site", "none", refold_none_native),
                    leaf(
                        "Ortholog of a gene whose human CDS is not built",
                        "none",
                        get("no arm: ortholog rule, no human member of the gene is built"),
                    ),
                ],
            ),
        ],
        width=2.2,
    )
    no = node(
        "No: the native signal peptide is lost at every candidate site",
        [
            leaf(
                "An Igκ-leader form keeps the leader's cleavage: the leader replaces the native one",
                "igk",
                get("igk: native signal peptide fails at every site (rule 4)")
                + get("igk: round-2 SignalP, Igκ form keeps the leader's cleavage"),
            ),
            leaf(
                "The Igκ form fails the refold check; the next Igκ form passes",
                "igk",
                get("igk: next SignalP-passing form after a failed refold"),
            ),
            leaf("The Igκ form fails the refold check; no other passes", "none", refold_none_igk),
            leaf(
                "No tagged form, native or Igκ, keeps a cleavage site",
                "none",
                get("no arm: no tagged form keeps a SignalP cleavage site")
                + get("no arm: replacement arm's refold pending"),
            ),
            leaf(
                "No tagged form keeps a cleavage site: returned to an internal site",
                "moved",
                get("routed back to internal: no form keeps a cleavage site"),
            ),
        ],
        width=2.2,
    )
    question = node(
        "SignalP 6 calls the signal peptide. Is there a candidate site, between the predicted "
        "cleavage site and the start of the folded domain, where SignalP 6 on the TAGGED "
        "protein still predicts cleavage upstream of the cassette?",
        [yes, no],
        width=2.9,
    )
    return node("CDSs placed after the signal peptide", [question, v1_only], width=1.7)


def wrap(s: str, width_in: float, fontsize: float) -> str:
    chars = max(12, int(width_in * 72 / (fontsize * 0.58)))
    return "\n".join(textwrap.wrap(s, chars))


def _leaf_totals(tree: dict, leaves: list[dict], rows: list[dict]) -> Counter:
    """CDSs per outcome over the leaves, checked to cover the table and match its arms."""
    ids = [r["design_id"] + "|" + r["candidate_id"] for lf in leaves for r in lf["rows"]]
    assert len(ids) == len(set(ids)) == len(rows) == count(tree), "tree does not cover the table"
    totals = Counter()
    for lf in leaves:
        totals[lf["outcome"]] += len(lf["rows"])
    arms = Counter(r["arm"] for r in rows)
    assert totals["native"] == arms["native_site"] and totals["igk"] == arms["igk_swap"]
    assert totals["none"] == arms["no_nterm_arm"], "tree outcomes do not match the table's arms"
    return totals


def draw_tree(tree: dict, by_rule: dict, rows: list[dict]) -> plt.Figure:
    paths = {r["candidate_id"]: r for r in _tsv(S4E9 / "l1_tag_paths.tsv")}

    # every leaf's rows, disjoint, covering the table
    leaves: list[dict] = []

    def collect(n: dict, depth: int) -> None:
        n["depth"] = depth
        if not n["children"]:
            leaves.append(n)
        for c in n["children"]:
            collect(c, depth + 1)

    collect(tree, 0)
    totals = _leaf_totals(tree, leaves, rows)
    # where the CDSs with no construct went next
    for lf in leaves:
        if lf["outcome"] == "none":
            lf["internal"] = sum(
                paths[r["candidate_id"]]["placement"] == "internal" for r in lf["rows"]
            )
    n_internal = sum(lf.get("internal", 0) for lf in leaves)
    print(f"  no-construct CDSs placed on the internal rung instead: {n_internal}")

    fs, fs_leaf = 10.5, 10.5
    col_x = {0: 0.0, 1: 2.0, 2: 5.2, 3: 7.8, 4: 10.4}
    leaf_w, row_h, gap = 5.3, 0.62, 0.22
    # leaf y positions, top to bottom, with a gap between top-level branches
    y = 0.0
    prev_branch = None
    for lf in leaves:
        branch = id(_branch_of(tree, lf))
        if prev_branch is not None and branch != prev_branch:
            y -= gap * 2
        lf["y"] = y
        y -= row_h
        prev_branch = branch

    def place(n: dict) -> float:
        if n["children"]:
            ys = [place(c) for c in n["children"]]
            n["y"] = (ys[0] + ys[-1]) / 2
        return n["y"]

    place(tree)

    fig = plt.figure(figsize=(15.5, 0.62 * len(leaves) + 2.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.grid(False)

    def box(x, yc, w, h, fc, ec, lw=1.2):
        ax.add_patch(
            FancyBboxPatch(
                (x, yc - h / 2),
                w,
                h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                fc=fc,
                ec=ec,
                lw=lw,
                zorder=2,
            )
        )

    def draw(n: dict) -> None:
        d = n["depth"]
        if not n["children"]:
            x = col_x[4]
            fc, ec, _ = OUTCOME_STYLE[n["outcome"]]
            h = row_h - 0.1
            box(x, n["y"], leaf_w, h, fc, ec)
            ax.text(
                x + 0.1,
                n["y"],
                wrap(n["label"], leaf_w - 1.0, fs_leaf),
                ha="left",
                va="center",
                fontsize=fs_leaf,
                zorder=3,
                linespacing=1.05,
            )
            ax.text(
                x + leaf_w - 0.1,
                n["y"],
                f"{len(n['rows']):,}",
                ha="right",
                va="center",
                fontsize=13,
                fontweight="bold",
                color=ec,
                zorder=3,
            )
            if n["outcome"] == "none":
                k = n["internal"]
                note = (
                    f"{k} tagged at an internal site instead"
                    if k == len(n["rows"])
                    else "held"
                    if k == 0
                    else f"{k} tagged at an internal site instead, {len(n['rows']) - k} held"
                )
                ax.text(
                    x + leaf_w + 0.1,
                    n["y"],
                    note,
                    ha="left",
                    va="center",
                    fontsize=9.5,
                    color=get_color("neutral"),
                    style="italic",
                )
            return
        x, w = col_x[d], n["width"]
        text = wrap(n["label"], w - 0.15, fs)
        h = 0.2 + 0.19 * (text.count("\n") + 1) + 0.24
        box(x, n["y"], w, h, "white", get_color("dark"), lw=1.0)
        ax.text(
            x + w / 2,
            n["y"] + 0.13,
            text,
            ha="center",
            va="center",
            fontsize=fs,
            zorder=3,
            linespacing=1.05,
        )
        ax.text(
            x + w / 2,
            n["y"] - h / 2 + 0.17,
            f"{count(n):,} CDSs",
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            zorder=3,
        )
        for c in n["children"]:
            x0 = x + w
            x1 = col_x[4] if not c["children"] else col_x[c["depth"]]
            xm = x0 + 0.18
            ax.plot(
                [x0, xm, xm, x1],
                [n["y"], n["y"], c["y"], c["y"]],
                color=get_color("neutral"),
                lw=1.0,
                zorder=1,
                solid_joinstyle="miter",
            )
            draw(c)

    draw(tree)
    handles = [
        Patch(
            fc=OUTCOME_STYLE[k][0],
            ec=OUTCOME_STYLE[k][1],
            label=f"{OUTCOME_STYLE[k][2]} ({totals[k]:,})",
        )
        for k in ("native", "igk", "none", "moved")
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.995),
        fontsize=11,
        frameon=False,
        title="Construct",
        title_fontsize=11.5,
        alignment="left",
    )
    ax.set_xlim(-0.1, col_x[4] + leaf_w + 2.6)
    ax.set_ylim(y - 0.1, 0.45)
    return fig


def _branch_of(tree: dict, target: dict) -> dict:
    """The child of the tree's top-level question node (or the root) holding ``target``."""

    def contains(n: dict) -> bool:
        return n is target or any(contains(c) for c in n["children"])

    question = tree["children"][0]
    for top in (*question["children"], *tree["children"][1:]):
        if contains(top):
            return top
    return tree


# ---------------------------------------------------------------------------------------------
# Figure 2: three worked examples
# ---------------------------------------------------------------------------------------------

EXAMPLES = (
    ("l1-p0587", "native signal peptide kept; site before the fold"),
    ("l1-p0654", "native signal peptide kept; the only cleavage-keeping site is inside the fold"),
    (
        "l1-p2333",
        "native signal peptide lost at every site; Igκ leader replaces residues 1\u201330",
    ),
)
TRACKS = (
    ("metapredict", "metapredict", "primary"),
    ("netsurfp-3.0", "NetSurfP-3.0", "secondary"),
    ("alphafold-disorder", "AlphaFold 1 \u2212 pLDDT", "quaternary"),
)
VERDICT_STYLE = {
    "preserved": ("o", "secondary", "tagged: cleavage kept upstream"),
    "cleavage_moved": ("o", "accent", "tagged: cleavage moved, still upstream"),
    "cleavage_crosses_junction": ("D", "quaternary", "tagged: cleavage at or past the junction"),
    "signal_peptide_lost": ("X", "primary", "tagged: no signal peptide called (plotted at 0)"),
}
XMAX = 60


def example_data(processed) -> list[dict]:
    import nterm_placement as npl

    members = {r["design_id"]: r for r in _tsv(NTERM / "nterm_members.tsv")}
    positions = _tsv(NTERM / "nterm_positions.tsv")
    decisions = {r["design_id"]: r for r in _tsv(NTERM / "nterm_decisions.tsv")}
    igk = _tsv(NTERM / "signalp_run_verdicts.tsv")
    snap = npl.load_snapshot(processed)
    out = []
    for did, story in EXAMPLES:
        m, dec = members[did], decisions[did]
        rt = {
            "route": m["route"],
            "source_key": m["source_key"],
            "first_difference": int(m["first_difference"]) if m["first_difference"] else None,
        }
        _, tracks = npl.evidence(rt, snap)
        lo = int(m["lo"])
        iv = ai.tag_interval(lo, tracks)
        # the recomputed interval must be the committed one
        assert (iv.lo, iv.hi) == (int(m["interval_lo"]), int(m["interval_hi"])), did
        assert ";".join(f"{k}={v}" for k, v in sorted(iv.onsets.items())) == m["onsets"], did
        assert set(tracks) == set(m["tracks_answered"].split(";")), did
        cands = ai.candidate_positions(iv, lo)
        assert ",".join(map(str, cands)) == m["candidates"], did
        out.append(
            {
                "design_id": did,
                "gene": m["gene_symbol"],
                "story": story,
                "member": m,
                "decision": dec,
                "tracks": tracks,
                "interval": iv,
                "median_onset": int(float(m["median_onset"])),
                "positions": [p for p in positions if p["design_id"] == did],
                "igk": [r for r in igk if r["design_id"] == did and r["arm"] == "igk_tagged"],
                "igk_native": [
                    r for r in igk if r["design_id"] == did and r["arm"] == "igk_native"
                ],
            }
        )
    return out


def draw_examples(examples: list[dict]) -> plt.Figure:
    n = len(examples)
    fig = plt.figure(figsize=(15, 5.0 * n))
    outer = fig.add_gridspec(n, 1, hspace=0.42)
    legend_handles: dict[str, object] = {}
    for i, ex in enumerate(examples):
        inner = outer[i].subgridspec(2, 1, height_ratios=[1, 1.5], hspace=0.08)
        a_sp = fig.add_subplot(inner[0])
        a_dis = fig.add_subplot(inner[1], sharex=a_sp)
        m, dec, iv = ex["member"], ex["decision"], ex["interval"]
        lo, onset = iv.lo, ex["median_onset"]
        native_site = int(float(dec["native_after_residue"]))
        for ax in (a_sp, a_dis):
            split = min(max(onset - 0.5, lo + 0.5), iv.hi + 0.5)
            ax.axvspan(lo + 0.5, split, color=get_color("accent"), alpha=0.22, lw=0, zorder=0)
            ax.axvspan(
                split,
                iv.hi + 0.5,
                fc="none",
                ec=get_color("warning"),
                alpha=0.5,
                hatch="///",
                lw=0,
                zorder=0,
            )
            ax.axvline(lo + 0.5, color=get_color("dark"), lw=1.6, zorder=2)
            ax.axvline(onset - 0.5, color=get_color("dark"), lw=2.2, ls=":", zorder=3)
            ax.axvline(native_site + 0.5, color=get_color("danger"), lw=1.6, ls="--", zorder=2)
            ax.set_xlim(0.5, XMAX + 0.5)
        legend_handles["interval"] = Patch(
            color=get_color("accent"), alpha=0.35, label="admissible interval, before the fold"
        )
        legend_handles["interval_fold"] = Patch(
            fc="white",
            ec=get_color("warning"),
            hatch="///",
            label="admissible interval, inside the fold",
        )
        legend_handles["cleavage"] = Line2D(
            [],
            [],
            color=get_color("dark"),
            lw=1.6,
            label="cleavage site (later of SignalP 6, DeepTMHMM v1)",
        )
        legend_handles["onset"] = Line2D(
            [],
            [],
            color=get_color("dark"),
            lw=1.4,
            ls=":",
            label="structure onset (median of the three tracks' onsets)",
        )
        legend_handles["chosen"] = Line2D(
            [], [], color=get_color("danger"), lw=1.6, ls="--", label="chosen insertion site"
        )

        # --- SignalP track
        a_sp.axhline(0.5, color=get_color("neutral"), lw=1.0, ls=(0, (4, 3)), zorder=1)
        ncp = float(m["signalp_cleavage_probability"])
        a_sp.plot(
            int(float(m["signalp_cleavage"])) + 0.5,
            ncp,
            marker="s",
            ms=10,
            color=get_color("dark"),
            ls="none",
            zorder=4,
        )
        legend_handles["native"] = Line2D(
            [],
            [],
            marker="s",
            ms=9,
            color=get_color("dark"),
            ls="none",
            label="untagged native protein, at its cleavage site",
        )
        for p in ex["positions"]:
            mk, col, lab = VERDICT_STYLE[p["verdict"]]
            yv = (
                float(p["tagged_cleavage_probability"]) if p["tagged_cleavage_probability"] else 0.0
            )
            a_sp.plot(
                int(p["after_residue"]) + 0.5,
                yv,
                marker=mk,
                ms=10,
                color=get_color(col),
                ls="none",
                zorder=4,
                mec="white",
                mew=0.8,
            )
            legend_handles[p["verdict"]] = Line2D(
                [], [], marker=mk, ms=9, color=get_color(col), ls="none", label=lab
            )
        if ex["igk"]:
            head = len(rn.IGK_MATURE_LEADER) + 1
            rep = dec["igk_replaces_residues"]
            r_end = int(rep.split("-")[1])
            a_sp.axvspan(0.5, r_end + 0.5, color=get_color("warning"), alpha=0.06, lw=0, zorder=0)
            a_sp.text(
                1.2,
                0.08,
                f"residues 1\u2013{r_end} replaced by the Igκ leader + Asp ({head} aa)",
                fontsize=11,
                color=get_color("warning"),
                va="bottom",
            )
            for r in ex["igk"]:
                a_sp.plot(
                    int(float(r["native_after_residue"])) + 0.5,
                    float(r["cleavage_probability"]),
                    marker="^",
                    ms=11,
                    color=get_color("warning"),
                    ls="none",
                    zorder=5,
                    mec="white",
                    mew=0.8,
                )
            legend_handles["igk"] = Line2D(
                [],
                [],
                marker="^",
                ms=9,
                color=get_color("warning"),
                ls="none",
                label="Igκ-leader form, tagged: leader cleavage kept",
            )
        a_sp.set_ylim(-0.08, 1.08)
        a_sp.set_yticks([0, 0.5, 1])
        a_sp.set_ylabel("SignalP 6\ncleavage-site p", fontsize=12)
        a_sp.tick_params(labelbottom=False)
        fold = "inside" if dec["inside_predicted_fold"] == "1" else "before"
        cs = int(float(m["signalp_cleavage"]))
        dt = m["deeptmhmm_v1_cleavage"]
        detail = (
            f"SignalP 6 cleaves after {cs} (p = {ncp:.2f}), DeepTMHMM v1 after {dt}; "
            f"tag after native residue {native_site}, {fold} the fold"
        )
        if ex["igk_native"]:
            r0 = ex["igk_native"][0]
            detail += (
                f"; untagged Igκ form cleaves after leader residue "
                f"{int(float(r0['cleavage_site']))} (p = {float(r0['cleavage_probability']):.2f})"
            )
        a_sp.text(
            0.0,
            1.22,
            f"{chr(97 + i)}   {ex['gene']}",  # panel label only; the story goes in the caption
            transform=a_sp.transAxes,
            fontsize=13.5,
            fontweight="bold",
            va="bottom",
        )
        a_sp.text(
            0.0,
            1.05,
            detail,
            transform=a_sp.transAxes,
            fontsize=11,
            va="bottom",
            color=get_color("neutral"),
        )

        # --- disorder tracks
        for key, lab, col in TRACKS:
            t = ex["tracks"].get(key)
            if not t:
                continue
            xs = list(range(1, min(XMAX, len(t)) + 1))
            a_dis.plot(xs, t[: len(xs)], color=get_color(col), lw=1.8, zorder=3)
            o = iv.onsets.get(key)
            if o is not None and o <= XMAX:
                a_dis.plot(
                    o,
                    1.04,
                    marker="v",
                    ms=8,
                    color=get_color(col),
                    ls="none",
                    zorder=4,
                    clip_on=False,
                )
            legend_handles[key] = Line2D([], [], color=get_color(col), lw=2, label=lab)
        legend_handles["onset_marker"] = Line2D(
            [],
            [],
            marker="v",
            ms=8,
            ls="none",
            color=get_color("neutral"),
            label="a track's structure onset (track colour)",
        )
        a_dis.axhline(0.5, color=get_color("neutral"), lw=1.0, ls=(0, (4, 3)), zorder=1)
        a_dis.set_ylim(0, 1.0)
        a_dis.set_yticks([0, 0.5, 1])
        a_dis.set_ylabel("disorder", fontsize=12)
        if i == n - 1:
            a_dis.set_xlabel(
                "residue of the native protein (sites are drawn between residues n and n+1)"
            )
        else:
            a_dis.tick_params(labelbottom=True)
    order = [
        "native",
        "preserved",
        "cleavage_crosses_junction",
        "signal_peptide_lost",
        "igk",
        "metapredict",
        "netsurfp-3.0",
        "alphafold-disorder",
        "onset_marker",
        "cleavage",
        "onset",
        "interval",
        "interval_fold",
        "chosen",
    ]
    hs = [legend_handles[k] for k in order if k in legend_handles]
    hs.append(
        Line2D(
            [],
            [],
            color=get_color("neutral"),
            lw=1.0,
            ls=(0, (4, 3)),
            label="0.5: weak-signal-peptide bar (top); order threshold (bottom)",
        )
    )
    fig.legend(
        handles=hs,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.065),
        ncol=3,
        fontsize=11,
        frameon=False,
    )
    fig.subplots_adjust(bottom=0.13, top=0.97)
    return fig


def _leaf_rows(n: dict) -> list[dict]:
    if not n["children"]:
        return [n]
    return [lf for c in n["children"] for lf in _leaf_rows(c)]


def draw_simple(tree: dict) -> plt.Figure:
    """The rule tree in one panel (Becca, 2026-11-22: "can we make an additional simpler
    version?"): one bar per first question the CDS meets, split by the construct it ends with.
    Counts come from the same tree as the full figure, so the two cannot disagree."""
    question, v1_only = tree["children"]
    yes, no = question["children"]
    groups = [
        ("SignalP 6 keeps the native\ncleavage at a candidate site", yes),
        ("SignalP 6: the native peptide is\nlost at every candidate site", no),
        ("Only DeepTMHMM calls a\nsignal peptide", v1_only),
    ]
    order = ["native", "igk", "none", "moved"]
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    for gi, (_label, n) in enumerate(groups):
        counts = Counter()
        for lf in _leaf_rows(n):
            counts[lf["outcome"]] += len(lf["rows"])
        left = 0
        for o in order:
            v = counts[o]
            if not v:
                continue
            fill, edge, _ = OUTCOME_STYLE[o]
            ax.barh(gi, v, left=left, color=fill, edgecolor=edge, linewidth=1, height=0.62)
            if v >= 25:
                ax.text(left + v / 2, gi, f"{v:,}", ha="center", va="center", fontsize=9)
            left += v
        ax.text(left + 12, gi, f"{sum(counts.values()):,}", va="center", fontsize=9)
    ax.set_yticks(range(len(groups)), [g[0] for g in groups])
    ax.invert_yaxis()
    ax.set_xlabel("CDSs placed after the signal peptide")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(
        handles=[
            Patch(
                facecolor=OUTCOME_STYLE[o][0],
                edgecolor=OUTCOME_STYLE[o][1],
                label=OUTCOME_STYLE[o][2],
            )
            for o in order
        ],
        loc="lower right",
        frameon=False,
        fontsize=8.5,
    )
    fig.tight_layout()
    return fig


def main(argv: list[str] | None = None) -> int:
    """Draw both figures, or one with ``--only tree`` / ``--only examples`` / ``--only simple``.

    Each figure is stamped with the commit it was drawn from; drawing one at a time lets the
    other stay committed, so neither stamp reads as drawn from a dirty tree.
    """
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=("tree", "simple", "examples"))
    only = ap.parse_args(argv).only
    if only in (None, "tree"):
        setup_plotting_style()
        by_rule, rows = rule_counts()
        tree = build_tree(by_rule)
        fig = draw_tree(tree, by_rule, rows)
        save_figure(fig, "signal_peptide_rule_tree", output_dir=FIGURES)
        plt.close(fig)
    if only in (None, "simple"):
        setup_plotting_style()
        by_rule, _ = rule_counts()
        fig = draw_simple(build_tree(by_rule))
        save_figure(fig, "signal_peptide_rule_simple", output_dir=FIGURES)
        plt.close(fig)
    if only in (None, "examples"):
        from deliverome_analysis.s3 import processed

        setup_plotting_style()
        fig = draw_examples(example_data(processed))
        save_figure(fig, "signal_peptide_example", output_dir=FIGURES)
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
