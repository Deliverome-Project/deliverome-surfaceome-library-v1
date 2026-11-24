"""Where to put SLC9A3's cassette: the built junction is the worst position tested.

SLC9A3 was the one arm of the eight whose signal-peptide call the cassette measurably
degraded without breaking (``canary_rp_1792371023``): p(SP) 0.867 -> 0.666 and the
cleavage-site probability 0.967 -> 0.650, at junction 26. It still predicts SP, so the
verdict reads ``preserved`` -- the verdict carries no confidence, so nothing in it says
this is the least-margin arm in the batch.

Junction 26 is ``sp_end + 1``, chosen by default rather than for SLC9A3. The protein has
room: DeepTMHMM v1 puts the signal peptide at 1-25, the extracellular N-terminal segment
at **26-56**, and TM1 at 57-80, so there are 31 residues to place a tag in, with no
N-glycosylation sequon and no cysteine anywhere in them. We used the first one.

This sweep re-predicts the same protein with the cassette at 26, 34, 40 and 46
(``canary_rp_1792374944``, SignalP 6.0+h; results archived at
``processed://sprint4/s4e9-orf-design-fg-panel/reprediction/canary_rp_1792374944/``).
Moving the junction off the cleavage site does not merely recover the call -- at 34 it
is **better than the untagged protein**, p(SP) 0.9998 against the native's 0.867.

The mechanism is the 70-residue SignalP window. At junction 26 the window holds ONE
native mature residue and then 44 residues of cassette, so almost the entire context
SignalP uses to place the cleavage site is synthetic. At 34 it holds nine, at 46 twenty-one.
That also explains the batch's pattern without needing a dose-response: the arms sitting
at sp_end+1 are the exposed ones, and the two with real offsets (FOLR1 +5, PTPRC +16)
lost 0.001 and 0.002.

**34 rather than 46**, on two grounds the probabilities do not show: it leaves 22 residues
between the cassette and TM1 where 46 leaves 10, and it sits in the Gly/Glu-rich stretch
26-45 rather than in the aromatic run 46-55 (TFEWAHV) that is more likely to be
structured. 46 scores marginally higher and is the riskier placement.

    uv run --locked python plot_slc9a3_junction_sweep.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

from deliverome_analysis.plotting import CATEGORICAL_PALETTE, save_figure, setup_plotting_style
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
KEY = "sprint4/s4e9-orf-design-fg-panel/reprediction/canary_rp_1792374944/results.json"
SP_END, TM1, BUILT = 25, 57, 26


def main() -> int:
    recs = json.loads(Path(processed.download(KEY)).read_text())
    native = next(r for r in recs if r["role"] == "native")
    tagged = sorted((r for r in recs if r["role"] == "tagged"), key=lambda r: r["after_residue"])

    setup_plotting_style()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    xs = [r["after_residue"] for r in tagged]
    site = [r["cleavage_probability"] for r in tagged]
    psp = [r["sp_probability"] for r in tagged]

    ax.axhline(
        native["cleavage_probability"],
        color=CATEGORICAL_PALETTE[3],
        ls="--",
        lw=1.4,
        label=f"untagged protein, cleavage-site p = {native['cleavage_probability']:.3f}",
    )
    ax.plot(
        xs,
        site,
        "o-",
        color=CATEGORICAL_PALETTE[0],
        ms=9,
        lw=2,
        label="cleavage-site probability, tagged",
    )
    ax.plot(
        xs,
        psp,
        "s-",
        color=CATEGORICAL_PALETTE[1],
        ms=7,
        lw=1.6,
        alpha=0.85,
        label="p(signal peptide), tagged",
    )
    ax.scatter(
        [BUILT],
        [next(r["cleavage_probability"] for r in tagged if r["after_residue"] == BUILT)],
        s=320,
        facecolors="none",
        edgecolors=CATEGORICAL_PALETTE[2],
        lw=2.4,
        zorder=5,
        label=f"junction {BUILT} — as built (p0025)",
    )

    for r in tagged:
        ax.annotate(
            f"{r['cleavage_probability']:.3f}",
            (r["after_residue"], r["cleavage_probability"]),
            textcoords="offset points",
            xytext=(0, 11),
            ha="center",
            fontsize=9,
        )

    ax.set_xlabel("cassette inserted after residue")
    ax.set_ylabel("SignalP 6 probability")
    ax.set_ylim(0.55, 1.06)
    ax.set_xlim(22, 50)
    ax.set_xticks(xs)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    save_figure(fig, HERE / "figures" / "s4e9_slc9a3_junction_sweep")

    print(f"native cleavage-site p = {native['cleavage_probability']:.4f}")
    for r in tagged:
        flag = "  <- as built" if r["after_residue"] == BUILT else ""
        print(
            f"  after {r['after_residue']:>3}: site {r['cleavage_probability']:.4f}  "
            f"p(SP) {r['sp_probability']:.4f}  "
            f"{TM1 - 1 - r['after_residue']:>2} residues to TM1{flag}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
