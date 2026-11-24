# Pooled internalization assay

The EndoNBV2 assay — how a pool of barcoded reporter cells is read for how much
of each surface receptor gets internalized, normalized for how much of it the
cell displays.

It is the method schematic for the proposal's internalization arm: the figure a
reviewer needs before any internalization number means anything, because the
whole readout is a ratio of two probes and that only reads if both probes are
drawn.

Rendered by [fig_internalization_assay.py](fig_internalization_assay.py):
`python3 library/l1/figures/fig_internalization_assay.py`

No `uv run` — the script is stdlib only, by design (see *House rules* below).

## Source art

The drawing is the Deliverome platform figure itself,
`20260522_Deliverome-platform_Functional Genomics Internalization Platform.svg`,
vendored at [`assets/`](assets/MANIFEST.md). **The script does not redraw
anything.** It places that file whole, drops its `text` group, and re-sets the
lettering.

The method it draws is the pooled extension of EndoNB,
[Lenaerts et al. 2025](https://doi.org/10.1101/2025.06.08.658482), bioRxiv
(Leonardo Almeida-Souza lab, a collaboration). The published assay reads one
endogenously ALFA-tagged protein at a time with a nanobody–3C–SNAPtag probe;
the pooled version adds the second tag and probe that carry the expression
reference, and the genomic barcode that makes it a screen.

**The platform master is not in Drive.** As of 2026-10-22 it exists only in a
local Downloads folder and in the copy here, as do the Effector and Trafficking
platform files. [`assets/MANIFEST.md`](assets/MANIFEST.md) records this; it
wants fixing before the figure is cited from anywhere durable.

## What the figure is for, and the one thing that was wrong with it

The master is a slide. At 6.5 in wide its own lettering lands between **4.4 and
8.4 pt** — below what any funder will accept, and not an export artefact: the
master's `font-size` really is 19–36 px on a 2011-unit canvas.

It cannot be fixed inside the artwork, because at 8 pt the two reagent
definitions are wider than the gaps the drawing leaves for them. So:

| | Master | Here |
|---|---|---|
| type at 6.5 in | 4.4 / 6.0 / 7.2 / 8.4 pt | 9.0 / 8.5 / 8.0 pt, and 7.1 pt for the two inline tag names |
| weight | Manrope ExtraLight (200) | Manrope 500–700 |
| families | Manrope + Avenir Book | Manrope throughout |
| the two reagent definitions | three- and four-line blocks inside the drawing | one line each, in a key strip across the top |
| sheet | 2011.58 × 970.12, with margin | 1980 × 1038, cropped to the ink plus the key strip |

3.4 in tall at 6.5 in wide, against 3.1 in for the master. That quarter-inch is
what legible type costs.

## Label positions

Labels are placed through `art(x, y)`, which maps the master's own coordinates
onto the sheet, so each sits where the original put it. Four had to move,
because bigger type does not fit the pocket the master left for it. All four
were placed by measuring the art around them, not by eye:

| Label | Master anchor | Here | Why |
|---|---|---|---|
| `SpyTag` | 469, 232 | 471, 252 | at 30 units it ran across the probe's own receptor and the next cell's green dot. Moved into the clear block under the orange bar (sheet x 415–525, y 300–345), bounded by the probe stem at 409 and the receptor below at 350 |
| `ALFATag` | 558, 235 | 616, 252 | on the same baseline as `SpyTag`, with a 35-unit gap; at the master's anchor the two names touched |
| `barcode` | 355, 503 | 415, 503 | centred on the nucleus and helix it names (sheet x 359) instead of left-aligned where the master's text box began |
| `Sort and sequence barcodes` | 1837, 723 | 1841, 706 | the third line ran through the right-hand sort arrow. Raised into the gap between the readout cells and the fork: text y 765–880, arrow at 890 |

`37 °C` keeps the master's position and is centred on its arrow.

If the master is ever redrawn these anchors go stale **silently** — the labels
still render, just in the wrong place. Re-render and look at it: the pockets
they sit in are only tens of units wide.

## Draft caption

> **Figure X. Pooled internalization assay.** Each cell in the library carries
> one surface receptor tagged in its own locus with both a SpyTag and an ALFA
> tag, and a genomic barcode identifying which receptor it is. SpyCatcher binds
> the SpyTag covalently and is not cleavable, so it reports how much receptor
> that cell displays; an anti-ALFA nanobody binds the ALFA tag through a
> cleavable tether. At 37 °C both probes traffic with the receptor, so a
> fraction is internalized into endosomes. A reducing agent or protease then
> removes the anti-ALFA signal from whatever remains at the surface and cannot
> reach what was internalized. Anti-ALFA signal divided by SpyCatcher signal is
> therefore internalization normalized for surface expression; cells are sorted
> on that ratio and the barcodes in each gate sequenced, scoring every receptor
> in the pool in one experiment. Schematic: no quantity in the figure is
> measured. Assay after Lenaerts et al., bioRxiv 2025
> (doi:10.1101/2025.06.08.658482), extended to a pooled, barcoded format.

## Caveats

| Caveat | Why it matters |
|---|---|
| **Nothing in the figure is measured.** | The three readout cells, the ladder and every glyph are illustrative. If any of it starts reading as a result, flatten it |
| **The citation is in the caption, not on the sheet.** | The master carried a Lenaerts credit line as a footer; it was removed so the figure sits in a document. The house rule that a figure citing published work carries that citation on the figure is aimed at replotted data — this is a schematic of our own assay — but the caption must carry it, and does above |
| **The pooled format is ours, not the published method.** | Lenaerts et al. read one protein at a time. The second tag, the second probe and the barcode are the extension. Text next to the figure should not imply the pooled version is published |
| **"Reducing agent or protease" is deliberately both.** | Our runs (S5E2) use TCEP on a reducible linker; the published method uses an HRV 3C site. The figure covers both because the design has not been fixed to one |
| **Tag placement is schematic.** | ALFA distal and Spy membrane-proximal, so the two probes do not compete for the same face. If the real construct differs, the drawing is wrong in a way no reviewer could catch |
| **Both dyes must be spectrally separable in the sorter.** | Assumed, not drawn. The ratio readout depends on it |

## House rules this figure does not follow, and why

This folder's rules assume a data figure. Two do not apply, and one is met
differently:

- **`setup_plotting_style()` / `save_figure()`** are not used, and there is no
  provenance stamp. This is vector art, not a plot; there is no matplotlib
  figure to pass. Nothing here is measured, so a stamp would certify only which
  commit drew a cartoon.
- **"Show individual observations"** has nothing to apply to.
- **"Regenerated by the script beside it"** holds in full: no hand-edited PNG,
  and the script is stdlib only (`re`, `math`, `subprocess`, `pathlib`) so it
  runs without the locked environment. It shells out to `qlmanage` and `sips`,
  so the PNG step is macOS-only; the SVG is written anywhere.

`figures/` carries both the PNG the convention asks for and the SVG, because
this is a vector figure and the SVG is what gets placed in the document.
