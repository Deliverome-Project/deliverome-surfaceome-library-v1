# Where the art comes from

Vendored so the figure regenerates from this folder alone. These are copies,
not the record. If one changes upstream, re-copy it rather than editing the
copy here.

| File | Source | Used for |
|---|---|---|
| `platform-internalization_20260522.svg` | `20260522_Deliverome-platform_Functional Genomics Internalization Platform.svg` | **the figure.** 2011.58 × 970.12, live text (Manrope ExtraLight + Avenir Book), top-level groups `cells` / `receptors` / `barcodes` / `DNA*` / `antibodies` / `text`. The script places this whole file and drops only its `text` group |
| `fluorophore-hexagons-orange_050.svg` | Drive `BrandAssets/svg/molecules/`, id `1eNNBExNfI-x72dDmjAbMavEzqmQlSuMd` | the SpyCatcher swatch in the key strip |
| `fluorophore-hexagons-green_051.svg` | Drive `BrandAssets/svg/fluorescence/`, id `1bbluyuVfMWS9vV_JCzUOhmTkeg2LlXAB` | the anti-ALFA nanobody swatch in the key strip |

A fourth file exists upstream — the same platform figure with its type
converted to outlines, `…_scaled_outlines_Functional Genomics Internalization
Platform.svg`, 1088 × 662.53. It is not vendored here because nothing reads it:
outlined type cannot be re-set, which is the whole point of this script. It is
worth opening if you need to check the build against the original lettering.

## The platform master is not in Drive

Searching Drive for `Deliverome-platform` and for `Internalization Platform`
returns nothing: as of 2026-10-22 the master exists only in a local Downloads
folder and in copies like this one. The companion files
`…Functional Genomics Effector Platform copy.svg` and
`…Functional Genomics Trafficking Platform.svg` are in the same position.

That is the failure mode [CLAUDE.md](../../../../CLAUDE.md) names for construct
maps: a record that lives on one machine is not a record. These are the source
art for the platform figures, so they want a durable, stated home — Drive
beside `BrandAssets/`, or `external://` if they are to be cited from an
experiment. Until that happens this vendored copy is the only one under version
control, and nothing should cite the Downloads path.

## The glyphs are in Drive, drawn the same way

The master is built from the same vocabulary as the BrandAssets set: grey cell
discs (`cells/cell-outline-circle_*`), coloured receptors
(`receptors/receptor-y-purple_197`, `receptor-y-crimson_057`,
`receptor-rods-teal_030`, `receptor-squiggle-blue_070`), barcode helices
(`nucleic-acids/dna-helix-*`), the hexagon fluorophores, and
`icons/crispr-scissors_GO-1`. Drive also holds the readout column of this
figure already traced, as
`diagrams/figure-internalization-barcode-screen_107.svg`. None are needed here,
because the master already contains that art; they are where to go to draw a
*new* panel in the same style.

## One gotcha when re-using the master

It embeds its glows as raster `<image>` elements with `xlink:href`. Any document
including its body must declare
`xmlns:xlink="http://www.w3.org/1999/xlink"` on the root `<svg>`, or the whole
sheet fails to render — not just the glows.
