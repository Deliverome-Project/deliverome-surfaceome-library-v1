"""Plotting configuration for deliverome-analysis figures.

Uses seaborn styling with a Deliverome-aligned palette for publication-quality
figures. Callers must invoke :func:`setup_plotting_style` explicitly —
importing this module has no side effects.

Copied (and lightly adapted) from accessible-surfaceome's
``audit/_plotting_config.py`` so figures across the two repos share a visual
identity.
"""

from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# Bundled brand fonts (committed under ./assets/fonts/). Registered with
# matplotlib's font manager on first call to setup_plotting_style so the
# rcParams below actually resolve to Manrope instead of falling back to
# DejaVu Sans.
#
# Path math: this file is src/deliverome_analysis/plotting.py, so the repo
# root is parents[2].
_REPO_ROOT = Path(__file__).resolve().parents[2]
_BUNDLED_FONTS_DIR = _REPO_ROOT / "assets" / "fonts"
_fonts_registered = False


def register_bundled_fonts() -> None:
    """Add ``./assets/fonts/*.ttf`` to matplotlib's font manager (idempotent)."""

    global _fonts_registered
    if _fonts_registered:
        return
    if not _BUNDLED_FONTS_DIR.is_dir():
        _fonts_registered = True
        return
    for ttf in sorted(_BUNDLED_FONTS_DIR.glob("*.ttf")):
        try:
            fm.fontManager.addfont(str(ttf))
        except Exception:
            continue
    _fonts_registered = True


# Core Deliverome palette tokens (from Deliverome Design System)
COLORS = {
    # Brand anchors
    "primary": "#BC3C4C",  # maroon-light / accent
    "secondary": "#3D6B60",  # teal-mid
    "accent": "#F4AA28",  # amber-bright
    "quaternary": "#8878C8",  # lavender-bright
    # Semantic / utility
    "success": "#2E7A55",
    "warning": "#8C4210",  # amber-dark
    "danger": "#6E1428",  # maroon-dark
    "info": "#5848A8",  # lavender-mid
    "neutral": "#6F5D5A",  # muted
    "light": "#F3ECE5",  # bg-warm
    "dark": "#1F1718",  # ink
    "line": "#E6DAD4",
    # Legacy aliases for downstream code compatibility
    "accent_blue": "#3D6B60",
    "accent_gold": "#F4AA28",
    "accent_green": "#2E7A55",
}

# Categorical plotting palette (recommended figure anchors)
CATEGORICAL_PALETTE = [
    "#BC3C4C",  # maroon-light
    "#3D6B60",  # teal-mid
    "#F4AA28",  # amber-bright
    "#8878C8",  # lavender-bright
    "#6E1428",  # maroon-dark
    "#7AAB9F",  # teal-light
]

# Sequential palettes by hue family (deep -> light)
SEQUENTIAL_PALETTES = {
    "maroon": ["#3E0A18", "#6E1428", "#922038", "#BC3C4C", "#F0A098", "#FDE8E6"],
    "teal": ["#152E28", "#244840", "#3D6B60", "#4D8A80", "#7AAB9F", "#CCE8E4"],
    "amber": ["#5A2608", "#8C4210", "#C07830", "#F4AA28", "#F4C070", "#FAECD4"],
    "lavender": ["#1E1450", "#3A2888", "#5848A8", "#8878C8", "#A090D4", "#E4E0F8"],
}

# Diverging palette (lavender <-> neutral <-> maroon)
DIVERGING_PALETTE = sns.blend_palette(
    ["#8878C8", "#FFFFFF", "#BC3C4C"],
    n_colors=11,
)

# Default bar styling configuration
ROUND_BARS_DEFAULT = False


def setup_plotting_style(
    style: str = "default", context: str = "notebook", font_scale: float = 2.0
) -> None:
    """Set up modern plotting style for all figures.

    Parameters
    ----------
    style : str
        Seaborn style ('default', 'white', 'whitegrid', 'dark', 'darkgrid', 'ticks')
    context : str
        Seaborn context ('paper', 'notebook', 'talk', 'poster')
    font_scale : float
        Scale factor for fonts (default: 2.0 for larger, more readable text)
    """

    register_bundled_fonts()

    if style == "default":
        sns.set_style("whitegrid")
    else:
        sns.set_style(style)

    sns.set_context(context, font_scale=font_scale)
    sns.set(font_scale=font_scale)
    sns.set_palette(CATEGORICAL_PALETTE)

    plt.rcParams.update(
        {
            # Figure
            "figure.figsize": (5, 4),
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
            "figure.facecolor": "none",
            "savefig.facecolor": "none",
            # Font
            "font.family": "sans-serif",
            "font.sans-serif": ["Manrope", "Outfit", "DejaVu Sans", "Liberation Sans", "Arial"],
            "font.serif": ["Playfair Display", "Georgia", "Times New Roman", "serif"],
            "font.size": 14,
            "font.weight": "medium",
            # Axes
            "axes.labelsize": 16,
            # Titles suppressed by project convention (figures get captions instead).
            # Existing ax.set_title(...) calls become visual no-ops. To re-enable for
            # a one-off plot, set axes.titlesize back to e.g. 16 via plt.rcParams.
            "axes.titlesize": 0,
            "axes.titlepad": 0,
            "axes.titleweight": "semibold",
            "axes.labelweight": "medium",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.edgecolor": COLORS["line"],
            "axes.labelcolor": COLORS["dark"],
            "axes.facecolor": "none",
            "text.color": COLORS["dark"],
            # Grid
            "grid.alpha": 0.35,
            "grid.linestyle": "-",
            "grid.linewidth": 0.7,
            "grid.color": COLORS["line"],
            # Ticks
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "xtick.color": COLORS["dark"],
            "ytick.color": COLORS["dark"],
            # Legend
            "legend.frameon": True,
            "legend.framealpha": 0.9,
            "legend.fancybox": True,
            "legend.fontsize": 13,
            "legend.title_fontsize": 14,
            # Lines
            "lines.linewidth": 2,
            "lines.markersize": 8,
            # Patches (bars, etc.)
            "patch.edgecolor": "none",
            "patch.linewidth": 0.0,
        }
    )

    sns.despine(top=True, right=True)
    _disable_titles_globally()


# Originals stashed once at module import time so this is idempotent —
# repeated setup_plotting_style() calls don't re-wrap an already-wrapped
# function.
_ORIGINAL_AXES_SET_TITLE = Axes.set_title
_ORIGINAL_FIGURE_SUPTITLE = Figure.suptitle


def _disable_titles_globally() -> None:
    """Replace Axes.set_title and Figure.suptitle with no-ops.

    Callers that want a title back can call ``set_title`` /  ``suptitle``
    via the ``_ORIGINAL_*`` references on this module, or set the
    rcParam ``axes.titlesize`` to a positive value and re-import.
    """

    def _no_title(self: object, *_args: Any, **_kwargs: Any) -> None:
        return None

    # Deliberate monkey-patch of the matplotlib classes (house style: titles go in captions).
    Axes.set_title = _no_title  # type: ignore[method-assign,assignment]  # ty: ignore[invalid-assignment]
    Figure.suptitle = _no_title  # type: ignore[method-assign,assignment]  # ty: ignore[invalid-assignment]


def save_figure(
    fig: Figure,
    filename: str,
    output_dir: str | Path = "figures",
    formats: Iterable[str] = ("png",),
    gist_url: str | None = None,
    dpi: float = 300,
    transparent: bool = False,
) -> None:
    """Save figure in multiple formats.

    Defaults to one white-background PNG so a single run creates one reviewable
    artifact. Request PDF explicitly for camera-ready vector output with
    ``formats=('png', 'pdf')``. SVG remains available when explicitly required,
    but is not generated by default and should not be uploaded to Notion.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save
    filename : str
        Base filename (without extension)
    output_dir : str or Path
        Output directory
    formats : iterable[str]
        Output formats. Defaults to PNG only.
    gist_url : str or None
        Optional reproduction-link URL. Embedded as PNG ``Source`` tEXt chunk
        and PDF ``Subject`` field. Useful for figures that ship with a
        figure-reproduction gist.

    Every PNG/PDF also carries the git commit, uncommitted-changes flag, script and lockfile
    hashes that produced it (see :mod:`deliverome_analysis.provenance`).
    dpi : int
        Resolution in dots per inch. Defaults to 300.
    transparent : bool
        PNG only: save with a transparent background instead of white. For slide overlays
        stacked on a white base figure (an on-click build); leave False for anything else.
    """
    import json

    from deliverome_analysis.provenance import PNG_KEY, provenance, warn_if_dirty

    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    # which commit and locked package versions made this figure (read back with
    # `python -m deliverome_analysis.provenance <file>`)
    stamp = provenance()
    stamp_json = json.dumps(stamp, separators=(",", ":"))

    for fmt in formats:
        filepath = output_path / f"{filename}.{fmt}"
        # bbox_inches="tight" crops to the artists (so long axis labels / outside
        # annotations are never clipped), and pad_inches keeps them off the very
        # edge so text doesn't sit flush against the image border.
        save_kwargs: dict[str, Any] = {
            "format": fmt,
            "dpi": dpi,
            "bbox_inches": "tight",
            "pad_inches": 0.15,
        }
        if fmt == "png":
            if transparent:
                save_kwargs.update(facecolor="none", transparent=True)
            else:
                save_kwargs.update(facecolor="white", transparent=False)
        if fmt == "png":
            save_kwargs["metadata"] = {PNG_KEY: stamp_json}
            if gist_url:
                save_kwargs["metadata"]["Source"] = gist_url
        elif fmt == "pdf":
            save_kwargs["metadata"] = {"Keywords": stamp_json}
            if gist_url:
                save_kwargs["metadata"]["Subject"] = gist_url
        fig.savefig(filepath, **save_kwargs)
        print(f"  Saved: {filepath}")
    warn_if_dirty(stamp, filename)


def create_figure(
    nrows: int = 1, ncols: int = 1, figsize: tuple[float, float] | None = None, **kwargs: Any
) -> tuple[Figure, Any]:
    """Create a figure with modern styling.

    If ``figsize`` is None, auto-sizes to (5*ncols, 4*nrows).
    """
    if figsize is None:
        width = 5 * ncols
        height = 4 * nrows
        figsize = (width, height)
    return plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, **kwargs)


def add_significance_bar(
    ax: Axes,
    x1: float,
    x2: float,
    y: float,
    height: float = 0.05,
    text: str = "***",
    **kwargs: Any,
) -> None:
    """Add a significance bar between two points."""
    bar_kwargs: dict[str, Any] = {"color": COLORS["dark"], "linewidth": 1.5}
    bar_kwargs.update(kwargs)
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], **bar_kwargs)
    ax.text(
        (x1 + x2) / 2,
        y + height,
        text,
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color=COLORS["dark"],
    )


def format_percentage_axis(ax: Axes, axis: str = "y") -> None:
    """Format an axis to show percentages."""
    from matplotlib.ticker import PercentFormatter

    if axis == "y":
        ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    else:
        ax.xaxis.set_major_formatter(PercentFormatter(1.0))


def get_color(name: str) -> str:
    """Get a color from the palette by name."""
    return COLORS.get(name, COLORS["primary"])


def get_categorical_colors(n: int | None = None) -> list[Any]:
    """Get categorical color palette."""
    if n is None:
        return CATEGORICAL_PALETTE
    return (
        CATEGORICAL_PALETTE[:n] if n <= len(CATEGORICAL_PALETTE) else sns.color_palette("husl", n)
    )


def get_diverging_palette(n: int = 11) -> list[Any]:
    """Get diverging color palette."""
    return cast(list[Any], sns.blend_palette(["#8878C8", "#FFFFFF", "#BC3C4C"], n_colors=n))


def get_sequential_palette(family: str = "maroon", n: int = 6) -> list[Any]:
    """Get a sequential palette by family.

    Parameters
    ----------
    family : str
        One of ``maroon``, ``teal``, ``amber``, ``lavender``.
    n : int
        Number of colors to return.
    """
    if family not in SEQUENTIAL_PALETTES:
        family = "maroon"
    palette = SEQUENTIAL_PALETTES[family]
    if n <= len(palette):
        return palette[:n]
    return cast(list[Any], sns.blend_palette(palette, n_colors=n))
