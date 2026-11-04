"""Read a Sony MA900 `.expdat` experiment export.

Why this exists. The MA900's FCS files carry no gain information at all -- not `$PnG`, not a
vendor keyword, nothing -- so for a long time the acquisition settings behind a sort were simply
unknown to us. They are not lost: they live in the `.expdat` the instrument can export, which is
a ZIP of SQL Server bcp tables whose settings columns are gzipped, base64- or hex-encoded XML.
This module reads that, so a sort run can record what it was actually acquired with.

It also recovers the worksheet gate hierarchy, which matters for a second reason: repo policy is
that gates must be recorded in some reviewable form, and for an MA900 run this is the only
machine-readable record there is.

How the gate geometry is recovered, and how far it has been checked.

Gates are stored as vertices in *plot* coordinates -- positions on the worksheet's own 350-unit
plot -- so turning them into data coordinates needs the plot extent, the axis scale, and the axis
range. All three are in the file, but which entry applies to which region is not documented, so
each step below was checked rather than assumed.

* **Plot extent** comes from `SetGateInfo`, not from `SetGridWidth`/`SetGridHeight`. The latter is
  snap-to-grid spacing -- 182 and 176 against a true extent of 350 -- and using it puts vertices
  at 21,000,000 on a channel whose maximum is 1,000,000.

* **Axis range** is stored per gate id, and each entry is complete: every channel, not just the
  two its plot shows. It describes the plot on which that gate is the *displayed* population, so
  a region converts against its PARENT's entry, the parent being whose plot it was drawn on.
  Confirmed by measuring the worksheet screenshot the `.expdat` carries: the tick spacing of the
  sort plot gives FL10-A 319..42,472 and FL1-A 100..114,690, against the stored parent entry's
  325..42,200 and 100..112,671. Agreement to 2%, which is the precision of reading pixels.

* **Y is not inverted.** Checked against events: the singlet gate resolves to 81.7% upright and
  0.0% inverted.

* **Ellipses** are stored as three control points with no documented meaning. The convention --
  centre = P1, radii = P2, rotation = P3.x degrees -- was recovered by testing candidates against
  the one property a hierarchy must satisfy for correlated scatter parameters: a child is drawn on
  a plot of its parent's population, so its events must sit inside the parent. It puts 97.9% of
  the singlet gate inside its scatter parent, against 37.3%, 10.3% and 2.7% for the alternatives,
  and rotation sign separates cleanly (+37.12 deg 97.9%, 0 deg 88.2%, -37.12 deg 72.9%). Reading
  the three points as a triangle instead gives 7.1% of events where a scatter gate should give far
  more. Vertices are generated in plot space and converted point by point, so an ellipse on a log
  axis stays a correct curve rather than being bent by converting centre and radii separately.

Drawing the result over the events reproduces the worksheet screenshot: the scatter ellipse on the
cell cloud, the singlet band along the diagonal, the viability rectangle at the low end of BV421,
and the three sort bands across the FL10/FL1 diagonal. See `scratch/becca/ma900_gate_overlay.py`.

Two things that are **not** checked, and should not be read as if they were.

Comparing gate *percentages* against that screenshot does not work. Its Gates and Statistics table
reports 100,000 events; no saved acquisition on the S3E9 plate holds 100,000 -- the sort tube holds
79,845, in the FCS and in the `.expdat`'s own per-channel event files alike -- and 100,000 is
exactly the acquisition stop condition. Those statistics come from a setup stream rather than from
the recorded sort, so the spread they show against the fluorescence gates reflects two different
samples. The sort itself was recorded in full. Closing this needs per-gate statistics computed on
an acquisition that was saved.

And parent containment is reported by `resolve_gates` for review but deliberately does not gate
anything. It is only meaningful where parent and child share correlated parameters. A child drawn
on unrelated axes has no geometric obligation to its parent at all -- on S3E9 gate B sits 97.9%
inside A while G sits 62.8% inside C -- because hierarchy there is applied by the analysis, not
implied by the regions.
"""

from __future__ import annotations

import base64
import gzip
import itertools
import math
import re
import subprocess
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Each plot's pixel extent, from its worksheet item. NOT SetGridWidth/Height, which is the
# snap-to-grid spacing (182/176 here) and produces vertices far outside the axis range.
_PLOT_EXTENT = re.compile(
    r"SetGateInfo, \$\$TargetWorkshettId\$\$, [0-9a-f-]+, [^,]*, \w+, [\d.]+, [\d.]+, ([\d.]+), ([\d.]+)"
)
# Sony names side scatter BSC; FCS files from other instruments and our own helpers say SSC.
CHANNEL_ALIASES = {"BSC-A": "SSC-A", "BSC-H": "SSC-H", "BSC-W": "SSC-W"}
# Points used to stand an ellipse up as a polygon. 64 keeps the area error under ~0.2%.
ELLIPSE_VERTICES = 64
# Axis configuration is stored per GATE id, and each entry is complete -- every channel, not
# just the two its plot shows. It describes the plot on which that gate is the *displayed*
# population, so a region is converted against its PARENT's entry: the parent is the population
# whose plot the region was drawn on. On S3E9 that is the difference between FL10-A running
# 100..1,000,000 and 129..82,264.
_AXIS = re.compile(
    r"SetAxis(Scale|MaximumValue|MinimumValue), \$\$TargetWorkshettId\$\$, ([0-9a-f-]+), "
    r"([\w-]+), ([^\n,]+)"
)
_REGION = re.compile(
    r"AddRegion, \$\$TargetWorkshettId\$\$, ([0-9a-f-]+), ([\w-]+), ([\w-]+), *, "
    r"(\w+), ([^,]+), ([0-9a-f-]+)"
)
_REGION_NAME = re.compile(r"SetRegionName, \$\$TargetWorkshettId\$\$, ([0-9a-f-]+), ([^\n,]+)")


def _inflate(value: bytes) -> str | None:
    """Settings columns are gzip under either base64 or hex; both appear in one file."""
    try:
        if value[:4] == b"H4sI":
            return gzip.decompress(base64.b64decode(value)).decode("utf-8", "replace")
        if value[:4] == b"1F8B":
            return gzip.decompress(bytes.fromhex(value.decode())).decode("utf-8", "replace")
    except Exception:
        return None
    return None


def _first(xml: str, tag: str) -> str | None:
    m = re.search(rf"<{tag}>([^<]*)</{tag}>", xml)
    return m.group(1) if m else None


def gating_ml_id(name: str) -> str:
    """A gate name that survives Gating-ML export.

    Sony names sort regions "L - 1" and "L - 2". Those go out as Gating-ML ids containing
    spaces, and the result does not parse back -- FlowKit rejects the whole document as "neither
    GatingML 2.0 compliant nor a FlowJo workspace", so the failure is total rather than local to
    the offending gate. The instrument's own name is kept on the Gate for display and for the
    join to sort destinations; only the exported id is rewritten.
    """
    return re.sub(r"\s+", "_", name.strip())


@dataclass
class Gate:
    name: str
    parent: str
    kind: str
    channels: tuple[str, ...]
    vertices_plot: list[tuple[float, float]]
    vertices_data: list[tuple[float, float]] | None = None
    #: True when the worksheet holds more than one axis range for a channel this gate uses, so

    @property
    def is_resolvable(self) -> bool:
        """False for shapes we cannot express without guessing at their geometry."""
        return self.kind in ("Polygon", "Rectangle", "Ellipse") and bool(self.vertices_data)

    def as_recipe_gate(self) -> dict[str, Any]:
        """agentflow recipe shape, for the kinds that map without guesswork."""
        if not self.is_resolvable:
            raise ValueError(
                f"gate {self.name!r} is a {self.kind}, which has no unambiguous recipe equivalent"
            )
        gate: dict[str, Any] = {
            "name": gating_ml_id(self.name),
            "parent": "root" if self.parent == "root" else gating_ml_id(self.parent),
            "channels": list(self.channels),
            "instrument_name": self.name,
        }
        if self.kind == "Rectangle":
            xs = [v[0] for v in self.vertices_data]
            ys = [v[1] for v in self.vertices_data]
            gate["kind"] = "rectangle"
            gate["bounds"] = [min(xs), max(xs), min(ys), max(ys)]
        else:
            gate["kind"] = "polygon"
            gate["vertices"] = [list(v) for v in self.vertices_data]
            if self.kind == "Ellipse":
                # A true ellipse has no recipe equivalent; this is a polygon through points on
                # it, which is lossy in a stated way rather than a different region.
                gate["approximated_from"] = "ellipse"
        return gate


@dataclass
class SortDestination:
    """One collection destination and what actually landed in it.

    A two-way sort writes several of these, and reporting only the first understates the run:
    on S3E9 the three destinations took 5,969, 15,794 and 2,258 cells, so the first row alone
    accounts for a quarter of what was sorted. They are kept separate rather than summed because
    which gate went where is the part a later analysis needs, and a total cannot be unmixed.
    """

    name: str
    gate: str
    row: int | None = None
    column: int | None = None
    mode: str = ""
    cell_size: str = ""
    cells_sorted: int = 0
    drops_sorted: int = 0
    cells_aborted: int = 0
    elapsed_ms: int = 0

    @property
    def abort_fraction(self) -> float | None:
        total = self.cells_sorted + self.cells_aborted
        return self.cells_aborted / total if total else None

    @property
    def sorted_per_second(self) -> float | None:
        return 1000 * self.cells_sorted / self.elapsed_ms if self.elapsed_ms else None


@dataclass
class Ma900Experiment:
    name: str
    sorter_serial: str | None
    #: when the sort started, from the instrument's own counters -- not the export timestamp,
    #: which moves every time the operator re-exports the same run
    sort_started: str | None = None
    detectors: dict[str, dict[str, Any]] = field(default_factory=dict)
    lasers: list[dict[str, str]] = field(default_factory=list)
    flow: dict[str, str] = field(default_factory=dict)
    chip: dict[str, str] = field(default_factory=dict)
    sort_destinations: list[SortDestination] = field(default_factory=list)
    compensation: dict[str, Any] = field(default_factory=dict)
    axes: dict[str, dict[str, str]] = field(default_factory=dict)
    #: every worksheet in the file, keyed by where it came from. A run carries more than one:
    #: S3E9 has five distinct worksheets across ten tubes, holding 6, 7 or 10 regions, and only
    #: the sort's own worksheet defines the L / L-1 / L-2 regions that fed its destinations.
    #: Reading one and applying it everywhere loses exactly the gates a sort is judged on.
    gate_sets: dict[str, list[Gate]] = field(default_factory=dict)
    gates_are_validated: bool = False
    #: child gate name -> fraction of its events that also fall inside its parent
    parent_containment: dict[str, float] = field(default_factory=dict)

    @property
    def gates(self) -> list[Gate]:
        """The worksheet to read by default: the sort's own if this run sorted, else the richest.

        A sort is judged on the gates that drove it, so those take precedence. Falling back to
        the worksheet with the most regions rather than the first is deliberate -- the first tube
        in the file is usually an unstained control, whose worksheet is a subset.
        """
        if "sort" in self.gate_sets:
            return self.gate_sets["sort"]
        return max(self.gate_sets.values(), key=len, default=[])

    @property
    def sort_gates(self) -> list[Gate]:
        """The regions the sort destinations drew from, empty when this run did not sort."""
        names = {d.gate for d in self.sort_destinations if d.gate}
        return [g for g in self.gate_sets.get("sort", []) if g.name in names]

    @property
    def unmatched_sort_gates(self) -> list[tuple[str, str]]:
        """(destination, gate) for destinations whose gate is not in the sort worksheet.

        The join is on the gate name that WellSortSettings records, so a name that does not
        appear in the worksheet leaves `sort_gates` quietly short rather than raising. That would
        mean the record cannot say what a destination collected, which is the one thing a sort
        record is for, so it is surfaced instead of being left to be noticed.
        """
        known = {g.name for g in self.gate_sets.get("sort", [])}
        return [(d.name, d.gate) for d in self.sort_destinations if d.gate and d.gate not in known]

    @property
    def is_sort(self) -> bool:
        """Whether this run sorted anything.

        This is the queryable fact, which is why there is no separate "tag this as a sort" step:
        a run either has destinations that took cells or it does not, and that is read from the
        instrument's own record rather than asserted by whoever filed it.
        """
        return any(d.cells_sorted for d in self.sort_destinations)

    @property
    def sort_totals(self) -> dict[str, Any]:
        """Summed across destinations. The per-destination rows stay the source of truth."""
        sorted_ = sum(d.cells_sorted for d in self.sort_destinations)
        aborted = sum(d.cells_aborted for d in self.sort_destinations)
        elapsed = max((d.elapsed_ms for d in self.sort_destinations), default=0)
        return {
            "destinations": len(self.sort_destinations),
            "cells_sorted": sorted_,
            "drops_sorted": sum(d.drops_sorted for d in self.sort_destinations),
            "cells_aborted": aborted,
            "elapsed_ms": elapsed,
            "abort_fraction": aborted / (sorted_ + aborted) if sorted_ + aborted else None,
            "sorted_per_second": 1000 * sorted_ / elapsed if elapsed else None,
        }

    @property
    def unresolvable_gates(self) -> list[str]:
        return [g.name for g in self.gates if not g.is_resolvable]

    def to_recipe(self) -> dict[str, Any]:
        """An agentflow recipe. Transforms follow the worksheet's own axis scales.

        Refuses when any gate cannot be expressed without guessing, rather than emitting a
        recipe that looks complete and gates on the wrong region. Dropping the gate is not an
        option either: its children name it as their parent, and a recipe whose parents are
        missing is invalid.
        """
        if self.unresolvable_gates:
            raise ValueError(
                "cannot build a recipe while these gates are unresolved: "
                + ", ".join(self.unresolvable_gates)
                + ". Redraw them as polygons in the instrument software."
            )
        transforms: dict[str, Any] = {}
        for gate in self.gates:
            for ch in gate.channels:
                scale = (self.axes.get(ch) or {}).get("Scale", "Linear")
                transforms[ch] = (
                    {"kind": "linear"}
                    if scale == "Linear"
                    else {
                        "kind": "asinh",
                        "cofactor": float(
                            (self.axes.get(ch) or {}).get("MinimumValue", 100) or 100
                        ),
                    }
                )
        # agentflow's compensation spec is {"mode": none|fcs|matrix}, not a bare label. "none" is
        # the honest value here rather than a default: the controls were defined and never
        # acquired, and the FCS spillover is the identity, so there is no compensation to apply.
        return {
            "version": 1,
            "compensation": {"mode": "fcs" if self.compensation.get("any_recorded") else "none"},
            "transforms": transforms,
            "gates": [g.as_recipe_gate() for g in self.gates],
        }


def _ellipse_points(control: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Expand Sony's three ellipse control points into points on the ellipse, in plot space.

    Convention recovered empirically (see the module docstring): centre, then the two radii as a
    point, then the rotation in degrees as the third point's x. Expanding here, in plot space,
    rather than after conversion is deliberate -- an ellipse is only an ellipse in the coordinates
    it was drawn in, so on a log axis converting the centre and radii separately would give the
    wrong curve.
    """
    (cx, cy), (rx, ry), (degrees, _) = control
    angle = math.radians(degrees)
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    points = []
    for i in range(ELLIPSE_VERTICES):
        t = 2 * math.pi * i / ELLIPSE_VERTICES
        dx, dy = rx * math.cos(t), ry * math.sin(t)
        points.append((cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a))
    return points


def _to_data(value: float, extent: float, axis: dict[str, str]) -> float:
    """Plot coordinate -> data coordinate, using the axis's own scale and range."""
    lo = float(axis.get("MinimumValue", 0) or 0)
    hi = float(axis.get("MaximumValue", 1) or 1)
    frac = value / extent if extent else 0.0
    if axis.get("Scale") == "Log":
        lo = lo or 1.0
        return lo * (hi / lo) ** frac
    return lo + frac * (hi - lo)


def read_experiment(path: str | Path) -> Ma900Experiment:
    """Parse an `.expdat`. Raises if it is not one."""
    path = Path(path)
    with zipfile.ZipFile(path) as zf:
        names = zf.namelist()
        if not any(n.startswith("Formats/") for n in names):
            raise ValueError(f"{path} does not look like an MA900 .expdat export")
        blobs = {n: zf.read(n) for n in names if "/" in n and not n.endswith("/")}

    def table(prefix: str) -> bytes:
        for n, data in blobs.items():
            if n.startswith(prefix):
                return data
        return b""

    experiment_row = table("Experiments/").split(b"\t")
    name = experiment_row[1].decode("latin-1") if len(experiment_row) > 1 else path.stem
    serial = None
    m = re.search(rb'"SorterSerialNumber":(\d+)', table("Experiments/"))
    if m:
        serial = m.group(1).decode()

    started = re.search(rb"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", table("SortResults/"))
    exp = Ma900Experiment(
        name=name,
        sorter_serial=serial,
        sort_started=started.group(0).decode() if started else None,
    )

    # --- acquisition settings, from the sample group's MeasurementSettings blob -------------
    for field_bytes in table("SampleGroups/").split(b"\t"):
        xml = _inflate(field_bytes)
        if not xml or "<Voltage>" not in xml:
            continue
        for m in re.finditer(r"<Voltage>([^<]*)</Voltage>", xml):
            ctx = xml[max(0, m.start() - 320) : m.end()]
            found = re.findall(r"<Name>([^<]*)</Name>", ctx)
            if found:
                exp.detectors[found[-1]] = {"voltage": m.group(1)}
        for m in re.finditer(r"<PDGain>([^<]*)</PDGain>", xml):
            ctx = xml[max(0, m.start() - 320) : m.end()]
            found = re.findall(r"<Name>([^<]*)</Name>", ctx)
            exp.detectors.setdefault(found[-1] if found else "FSC", {})["pd_gain"] = m.group(1)
        for m in re.finditer(r"<Name>(LD\d+)</Name><Power>([^<]*)</Power>", xml.replace("\n", "")):
            exp.lasers.append({"name": m.group(1), "power": m.group(2)})
        for tag in ("SheathPressure", "SamplePressure", "SamplePressureLevel"):
            value = _first(xml, tag)
            if value:
                exp.flow[tag] = value
        break

    # --- chip, from any tube's ChipInfo blob ------------------------------------------------
    for field_bytes in table("Tubes/").split(b"\t"):
        xml = _inflate(field_bytes)
        if xml and "OrificeSize" in xml:
            for tag in ("OrificeSize", "ChipId", "Generation", "AppType", "LotNumber"):
                value = _first(xml, tag)
                if value:
                    exp.chip[tag] = value
            break

    # --- sort destinations and what landed in each -------------------------------------------
    # Results and settings are separate tables joined on WellSortSetting_Id, which is the only
    # way to say *what* was sorted rather than just how many: the settings row carries the
    # destination's name and the gate feeding it.
    settings: dict[str, dict[str, Any]] = {}
    for line in table("WellSortSettings/").decode("latin-1").replace("\r\n", "\n").split("\n"):
        parts = line.split("\t")
        if len(parts) < 12 or not re.fullmatch(r"[0-9A-Fa-f-]{36}", parts[0].strip()):
            continue
        settings[parts[0].strip()] = {
            "name": parts[1],
            "row": int(parts[2]) if parts[2].isdigit() else None,
            "column": int(parts[3]) if parts[3].isdigit() else None,
            "gate": parts[5],
            "mode": parts[7],
            "cell_size": parts[8],
        }

    for line in table("WellSortResults/").decode("latin-1").replace("\r\n", "\n").split("\n"):
        parts = line.split("\t")
        if len(parts) < 5 or not re.fullmatch(r"[0-9A-Fa-f-]{36}", parts[0].strip()):
            continue
        try:
            elapsed, sorted_, drops, aborted = (int(parts[i]) for i in (1, 2, 3, 4))
        except ValueError:
            continue
        setting = next(
            (settings[p.strip()] for p in parts[::-1] if p.strip() in settings),
            {"name": "unknown", "gate": ""},
        )
        exp.sort_destinations.append(
            SortDestination(
                name=setting["name"],
                gate=setting["gate"],
                row=setting.get("row"),
                column=setting.get("column"),
                mode=setting.get("mode", ""),
                cell_size=setting.get("cell_size", ""),
                cells_sorted=sorted_,
                drops_sorted=drops,
                cells_aborted=aborted,
                elapsed_ms=elapsed,
            )
        )

    # --- compensation: controls are often defined but never acquired -------------------------
    cdf = table("CompensationDataFiles/").decode("latin-1")
    controls = []
    for m in re.finditer(r"([0-9A-Fa-f-]{36})\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t", cdf):
        controls.append(
            {
                # bcp writes a NUL for an empty column, which .strip() does not remove.
                "control": m.group(2).strip("\x00 \t"),
                "sensor": m.group(4),
                "recorded": bool(m.group(3).strip("\x00")),
            }
        )
    exp.compensation = {
        "controls": controls,
        "any_recorded": any(c["recorded"] for c in controls),
        "matrix": "identity",
    }

    # --- gates, from the worksheet script ----------------------------------------------------
    # --- gates: every worksheet in the file, not just the first ------------------------------
    # A run carries several. S3E9 has five distinct worksheets across ten tubes holding 6, 7 or
    # 10 regions, and only the sort's own defines the L / L-1 / L-2 regions its destinations drew
    # from. Reading one and applying it everywhere loses exactly the gates a sort is judged on.
    def _worksheets():
        for label, prefix in (("tube", "Tubes/"), ("sort", "SortResults/")):
            index = 0
            for chunk in table(prefix).split(b"\t"):
                xml = _inflate(chunk)
                if not xml or "<Script>" not in xml:
                    continue
                yield ("sort" if label == "sort" else f"tube[{index}]"), xml
                index += 1

    for key, xml in _worksheets():
        body = re.search(r"<Script>(.*?)</Script>", xml, re.S).group(1).replace("&#xD;", "")
        extents = {float(w) for w, _ in _PLOT_EXTENT.findall(body)} | {
            float(h) for _, h in _PLOT_EXTENT.findall(body)
        }
        if len(extents) != 1:
            continue
        extent = extents.pop()
        # Axis ranges are stored per GATE id, but several plots can display the same population
        # on different channel pairs with different zooms -- S3E9 shows C three times, on
        # FL1-A x FSC-A, FL10-A x FSC-A and FL10-A x FL1-A. One entry cannot describe three
        # plots, so when a channel carries more than one range in a worksheet, which one a given
        # region was drawn against is not determinable from the file. Collect them all and mark
        # the affected gates rather than silently picking the last one seen.
        plot_axes: dict[str, dict[str, dict[str, str]]] = {}
        for m in _AXIS.finditer(body):
            plot_axes.setdefault(m.group(2), {}).setdefault(m.group(3), {})[m.group(1)] = m.group(
                4
            ).strip()
            # Flattened per channel for reporting; conversion uses the per-plot entries above.
            exp.axes.setdefault(m.group(3), {})[m.group(1)] = m.group(4).strip()
        names = dict(_REGION_NAME.findall(body))
        by_id: dict[str, str] = {}
        collected: list[Gate] = []
        for m in _REGION.finditer(body):
            parent_id, x_ch, y_ch, kind, verts, gate_id = m.groups()
            points = []
            for pair in verts.split("|"):
                if ":" not in pair:
                    continue
                a, b = pair.split(":")[:2]
                try:
                    points.append((float(a), float(b)))
                except ValueError:
                    continue
            if kind == "Ellipse" and len(points) == 3:
                points = _ellipse_points(points)
            gate_name = names.get(gate_id, gate_id[:8]).strip()
            by_id[gate_id] = gate_name
            drawn_on = plot_axes.get(parent_id, {})
            x_axis = drawn_on.get(x_ch) or exp.axes.get(x_ch, {})
            y_axis = drawn_on.get(y_ch) or exp.axes.get(y_ch, {})
            collected.append(
                Gate(
                    name=gate_name,
                    parent=by_id.get(parent_id, "root"),
                    kind=kind,
                    channels=(x_ch, y_ch),
                    vertices_plot=points,
                    vertices_data=[
                        (_to_data(px, extent, x_axis), _to_data(py, extent, y_axis))
                        for px, py in points
                    ],
                )
            )
        exp.gate_sets[key] = collected

    return exp


def resolve_gates(experiment: Ma900Experiment, events: dict[str, Any]) -> dict[str, float]:
    """Resolve each gate against real events and record that the transform holds.

    `events` maps FCS channel name to a 1-D array. Sony's BSC is mapped to SSC on the way in.
    Sets `gates_are_validated` when every gate whose channels are present resolves to a
    non-degenerate fraction -- an inverted or mis-scaled transform collapses gates to ~0%, which
    is what makes this worth running rather than assuming.
    """
    import numpy as np
    from matplotlib.path import Path as MplPath

    fractions: dict[str, float] = {}
    masks: dict[str, Any] = {}
    resolved = 0
    for gate in experiment.gates:
        channels = [CHANNEL_ALIASES.get(c, c) for c in gate.channels]
        if not all(c in events for c in channels) or not gate.is_resolvable:
            continue
        x = np.asarray(events[channels[0]], dtype=float)
        y = np.asarray(events[channels[1]], dtype=float)
        if gate.kind == "Rectangle":
            xs = [v[0] for v in gate.vertices_data]
            ys = [v[1] for v in gate.vertices_data]
            inside = (x >= min(xs)) & (x <= max(xs)) & (y >= min(ys)) & (y <= max(ys))
        else:
            inside = MplPath(gate.vertices_data).contains_points(np.column_stack([x, y]))
        fractions[gate.name] = float(inside.mean())
        masks[gate.name] = inside
        resolved += 1

    # Containment of a child's region inside its parent's is reported, but it is only meaningful
    # when the two gates share correlated parameters. A child drawn on entirely different axes --
    # G on FL10 x FL1 under C on FL6 x FSC -- has no geometric obligation to sit inside its parent
    # at all: hierarchy there is applied by the analysis, not implied by the regions. On S3E9 that
    # shows plainly, B sits 97.9% inside A while G sits 62.8% inside C, and neither number says
    # anything about whether the transform is right. So this informs review; it does not gate.
    containment: dict[str, float] = {}
    for gate in experiment.gates:
        parent = masks.get(gate.parent)
        child = masks.get(gate.name)
        if parent is None or child is None or not child.any():
            continue
        containment[gate.name] = float((child & parent).sum() / child.sum())
    experiment.parent_containment = containment
    experiment.gates_are_validated = resolved > 0 and any(f > 0.01 for f in fractions.values())
    return fractions


# Heuristics for the run-quality verdict. These are review triggers chosen from one plate, not
# measured thresholds, which is why the verdict names what it saw rather than asserting a cause.
RATE_DROP_FLAG = 0.5  # last window below half the first counts as a material fall
MEDIAN_DRIFT_FLAG = 0.10  # 10% shift in a scatter median
# A median from few events is noisy on its own, so a fixed percentage flags sparse windows for
# nothing: S3E9's purity check has ~42 cells per window and tripped a 3% drift that is pure
# sampling. The standard error of a median is about 1.253 * sigma / sqrt(n), so the threshold
# also has to clear a multiple of that before a drift means anything.
MEDIAN_SE_FACTOR = 1.253
DRIFT_SIGMA_MULTIPLE = 3.0
CV_DRIFT_FLAG = 5.0  # percentage points of FSC CV
GATE_FRACTION_FLAG = 0.80  # root gate keeping under 80% of its starting share


def time_series_qc(
    experiment: Ma900Experiment, events: dict[str, Any], bins: int = 6
) -> dict[str, Any]:
    """Split a run by its TIME channel and report whether it held up.

    The point is the pairing, not either number alone. Event rate falling on its own says little:
    a sample settles, a tube runs low, an operator changes the pressure. What separates a clog
    from that is whether the cells still look the same on the way through -- a clog shears cells
    and destabilises the stream, so the scatter medians shift, the spread widens and the debris
    fraction climbs alongside the falling rate. On the S3E9 plate the rate falls 4.3-fold while
    every distributional metric stays flat, which is settling, not a blockage.

    One case this cannot separate: a partial clog that restricts flow without damaging cells looks
    the same. Aborts over time would settle it, and the .expdat records only the total. So the
    verdict flags what it saw for review; it does not diagnose.
    """
    import numpy as np

    if "TIME" not in events:
        raise ValueError("no TIME channel; the run cannot be split over time")

    time = np.asarray(events["TIME"], dtype=float)
    scatter = {
        name: np.asarray(events[name], dtype=float) for name in ("FSC-A", "SSC-A") if name in events
    }
    gate_masks: dict[str, Any] = {}
    for gate in experiment.gates:
        channels = [CHANNEL_ALIASES.get(c, c) for c in gate.channels]
        if not gate.is_resolvable or not all(c in events for c in channels):
            continue
        from matplotlib.path import Path as MplPath

        x = np.asarray(events[channels[0]], dtype=float)
        y = np.asarray(events[channels[1]], dtype=float)
        if gate.kind == "Rectangle":
            xs = [v[0] for v in gate.vertices_data]
            ys = [v[1] for v in gate.vertices_data]
            gate_masks[gate.name] = (
                (x >= min(xs)) & (x <= max(xs)) & (y >= min(ys)) & (y <= max(ys))
            )
        else:
            gate_masks[gate.name] = MplPath(gate.vertices_data).contains_points(
                np.column_stack([x, y])
            )

    # Debris and aggregates swamp the scatter spread, so measure the cells: the root-most
    # resolvable gate, which is the scatter gate in every MA900 worksheet seen so far.
    population = np.ones(len(time), dtype=bool)
    for gate in experiment.gates:
        if gate.parent == "root" and gate.name in gate_masks:
            population = gate_masks[gate.name]
            break

    edges = np.linspace(time.min(), time.max(), bins + 1)
    windows = []
    for low, high in itertools.pairwise(edges):
        inside = (time >= low) & (time < high)
        if not inside.any():
            continue
        window: dict[str, Any] = {
            "from_s": float(low),
            "to_s": float(high),
            "events": int(inside.sum()),
            "cells": int((inside & population).sum()),
            "events_per_second": float(inside.sum() / (high - low)) if high > low else None,
        }
        for name, values in scatter.items():
            selected = values[inside & population]
            if not selected.size:
                continue
            window[f"{name}_median"] = float(np.median(selected))
            mean = selected.mean()
            window[f"{name}_cv_percent"] = float(100 * selected.std() / mean) if mean else None
        for name, mask in gate_masks.items():
            window[f"gate_{name}_fraction"] = float(mask[inside].mean())
        windows.append(window)

    verdict, notes = "stable", []
    if len(windows) >= 2:
        first, last = windows[0], windows[-1]
        rate_ratio = (
            last["events_per_second"] / first["events_per_second"]
            if first["events_per_second"]
            else 1.0
        )
        drifted = False

        def noise_floor(window: dict[str, Any], name: str) -> float:
            """Relative sampling noise on this window's median, from its own n and spread."""
            cells = window.get("cells") or window.get("events") or 0
            cv = window.get(f"{name}_cv_percent")
            if not cells or not cv:
                return float("inf")
            return MEDIAN_SE_FACTOR * (cv / 100) / math.sqrt(cells)

        for name in scatter:
            a, b = first.get(f"{name}_median"), last.get(f"{name}_median")
            if not (a and b):
                continue
            observed = abs(b - a) / a
            combined = math.hypot(noise_floor(first, name), noise_floor(last, name))
            if observed > MEDIAN_DRIFT_FLAG and observed > DRIFT_SIGMA_MULTIPLE * combined:
                drifted = True
                notes.append(f"{name} median moved {100 * (b - a) / a:+.0f}%")
        a, b = first.get("FSC-A_cv_percent"), last.get("FSC-A_cv_percent")
        if a and b and abs(b - a) > CV_DRIFT_FLAG:
            drifted = True
            notes.append(f"FSC-A CV moved {b - a:+.1f} points")
        root = next(
            (g.name for g in experiment.gates if g.parent == "root" and g.name in gate_masks), None
        )
        if root is not None:
            a = first.get(f"gate_{root}_fraction")
            b = last.get(f"gate_{root}_fraction")
            if a and (b is None or b / a < GATE_FRACTION_FLAG):
                drifted = True
                share = 0.0 if b is None else b / a
                notes.append(f"gate {root} kept {share:.2f} of its starting share of events")
        if rate_ratio < RATE_DROP_FLAG:
            notes.insert(0, f"rate fell to {rate_ratio:.2f} of its starting value")
            verdict = "rate fell with degrading scatter" if drifted else "rate fell, scatter stable"
        elif drifted:
            verdict = "scatter drifted at a steady rate"
    return {"windows": windows, "verdict": verdict, "notes": notes}


def _sort_run_id(experiment: Ma900Experiment, started: str | None) -> str:
    """A run id tied to the sort itself, not to when someone happened to export the file.

    The `.expdat` filename carries its export timestamp, which moves every time the operator
    re-exports the same run. The sort's own start time does not.
    """
    if started:
        return "ma900-" + re.sub(r"[^0-9T]", "", started.replace(" ", "T"))[:15] + "Z"
    return "ma900-" + (experiment.sorter_serial or "unknown")


def write_run_record(
    experiment: Ma900Experiment,
    *,
    analysis_path: str | Path,
    expdat_key: str,
    expdat_local: str | Path,
    input_keys: list[str],
    started: str | None = None,  # defaults to the sort's own start time
    record_dir: str | Path | None = None,
) -> Path:
    """Write the `run-records/<id>.yaml` for a sort the instrument ran.

    Same schema as the agentflow wrapper's records, because a reviewer should not have to learn
    two shapes. What differs is where provenance comes from: an agentflow run proves itself with
    its own `run.json`, while an instrument sort has no such artefact, so the immutable evidence
    is the `.expdat` on S3 plus its checksum. Everything derived from it -- settings, destinations,
    gates -- is recorded with the limitations that apply, because a record that omits them reads
    as stronger than it is.
    """
    import yaml

    from .flow_records import RUN_RECORD_SCHEMA_VERSION, digest, relative, repository

    analysis_path = Path(analysis_path)
    root = repository(analysis_path)
    run_id = _sort_run_id(experiment, started or experiment.sort_started)
    directory = Path(record_dir) if record_dir else analysis_path / "run-records"
    directory.mkdir(parents=True, exist_ok=True)

    module = Path(__file__).resolve()
    source_files = [
        {
            "role": "instrument_export",
            "s3_key": f"raw://{expdat_key}",
            "sha256": digest(expdat_local),
        },
        {"role": "reader", "path": relative(module, root), "sha256": digest(module)},
    ]
    if root:
        for name in ("uv.lock", "pyproject.toml"):
            if (root / name).is_file():
                source_files.append({"role": name, "path": name, "sha256": digest(root / name)})

    git_commit = None
    if root:
        result = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        )
        git_commit = result.stdout.strip() or None
    dirty = None
    if root:
        result = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            capture_output=True,
            text=True,
            check=False,
        )
        dirty = bool(result.stdout.strip())

    totals = experiment.sort_totals
    document = {
        "schema_version": RUN_RECORD_SCHEMA_VERSION,
        "experiment": relative(analysis_path, root) if root else str(analysis_path),
        "run_id": run_id,
        "sort_started": experiment.sort_started,
        "produced_by": "deliverome_analysis.ma900",
        "instrument": {
            "model": "Sony MA900",
            "serial": experiment.sorter_serial,
            "experiment_name": experiment.name,
            "detectors": experiment.detectors,
            "lasers": experiment.lasers,
            "flow": experiment.flow,
            "chip": experiment.chip,
        },
        "compensation": {
            "controls_defined": [
                c["control"] for c in experiment.compensation["controls"] if c["control"].strip()
            ],
            "any_recorded": experiment.compensation["any_recorded"],
            "matrix": "identity",
        },
        "sort": {
            "destinations": [
                {
                    "name": d.name,
                    "gate": d.gate,
                    "row": d.row,
                    "column": d.column,
                    "mode": d.mode,
                    "cell_size": d.cell_size,
                    "cells_sorted": d.cells_sorted,
                    "drops_sorted": d.drops_sorted,
                    "cells_aborted": d.cells_aborted,
                    "abort_fraction": round(d.abort_fraction, 4) if d.abort_fraction else None,
                }
                for d in experiment.sort_destinations
            ],
            "totals": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in totals.items()},
        },
        "gates": {
            "worksheets": {k: len(v) for k, v in experiment.gate_sets.items()},
            "used": "sort" if "sort" in experiment.gate_sets else "largest",
            "names": [g.name for g in experiment.gates],
            "unresolvable": experiment.unresolvable_gates,
            "destinations_without_a_matching_gate": experiment.unmatched_sort_gates,
        },
        "source_files": source_files,
        "inputs": [f"raw://{k}" for k in input_keys],
        "input_count": len(input_keys),
        "analysis_source": {"git_commit": git_commit, "working_tree_dirty": dirty},
        "warnings": _record_warnings(experiment),
    }
    path = directory / f"{run_id}.yaml"
    path.write_text(yaml.safe_dump(document, sort_keys=False, allow_unicode=True))
    return path


def _record_warnings(experiment: Ma900Experiment) -> list[str]:
    """State the limits in the record itself; a silent record reads as a stronger claim."""
    warnings = [
        "Gate geometry is recovered from undocumented fields. The coordinate transform is "
        "confirmed against the worksheet rendering the .expdat carries, to about 2%, but gate "
        "percentages have not been checked against instrument-computed statistics: the only "
        "statistics the file stores are on a setup stream that was never saved.",
        "The ellipse control-point convention is inferred, not documented. See the module "
        "docstring for the evidence.",
    ]
    if not experiment.compensation["any_recorded"]:
        warnings.append(
            "Compensation controls were defined but none were acquired, and the spillover matrix "
            "is the identity. Any spectral overlap in this run is uncorrected."
        )
    if experiment.unmatched_sort_gates:
        warnings.append(
            "Some sort destinations name a gate absent from the worksheet, so what they "
            "collected cannot be stated: "
            + ", ".join(f"{d} -> {g}" for d, g in experiment.unmatched_sort_gates)
        )
    return warnings
