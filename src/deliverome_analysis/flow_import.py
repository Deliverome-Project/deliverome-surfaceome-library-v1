"""Publish flow analyses that Agentflow did not produce, in the Agentflow run layout.

The flow index (:mod:`deliverome_analysis.qcdb`) and the event-level SQL in
:mod:`deliverome_analysis.flowdb` read published Agentflow runs. Two other kinds of analysis
exist and belong in the same place:

* **FlowJo workspaces.** :func:`import_flowjo` applies a ``.wsp``'s gates to its FCS files with
  FlowKit (vendored in Agentflow), checks every population count against the count FlowJo saved
  in the workspace, and publishes only when they agree within :data:`TOLERANCE`.
* **Custom re-gating** (threshold or mask gates computed in a script). Build the per-sample
  gate masks yourself and pass them to :func:`write_run`; the experiment's committed script is
  the recipe.

Both write a run folder that :func:`deliverome_analysis.flow.publish_run` accepts and the index
reads like any other: ``events.parquet`` (``raw:``/``signal:``/``gate:`` columns),
``summary.csv`` (the contract in :data:`qcdb.FLOW_SUMMARY_COLUMNS`), ``run.json``,
``reproducibility.yaml`` and ``samples.csv``. ``run.json`` carries ``producer`` (``flowjo`` or
``script``) so imported runs are never mistaken for Agentflow gating.

    uv run --group flow python -m deliverome_analysis.flow_import flowjo \\
        --experiment sprint2/s2e17-lifting-method-surface-proteome \\
        --wsp experiments/sprint2/s2e17/11-Sep-2026_Yan_live-dead.wsp
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as md
import json
import platform
import re
import shutil
import tempfile
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from deliverome_analysis import qcdb

# FlowKit reproduces FlowJo gates closely but not to the event: polygon gates are drawn in each
# axis's transformed space and FlowKit's transforms differ slightly from FlowJo's. Measured on our
# workspaces (S2E10, S2E17, S3E9, S4E1, S4E6, S4E7), percent of parent agrees within ~1 point
# (median < 0.1), while raw counts can be offset by up to ~2% (S4E7's scatter gate is uniformly
# ~1.4% low and the offset carries down the hierarchy). A population passes when BOTH hold:
#   count:   |FlowKit - FlowJo| <= max(TOLERANCE x FlowJo count, TOLERANCE_EVENTS)
#   percent: |% of parent, FlowKit - FlowJo| <= TOLERANCE_POINTS, for parents of >= 200 events
# Wrong files, wrong gates or a misread compensation fail both by a wide margin.
TOLERANCE = 0.025
TOLERANCE_EVENTS = 10
TOLERANCE_POINTS = 1.5
_MIN_PARENT_FOR_PERCENT = 200
STATS = ("median", "mean", "min", "max", "p05", "p25", "p75", "p95", "sd")
_SKIP_STATS = {"Time"}


def sha256(path: str | Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


@dataclass
class GatedSample:
    """One sample's events and gate membership, ready to write.

    ``raw``: events as acquired (columns = detector names). ``signal``: the same events in the
    run's signal space (compensated when the analysis compensated; else ``raw``). ``gates``:
    gate name -> boolean mask over all events, in hierarchy order (parents before children).
    """

    sample_id: str
    raw: pd.DataFrame
    signal: pd.DataFrame
    gates: dict[str, np.ndarray]
    sha256: str
    source_uri: str
    metadata: dict = field(default_factory=dict)
    signal_space: str = "raw"


def _stats(frame: pd.DataFrame) -> dict:
    out = {}
    for det in frame.columns:
        if det in _SKIP_STATS:
            continue
        v = frame[det].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        out[f"finite_signal_count:{det}"] = int(v.size)
        if v.size:
            q = np.percentile(v, [5, 25, 75, 95])
            vals = [
                np.median(v),
                v.mean(),
                v.min(),
                v.max(),
                *q,
                v.std(ddof=1) if v.size > 1 else 0,
            ]
        else:
            vals = [np.nan] * len(STATS)
        out.update({f"{s}_signal:{det}": float(x) for s, x in zip(STATS, vals, strict=True)})
    return out


def summarize(samples: Sequence[GatedSample], parents: Mapping[str, str]) -> pd.DataFrame:
    """The contract summary: one row per sample x gate (plus ``root``), counts and statistics."""
    rows = []
    for s in samples:
        total = len(s.raw)
        masks = {"root": np.ones(total, dtype=bool), **s.gates}
        for gate, mask in masks.items():
            parent = parents.get(gate, "root")
            count, parent_count = int(mask.sum()), int(masks[parent].sum())
            row = {
                "sample_id": s.sample_id,
                "gate": gate,
                "parent": parent,
                "count": count,
                "parent_count": parent_count,
                "percent_parent": 100 * count / parent_count if parent_count else np.nan,
                "percent_total": 100 * count / total if total else np.nan,
                **{f"metadata:{k}": v for k, v in s.metadata.items()},
                **_stats(s.signal[mask]),
            }
            rows.append(row)
    summary = pd.DataFrame(rows)
    qcdb.check_flow_summary(summary)
    return summary


def _events(samples: Sequence[GatedSample], gate_names: Sequence[str]) -> pd.DataFrame:
    parts = []
    for s in samples:
        n = len(s.raw)
        cols = {
            "sample_id": np.repeat(s.sample_id, n),
            "event_index": np.arange(n, dtype=np.int64),
            "input_sha256": np.repeat(s.sha256, n),
            "signal_space": np.repeat(s.signal_space, n),
        }
        cols.update({f"metadata:{k}": np.repeat(str(v), n) for k, v in s.metadata.items()})
        cols.update({f"raw:{c}": s.raw[c].to_numpy(np.float32) for c in s.raw.columns})
        cols.update({f"signal:{c}": s.signal[c].to_numpy(np.float32) for c in s.signal.columns})
        cols["gate:root"] = np.ones(n, dtype=bool)
        cols.update({f"gate:{g}": s.gates.get(g, np.zeros(n, dtype=bool)) for g in gate_names})
        parts.append(pd.DataFrame(cols))
    return pd.concat(parts, ignore_index=True)


def _versions(*names: str) -> dict:
    out = {"python": platform.python_version()}
    for name in ("deliverome-analysis", "agentflow-cytometry", "flowio", "numpy", "pandas", *names):
        try:
            out[name] = md.version(name)
        except md.PackageNotFoundError:
            pass
    return out


def write_run(
    out: str | Path,
    samples: Sequence[GatedSample],
    parents: Mapping[str, str],
    *,
    producer: str,
    recipe: dict,
    recipe_sha256: str,
    signal_space: str,
    extra: Mapping[str, str | Path] | None = None,
    notes: str = "",
) -> Path:
    """Write a publishable run folder. ``parents``: gate -> parent gate (``root`` at the top).

    ``recipe`` is the human-readable gating definition saved as ``reproducibility.yaml``
    (published as ``gates.yaml``); ``recipe_sha256`` identifies the exact recipe source (the
    ``.wsp`` file, or the committed script). ``extra`` files are copied in by name.
    """
    out = Path(out)
    run = out / "run"
    run.mkdir(parents=True, exist_ok=False)
    gate_names = list(parents)
    if not samples:
        raise ValueError("no samples to write")
    for s in samples:
        missing = set(s.gates) - set(gate_names)
        if missing:
            raise ValueError(f"{s.sample_id}: gates without a parent entry {sorted(missing)}")
        if len(s.signal) != len(s.raw):
            raise ValueError(f"{s.sample_id}: signal and raw event counts differ")
    _events(samples, gate_names).to_parquet(run / "events.parquet", index=False)
    summarize(samples, parents).to_csv(run / "summary.csv", index=False)
    sheet = pd.DataFrame(
        [{"sample_id": s.sample_id, "source_uri": s.source_uri, "sha256": s.sha256, **s.metadata}
         for s in samples]
    )  # fmt: skip
    sheet.to_csv(run / "samples.csv", index=False)
    (run / "reproducibility.yaml").write_text(
        yaml.safe_dump({"producer": producer, **recipe}, sort_keys=False, allow_unicode=True)
    )
    versions = _versions()
    provenance = {
        "producer": producer,
        "inputs": [
            {
                "sample_id": s.sample_id,
                "sha256": s.sha256,
                "signal_space": s.signal_space,
                "source_uri": s.source_uri,
                "metadata": s.metadata,
                "event_count": len(s.raw),
            }
            for s in samples
        ],
        # Agentflow did not gate these runs; its vendored FlowKit applied FlowJo gates when the
        # producer is flowjo. Recorded so the index's agentflow_version column is truthful.
        "agentflow": {"version": versions.get("agentflow-cytometry"), "role": "library only"},
        "recipe_sha256": recipe_sha256,
        "samples_sha256": sha256(run / "samples.csv"),
        "signal_space": signal_space,
        "versions": versions,
        "summary": notes or "Statistics in the signal space over all events of each gate",
    }
    (run / "run.json").write_text(json.dumps(provenance, indent=2, default=str) + "\n")
    for name, src in (extra or {}).items():
        shutil.copy2(src, run / name)
    return out


# ----------------------------------------------------------------------------------------------
# FlowJo
# ----------------------------------------------------------------------------------------------


def flowjo_saved_counts(wsp: str | Path) -> dict[tuple[str, tuple[str, ...]], int]:
    """FlowJo's own count per (sample name, gate path) as saved in the workspace XML."""
    counts = {}
    for sample in ET.parse(wsp).getroot().iter("Sample"):
        node = sample.find("SampleNode")
        if node is None:
            continue

        def walk(el, path, name=node.get("name")):
            for pop in el.findall("Subpopulations/Population"):
                p = (*path, pop.get("name"))
                if pop.get("count") is not None:
                    counts[(name, p)] = int(pop.get("count"))
                walk(pop, p)

        walk(node, ())
    return counts


def patch_comp_infix(text: str) -> str:
    """Map FlowJo ``userProvidedCompInfix`` channel names (Comp-FL5-CF488-A) back to FL5-A.

    FlowJo lets a user rename compensated channels; the vendored FlowKit does not resolve those
    names and fails to find the channel. Gate geometry is untouched.
    """
    pairs = re.findall(r'data-type:name="([^"]+)"\s+userProvidedCompInfix="Comp-([^"]+)"', text)
    for raw, named in dict(pairs).items():
        text = text.replace(named, raw)
    return text


def _unique_names(ids: Sequence[tuple[str, tuple[str, ...]]]) -> dict:
    """FlowJo gate id -> index name: the gate name, or its path when the name repeats."""
    names = [g for g, _ in ids]
    out = {}
    for g, path in ids:
        out[(g, path)] = g if names.count(g) == 1 else "/".join([*path[1:], g])
    return out


def gate_flowjo(
    wsp: str | Path,
    fcs: Mapping[str, str | Path],
    *,
    group: str | None = None,
    samples: Sequence[str] | None = None,
    source_uri: Mapping[str, str] | None = None,
    metadata: Mapping[str, dict] | None = None,
) -> tuple[list[GatedSample], dict[str, str], pd.DataFrame, dict]:
    """Apply a FlowJo workspace's gates; returns samples, parents, validation table, info.

    ``fcs``: FCS file name as FlowJo knows it -> local path. ``group``: FlowJo sample group
    (default: every sample in the workspace that has a file). ``samples``: restrict to these
    FCS names. Validation compares every FlowKit population count with FlowJo's saved count.
    """
    from agentflow._vendor import flowkit as fk

    with tempfile.TemporaryDirectory(prefix="flowjo-import-") as tmp:
        tmp = Path(tmp)
        text = Path(wsp).read_text()
        patched = patch_comp_infix(text)
        (tmp / "workspace.wsp").write_text(patched)
        loaded, offset_fixed = [], []
        for name, path in fcs.items():
            try:
                loaded.append(fk.Sample(str(path), sample_id=name))
            except Exception as exc:  # CytExpert writes a data offset one byte off
                if "off by 1" not in str(exc):
                    raise
                loaded.append(fk.Sample(str(path), sample_id=name, ignore_offset_error=True))
                offset_fixed.append(name)
        ws = fk.Workspace(str(tmp / "workspace.wsp"), fcs_samples=loaded)
    ids = ws.get_sample_ids(group_name=group) if group else ws.get_sample_ids()
    ids = [i for i in ids if i in fcs and (samples is None or i in set(samples))]
    if not ids:
        raise ValueError("no workspace samples matched the FCS files")
    saved = flowjo_saved_counts(wsp)
    gated, checks, parents = [], [], {}
    for sid in ids:
        ws.analyze_samples(sample_id=sid, use_mp=False)
        gate_ids = ws.get_gate_ids(sid)
        names = _unique_names(gate_ids)
        sample = ws.get_sample(sid)
        raw = pd.DataFrame(sample.get_events(source="raw"), columns=sample.pnn_labels)
        comp = ws.get_comp_matrix(sid)
        if comp is not None:
            sample.apply_compensation(comp)
            signal = pd.DataFrame(sample.get_events(source="comp"), columns=sample.pnn_labels)
            space = "compensated"
        else:
            signal, space = raw, "raw"
        masks = {}
        for g, path in sorted(gate_ids, key=lambda gp: len(gp[1])):
            name = names[(g, path)]
            parent = "root" if len(path) == 1 else names[(path[-1], tuple(path[:-1]))]
            if parents.setdefault(name, parent) != parent:
                raise ValueError(f"gate {name} has different parents across samples")
            masks[name] = np.asarray(ws.get_gate_membership(sid, g, gate_path=path), dtype=bool)
            flowjo = saved.get((sid, (*path[1:], g)))
            ours = int(masks[name].sum())
            checks.append(
                {
                    "sample_id": sid,
                    "gate": name,
                    "parent": parent,
                    "flowkit_count": ours,
                    "flowjo_count": flowjo,
                }
            )
        meta = {"fcs_name": sid, **(metadata or {}).get(sid, {})}
        well = re.search(r"Well-([A-P]\d{1,2})", sid)
        if well and "well" not in meta:
            meta["well"] = well[1]
        if group:
            meta.setdefault("flowjo_group", group)
        gated.append(
            GatedSample(
                sample_id=Path(sid).stem,
                raw=raw,
                signal=signal,
                gates=masks,
                sha256=sha256(fcs[sid]),
                source_uri=(source_uri or {}).get(sid, str(fcs[sid])),
                metadata=meta,
                signal_space=space,
            )
        )
    checks = validate_counts(pd.DataFrame(checks))
    spaces = {s.signal_space for s in gated}
    info = {
        "comp_infix_patched": patched != text,
        "offset_error_ignored": offset_fixed,
        "signal_space": spaces.pop() if len(spaces) == 1 else "mixed",
    }
    return gated, parents, checks, info


def validate_counts(checks: pd.DataFrame) -> pd.DataFrame:
    """Add percent-of-parent columns and ``within_tolerance`` (see :data:`TOLERANCE`).

    ``checks``: one row per sample x gate with ``parent``, ``flowkit_count``, ``flowjo_count``.
    """
    c = checks.copy()
    look = c.set_index(["sample_id", "gate"])
    for who in ("flowkit", "flowjo"):
        parent = [
            look[f"{who}_count"].get((sid, par)) if par != "root" else None
            for sid, par in zip(c.sample_id, c.parent, strict=True)
        ]
        # FlowJo saves no count for all events, so top-level gates are judged on count only
        c[f"{who}_parent"] = pd.Series(parent, index=c.index, dtype=float)
    count_diff = (c.flowkit_count - c.flowjo_count).abs()
    count_ok = count_diff <= np.maximum(TOLERANCE * c.flowjo_count, TOLERANCE_EVENTS)
    c["percent_parent_diff"] = 100 * (
        c.flowkit_count / c.flowkit_parent.clip(lower=1)
        - c.flowjo_count / c.flowjo_parent.clip(lower=1)
    )
    judged = c.flowjo_parent >= _MIN_PARENT_FOR_PERCENT
    percent_ok = ~judged | (c.percent_parent_diff.abs() <= TOLERANCE_POINTS)
    c["relative_difference"] = (c.flowkit_count - c.flowjo_count) / c.flowjo_count.clip(lower=1)
    c["within_tolerance"] = c.flowjo_count.notna() & count_ok & percent_ok
    return c


def read_sample_sheet(path: str | Path) -> dict[str, dict]:
    """A sample sheet (CSV, ``fcs_name`` + metadata columns) as ``fcs_name -> metadata``.

    Values are kept verbatim as text; blank cells are dropped rather than stored as "".
    """
    sheet = pd.read_csv(path, dtype=str).fillna("")
    if "fcs_name" not in sheet or not sheet["fcs_name"].is_unique:
        raise ValueError(f"{path}: needs a unique fcs_name column")
    return {
        r.pop("fcs_name"): {k: v for k, v in r.items() if v != ""} for r in sheet.to_dict("records")
    }


def import_flowjo(
    experiment: str,
    wsp_key: str,
    *,
    group: str | None = None,
    samples: Sequence[str] | None = None,
    metadata: Mapping[str, dict] | str | Path | None = None,
    fcs_prefix: str | None = None,
    run_id: str | None = None,
    publish: bool = True,
    out: str | Path | None = None,
    record_dir: str | Path | None = None,
    accept_mismatch: str | None = None,
) -> dict:
    """Gate a FlowJo workspace on ``external`` with FlowKit, validate, and publish the run.

    FCS files are the ``.fcs`` keys beside the workspace (or under ``fcs_prefix``). Refuses to
    publish when any population disagrees with FlowJo beyond :data:`TOLERANCE` unless
    ``accept_mismatch`` gives the reason; the reason and the disagreeing populations are then
    recorded in the run (``gates.yaml`` validation block and ``validation.csv``).
    """
    from deliverome_analysis import flow
    from deliverome_analysis.s3 import external

    sheet_path = None
    if isinstance(metadata, (str, Path)):
        sheet_path, metadata = Path(metadata), read_sample_sheet(metadata)
    prefix = fcs_prefix or wsp_key.rsplit("/", 1)[0] + "/"
    keys = [k for k in external.list_keys(prefix=prefix) if k.lower().endswith(".fcs")]
    by_name = {k.rsplit("/", 1)[1]: k for k in keys}
    if samples is not None:
        by_name = {n: k for n, k in by_name.items() if n in set(samples)}
    work = Path(out) if out else Path(tempfile.mkdtemp(prefix="flowjo-run-"))
    files = work / "inputs"
    files.mkdir(parents=True, exist_ok=True)
    fcs = {}
    for name, key in by_name.items():
        fcs[name] = files / name
        if not fcs[name].exists():  # link the cached download: FCS sets run to gigabytes
            fcs[name].symlink_to(Path(external.download(key)).resolve())
    wsp = files / wsp_key.rsplit("/", 1)[1]
    shutil.copy(external.download(wsp_key), wsp)
    bucket = external.name
    gated, parents, checks, info = gate_flowjo(
        wsp,
        fcs,
        group=group,
        samples=samples,
        metadata=metadata,
        source_uri={n: f"s3://{bucket}/{k}" for n, k in by_name.items()},
    )
    bad = checks[~checks.within_tolerance]
    rel = checks.relative_difference
    validation = work / "validation.csv"
    checks.to_csv(validation, index=False)
    result = {
        "experiment": experiment,
        "wsp": wsp_key,
        "group": group,
        "samples": len(gated),
        "populations": len(checks),
        "mismatches": len(bad),
        "samples_without_sheet_row": sorted(
            s.metadata["fcs_name"]
            for s in gated
            if metadata and s.metadata["fcs_name"] not in metadata
        ),
        "median_relative_difference": float(rel.median()),
        "max_abs_relative_difference": float(rel.abs().max()),
        "validation": str(validation),
        **info,
    }
    if len(bad) and not (accept_mismatch or "").strip():
        result["status"] = "refused: FlowKit counts disagree with FlowJo"
        result["worst"] = bad.head(10).to_dict("records")
        return result
    recipe = {
        "flowjo_workspace": {
            "uri": f"s3://{bucket}/{wsp_key}",
            "sha256": sha256(wsp),
            "group": group,
        },
        "gates": [{"name": g, "parent": p} for g, p in parents.items()],
        "validation": {
            "compared_to": "population counts FlowJo saved in the workspace",
            "tolerance": (
                f"count within max({TOLERANCE:.1%} of FlowJo, {TOLERANCE_EVENTS} events) and "
                f"percent of parent within {TOLERANCE_POINTS} points"
            ),
            "max_abs_percent_parent_difference": float(checks.percent_parent_diff.abs().max()),
            "accepted_mismatch": accept_mismatch if len(bad) else None,
            "outside_tolerance_populations": bad[["sample_id", "gate"]].to_dict("records"),
            "populations": len(checks),
            "outside_tolerance": len(bad),
            "median_abs_relative_difference": float(rel.abs().median()),
            "max_abs_relative_difference": float(rel.abs().max()),
        },
        "sample_sheet": {"path": str(sheet_path), "sha256": sha256(sheet_path)}
        if sheet_path
        else None,
        "notes": {
            "comp_infix_patched": info["comp_infix_patched"],
            "offset_error_ignored": info["offset_error_ignored"],
        },
    }
    write_run(
        work / "output",
        gated,
        parents,
        producer="flowjo",
        recipe=recipe,
        recipe_sha256=sha256(wsp),
        signal_space=info["signal_space"],
        extra={"validation.csv": validation, "workspace.wsp": wsp},
        notes="FlowJo gates applied with FlowKit; counts validated against FlowJo's saved counts",
    )
    result["output"] = str(work / "output")
    if publish:
        result["published"] = flow.publish_run(
            work / "output", experiment, run_id, record_dir=record_dir
        )
    result["status"] = "published" if publish else "written"
    return result


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python -m deliverome_analysis.flow_import")
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("flowjo", help="gate a FlowJo workspace with FlowKit and publish")
    f.add_argument("--experiment", required=True, help="sprintN/sNeM-slug")
    f.add_argument("--wsp", required=True, help="external key of the .wsp")
    f.add_argument("--group", help="FlowJo sample group (default: all samples)")
    f.add_argument("--samples", nargs="*", help="only these FCS file names")
    f.add_argument("--sample-sheet", help="CSV: fcs_name + metadata columns (cell_line, ...)")
    f.add_argument("--run-id")
    f.add_argument("--record-dir", help="where the run record goes (default: experiment folder)")
    f.add_argument("--dry-run", action="store_true", help="validate and write, do not publish")
    f.add_argument(
        "--accept-mismatch",
        metavar="REASON",
        help="publish despite populations outside tolerance; the reason is recorded in the run",
    )
    a = p.parse_args(argv)
    r = import_flowjo(
        a.experiment,
        a.wsp,
        group=a.group,
        samples=a.samples,
        metadata=a.sample_sheet,
        run_id=a.run_id,
        publish=not a.dry_run,
        record_dir=a.record_dir,
        accept_mismatch=a.accept_mismatch,
    )
    print(json.dumps(r, indent=2, default=str))
    return 0 if r["status"] in {"published", "written"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
