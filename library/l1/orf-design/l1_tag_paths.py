"""Which tag-placement rung every L1 library MEMBER lands on.

This redoes `library/l1/tag-site-architecture/derive_tag_paths.py` on a different unit.
That analysis assigns a rung per **canonical gene** (2,367 of them, which is where the 515
internal-site figure comes from). But we do not order canonical genes, we order 2,700 L1
members -- 2,360 library references plus isoforms, orthologs and rescue arms -- and a
member's topology is not its gene's. SLC9D1's ordered isoform has 3 TM helices where the
canonical has 13; TRPV1's has 3 against 6; KCNN2's has 1 against 6. A rung assigned from
the canonical is therefore the wrong rung for a substantial number of the things we build.

Every member is classified from the v1 topology of its **own ordered sequence**, matched by
sequence hash through `deliverome_analysis.surfaceome.topology`, so the rung describes the protein the
construct encodes. Accession-keyed lookups cannot do this: an Ensembl-transcript ORF often
corresponds to a UniProt isoform record, so the canonical accession returns a different
protein's topology.

The ladder is unchanged from `derive_tag_paths.classify`, and deliberately so -- the point
is to change the unit, not the rule:

    ecto-N (at residue 1, or after the signal peptide) -> ecto-C -> internal

Two corrections on top of it, both of which the prior analyses name as known gaps:

* **GPI.** A GPI-anchored protein's C-terminus is removed at the omega site, so an ecto-C
  call for one is wrong -- it needs an internal site. `internal-tag-targets/README.md`
  flags exactly this ("Known gap: GPI ... BST2 and FOLR1 are the confirmed cases") and
  excludes them. Here they are rerouted instead of dropped.
* **Shedding** is recorded but never changes the rung. `appris-principal/SHEDDING_AND_TAGS.md`
  established that of 122 members whose N-terminal tag is at proteolytic risk, **zero** have
  a usable inter-TM loop -- so there is no internal fallback to route them to. Flagging it
  is the honest action; rerouting would invent an option that does not exist.

**Lane overrides** (``lane_overrides.tsv``, keyed on design id) are a documented layer on
top, like the GPI reroute: DeepTMHMM v1 stays the topology of record and ``rung_before_gpi``
stays its call, but where structures or homology contradict v1 (MFSD2A, AQP12B, GLP2R; decided
2026-11-12 by l1-reviewer-package under Becca's delegation, evidence from the PR #296 lane
check) the override sets ``placement`` (and ``signal_peptide_aa`` where given), and
``lane_override`` carries its reason. Each override names the protein it was decided on
(``sequence_sha256``); a member whose designed protein has changed stops the run.

    uv run --locked python l1_tag_paths.py
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

from deliverome_analysis.cloning import library_design as ld
from deliverome_analysis.surfaceome import topology

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
#: Authoritative GPI set from build_gpi_set.py -- NetGPI 1.1 over reference AND
#: alternative forms, unioned with curation under the CONFIDENT rule. NOT
#: gpi_anchored_genes.tsv, which is the curated input alone and misses the 27 genes
#: anchored only in a non-canonical form (NCAM1, BCAN, TFPI, CACNA2D3, ...).
GPI_FILE = HERE / "gpi_set.tsv"
SHED_FILE = HERE / "shed_form_risk.tsv"
LANE_OVERRIDES = HERE / "lane_overrides.tsv"

#: Extracellular state in the DeepTMHMM v1 alphabet. "I" is cytoplasmic, "M" a TM helix,
#: "S" the signal peptide, "B" a beta strand.
EXTRACELLULAR = "O"


def signal_peptide_end(topo: str) -> int:
    """Index one past the last signal-peptide residue, or 0 when there is none."""
    n = 0
    while n < len(topo) and topo[n] == "S":
        n += 1
    return n


def classify(topo: str) -> tuple[str, int]:
    """(rung, signal-peptide length) -- identical rule to derive_tag_paths.classify."""
    sp = signal_peptide_end(topo)
    if topo[0] == EXTRACELLULAR:
        return "N at residue 1", 0
    if sp and sp < len(topo) and topo[sp] == EXTRACELLULAR:
        return "N after SP", sp
    return ("C-terminal" if topo[-1] == EXTRACELLULAR else "internal"), sp


def longest_inter_tm_extracellular_run(topo: str) -> int:
    """Longest extracellular run that lies BETWEEN two TM helices.

    An N-terminal ectodomain is not an inter-TM loop, and counting it is the error
    SHEDDING_AND_TAGS.md records (ADGRB1's 918-residue N-terminal segment was once read as
    a loop candidate). Only runs with a TM on both sides qualify as insertion sites.
    """
    best = 0
    first_tm = topo.find("M")
    last_tm = topo.rfind("M")
    if first_tm < 0 or last_tm <= first_tm:
        return 0
    run = 0
    for i in range(first_tm, last_tm + 1):
        if topo[i] == EXTRACELLULAR:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


def _load_set(path: Path, column: str) -> set[str]:
    if not path.exists():
        return set()
    with path.open() as fh:
        return {r[column].strip() for r in csv.DictReader(fh, delimiter="\t") if r[column].strip()}


def longest_extracellular_run(topo: str) -> tuple[int, int, str]:
    """Longest extracellular run anywhere: (length, 1-based start, where it sits).

    Complements :func:`longest_inter_tm_extracellular_run`, which deliberately counts only
    runs flanked by a TM on both sides and so returns 0 for any single-pass protein. That
    is right for a multipass loop inventory and useless for the GPI-anchored single-pass
    members, whose site goes somewhere in the ectodomain rather than in a loop between two
    helices -- and reading their 0 as "no site available" is simply wrong. BST2 has a
    137-residue ectodomain, RGMB 394, GPC1 535, SPRN 127.

    ``where`` is ``n_terminal`` (before the first TM), ``inter_tm`` (between two), or
    ``c_terminal`` (after the last). With no TM at all the whole extracellular region is
    reported as ``no_tm`` -- a GPI-anchored protein where the anchor is the only membrane
    attachment, so position relative to a helix is not a meaningful question.
    """
    best = (0, 0, "none")
    first_tm, last_tm = topo.find("M"), topo.rfind("M")
    run = 0
    for i, c in enumerate(topo):
        if c == EXTRACELLULAR:
            run += 1
            if run > best[0]:
                start = i - run + 2  # 1-based start of this run
                if first_tm < 0:
                    where = "no_tm"
                elif i < first_tm:
                    where = "n_terminal"
                elif start - 1 > last_tm:
                    where = "c_terminal"
                else:
                    where = "inter_tm"
                best = (run, start, where)
        else:
            run = 0
    return best


def _load_gpi(path: Path) -> set[str]:
    """Genes whose C-terminus is removed at a GPI omega site, so ecto-C is unavailable."""
    if not path.exists():
        raise SystemExit(f"{path.name} missing -- run build_gpi_set.py first")
    with path.open() as fh:
        return {
            r["gene_symbol"].strip()
            for r in csv.DictReader(fh, delimiter="\t")
            if r["gpi_anchored"] == "1"
        }


def _load_lane_overrides(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as fh:
        return {r["design_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _load_shed(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open() as fh:
        return {r["gene_symbol"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _print_report(rows: list[dict]) -> None:
    """Placement counts by candidate set, GPI reroutes and thin internal members."""
    counts = collections.Counter(r["placement"] for r in rows)
    print(f"\n{'placement':46s} {'members':>8} {'pct':>6}")
    for p, n in counts.most_common():
        print(f"{p:46s} {n:>8} {100 * n / len(rows):>5.1f}%")

    print(
        f"\n{'placement':46s} "
        + " ".join(f"{s:>12}" for s in ("library_ref", "isoform", "ortholog"))
    )
    sets = ("library_reference", "isoform", "ortholog")
    for p, _ in counts.most_common():
        cells = [
            sum(1 for r in rows if r["placement"] == p and r["candidate_set"] == s) for s in sets
        ]
        print(f"{p:46s} " + " ".join(f"{c:>12}" for c in cells))

    rerouted = [r for r in rows if r["gpi_rerouted"]]
    print(
        f"\nGPI reroutes (ecto-C -> internal): {len(rerouted)}"
        f"  genes: {sorted({r['gene_symbol'] for r in rerouted})[:8]}"
    )
    internal = [r for r in rows if r["placement"] == "internal"]
    thin = [r for r in internal if r["longest_inter_tm_ecto_run_aa"] < 8]
    print(
        f"internal-site members: {len(internal)}; of those with no inter-TM loop >= 8 aa: {len(thin)}"
    )
    # A 0 in the inter-TM column reads as "nowhere to put a tag" and for a single-pass
    # member that is wrong, so report the ectodomain these actually have.
    for r in sorted(thin, key=lambda x: -x["longest_ecto_run_aa"]):
        print(
            f"    {r['candidate_id']:12s} {r['gene_symbol']:10s} tm={r['n_tm']:>2} "
            f"inter_tm={r['longest_inter_tm_ecto_run_aa']:>3}  "
            f"longest_ecto={r['longest_ecto_run_aa']:>4} aa at {r['longest_ecto_run_start']:>4} "
            f"({r['longest_ecto_run_position']})"
        )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", choices=ld.LIBRARIES, default=ld.LIBRARY_L1)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    out = args.out or HERE / f"{args.library}_tag_paths.tsv"

    members = ld.designed_proteins(args.manifest, library=args.library)
    report = ld.topology_report(members, tool=topology.TOOL_V1)
    gpi = _load_gpi(GPI_FILE)
    shed = _load_shed(SHED_FILE)
    print(f"{args.library}: {len(members)} members; GPI genes {len(gpi)}; shed genes {len(shed)}")

    by_cid = {r["candidate_id"]: r for r in report.to_dict("records")}
    overrides = _load_lane_overrides(LANE_OVERRIDES)
    mane_file = HERE / "mane_gtf_changes_2026-11-22.tsv"
    mane_switched: set[str] = set()
    if mane_file.exists():
        with mane_file.open() as fh:
            mane_switched = {r["design_id"] for r in csv.DictReader(fh, delimiter="\t")}
    design_id = ld.design_ids(args.manifest)
    applied: list[str] = []
    stale: list[str] = []
    rows: list[dict] = []
    no_topology: list[str] = []
    for cid, m in members.items():
        r = by_cid[cid]
        if not r["found"]:
            no_topology.append(cid)
            continue
        topo = r["per_residue_topology"]
        rung, sp = classify(topo)
        ecto_run, ecto_start, ecto_where = longest_extracellular_run(topo)
        n_tm = topo.count("M") and len(
            [1 for i, c in enumerate(topo) if c == "M" and (i == 0 or topo[i - 1] != "M")]
        )
        gene = m["gene_symbol"]

        # Split the residue-1 rung the way derive_tag_paths does: multipass is Tedman's
        # validated GPCR case; a single TM is type III, where that one helix carries the
        # whole topogenic decision and an N-terminal cassette shifts its charge balance.
        placement = rung
        if rung == "N at residue 1":
            placement = (
                "N at residue 1, multipass"
                if n_tm > 1
                else "N at residue 1, single-pass (type III)"
            )

        # GPI reroute: the omega site deletes the C-terminus, so an ecto-C call is void.
        gpi_rerouted = False
        if rung == "C-terminal" and gene in gpi:
            placement = "internal"
            gpi_rerouted = True

        lane_override = ""
        ov = overrides.get(design_id[cid])
        if ov is not None:
            if ov["sequence_sha256"] != r["sequence_sha256"]:
                # A member whose reference moved (the GTF MANE rerun, 2026-11-22) has a new
                # protein; an override decided on the old one does not carry over. It is placed
                # by the default rule and listed, for the override to be re-decided.
                if design_id[cid] not in mane_switched:
                    raise SystemExit(
                        f"{design_id[cid]}: lane override was decided on another protein"
                    )
                print(f"stale lane override (protein changed by the MANE rerun): {design_id[cid]}")
                stale.append(design_id[cid])
                ov = None
        if ov is not None:
            placement, lane_override = ov["placement"], ov["reason"]
            if ov["signal_peptide_aa"]:
                sp = int(ov["signal_peptide_aa"])
            applied.append(design_id[cid])

        s = shed.get(gene, {})
        rows.append(
            {
                "candidate_id": cid,
                "library": args.library,
                "candidate_set": m["candidate_set"],
                "species": m["species"],
                "gene_symbol": gene,
                "placement": placement,
                "rung_before_gpi": rung,
                "gpi_rerouted": int(gpi_rerouted),
                "lane_override": lane_override,
                "n_tm": n_tm,
                "signal_peptide_aa": sp,
                "length_aa": len(topo),
                "longest_inter_tm_ecto_run_aa": longest_inter_tm_extracellular_run(topo),
                # The ectodomain measure, for members whose site is not an inter-TM loop.
                "longest_ecto_run_aa": ecto_run,
                "longest_ecto_run_start": ecto_start,
                "longest_ecto_run_position": ecto_where,
                # Informational only -- see the module docstring on why shedding does not
                # change the rung.
                "shed_risk": s.get("shed_severity", ""),
                "shed_evidence": s.get("shed_evidence", ""),
                "topology_cohort": r["cohort"],
                "topology_version": r["topology_version"],
                "sequence_sha256": r["sequence_sha256"],
            }
        )

    if args.library == ld.LIBRARY_L1:
        # An override for a member since retired (e.g. the 2026-11-22 isoform redraw) no longer
        # applies; one that names no member at all is still an error.
        with MANIFEST.open() as fh:
            retired = {
                r["candidate_id"] for r in csv.DictReader(fh, delimiter="\t") if r["retired_reason"]
            }
        gone = sorted(set(overrides) - set(applied) - retired - set(stale))
        if gone:
            raise SystemExit(f"lane overrides name no active L1 member: {gone}")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {out.name}: {len(rows)} members classified")
    if no_topology:
        print(
            f"  {len(no_topology)} with no v1 topology for their exact sequence: {no_topology[:5]}"
        )

    _print_report(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
