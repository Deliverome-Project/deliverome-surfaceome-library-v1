"""The surface-protein ORF library, from gene list to ordered proteoforms, as one pipeline.

Every stage is listed here in run order with what it decides, the script that decides it,
where its inputs come from, and whether it is **verified** (rerun from pinned inputs and
compared with its record) or **pinned** (its inputs are lost or live, so its published
output is the frozen record). Freeze-then-verify (Becca, 2026-11-01): nothing is rebuilt
from today's sources; each rerunnable stage must reproduce what the library actually used.

Inputs that are not already immutable S3 keys live in one dated snapshot,
``processed://sprint4/s4e9-orf-design-fg-panel/library-design/pipeline-inputs/2026-11-01/``,
with a SHA-256 per file in ``library_pipeline_inputs.json``. The panel draws keep their own
snapshot (``scratch/becca/appris-principal/panel_inputs.json``).

    uv run --locked python library_pipeline.py plan              # the stage table
    uv run --locked python library_pipeline.py verify            # rerun and compare all
    uv run --locked python library_pipeline.py verify --stage reference
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[2]
AP = REPO / "scratch/becca/appris-principal"
INPUTS = HERE / "library_pipeline_inputs.json"
ADDITIONS = HERE / "gene_list_additions.tsv"
EXCLUSIONS = HERE / "gene_list_exclusions.tsv"
S1E3 = REPO / "analysis/sprint1/s1e3-surface-library-filters"
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stage_inputs(names: list[str], dest: Path) -> dict[str, Path]:
    """Fetch named inputs from the snapshot (or their pinned S3 key) and verify each."""
    from deliverome_analysis import s3

    spec = json.loads(INPUTS.read_text())
    out = {}
    for name in names:
        entry = spec["files"][name]
        bucket, _, key = entry["source"].partition(":")
        if bucket == "snapshot":
            bucket, key = "processed", spec["prefix"] + name
        path = dest / name
        path.parent.mkdir(parents=True, exist_ok=True)
        getattr(s3, bucket).download(key, path)
        if sha256(path) != entry["sha256"]:
            raise SystemExit(f"{name}: SHA-256 does not match library_pipeline_inputs.json")
        out[name] = path
    return out


def record(key: str, dest: Path) -> Path:
    from deliverome_analysis.s3 import processed

    return Path(processed.download(PREFIX + key, dest))


def run(cmd: list[str], cwd: Path) -> None:
    env = {**os.environ, "PYTHONPATH": str(AP)}
    proc = subprocess.run([sys.executable, *cmd], cwd=cwd, env=env, capture_output=True, text=True)
    if proc.returncode:
        raise SystemExit(f"{' '.join(cmd[:2])} failed:\n{proc.stderr[-2000:]}")


def same_bytes(a: Path, b: Path) -> bool:
    return Path(a).read_bytes().replace(b"\r\n", b"\n") == Path(b).read_bytes().replace(
        b"\r\n", b"\n"
    )


# --------------------------------------------------------------------------- verifiers


def verify_gpi_set(work: Path) -> str:
    out = work / "gpi_set.tsv"
    run([str(HERE / "build_gpi_set.py"), "--out", str(out)], work)
    assert same_bytes(out, HERE / "gpi_set.tsv"), "gpi_set.tsv does not reproduce"
    return "gpi_set.tsv rebuilt byte for byte from committed NetGPI predictions + curation"


GENE_LIST_INPUTS = {
    # name in library_pipeline_inputs.json -> path inside the staged S1E3 folder
    "deep_dives_1.3.0.tar.gz": None,
    "deep-dive-preset-membership.tsv": None,
    "genomewide_expression_tpm.csv": "external/expression_tpm/genomewide_expression_tpm.csv",
    "genomewide_expression_tpm_rederived.csv": (
        "external/expression_tpm/genomewide_expression_tpm_rederived.csv"
    ),
    "hpa_celline_per_gene.csv": "external/hpa_rna_celline/hpa_celline_per_gene.csv",
    "mygene_response.json": "external/mygene_symbol_resolution/mygene_response.json",
    "mygene_response_full.json": "external/mygene_symbol_resolution/full/mygene_response.json",
    "coreceptor_calls.tsv": "coreceptor_calls.tsv",
    "oe_calls_136.tsv": "oe_calls_136.tsv",
    "willow_screen_rescue_genes.json": "willow_screen_rescue_genes.json",
    "willow_screen_rescue_isoforms.json": "willow_screen_rescue_isoforms.json",
    "coreceptor_rationales.json": "coreceptor_rationales.json",
}


def run_gene_filter(work: Path) -> set[str]:
    """The S1E3 cohort filter over accessible-surfaceome release 1.3.0, offline.

    Per-gene records come from the release's deep-dive tarball (each file checked against
    the release manifest's SHA-256) and are compacted exactly as fetch_cohort.py compacts
    the API's; the low-literature UniProt gate reads the release's ``low_lit_uniprot``
    badge, which equals (UniProt bit and n_papers_found < 100), the pair the gate tests.
    MyGene symbol resolution replays from its pinned responses, and the run fails if the
    resolver had to ask the network.
    """
    import csv
    import tarfile

    got = stage_inputs(list(GENE_LIST_INPUTS), work / "in")
    s1 = work / "s1e3"
    s1.mkdir()
    for py in S1E3.glob("*.py"):
        shutil.copy(py, s1 / py.name)
    for name, rel in GENE_LIST_INPUTS.items():
        if rel:
            (s1 / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(got[name], s1 / rel)
    rel_dir = work / "release"
    with tarfile.open(got["deep_dives_1.3.0.tar.gz"]) as tar:
        tar.extractall(
            rel_dir,
            members=[
                m
                for m in tar.getmembers()
                if m.name == "manifest.tsv"
                or (m.name.endswith(".json") and not m.name.endswith(".evidence.json"))
            ],
            filter="data",
        )
    sys.path.insert(0, str(s1))
    import fetch_cohort

    with (rel_dir / "manifest.tsv").open() as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    with (s1 / "cohort.jsonl").open("w") as out:
        for m in manifest:
            path = rel_dir / "genes" / f"{m['gene_symbol']}.json"
            assert sha256(path) == m["json_sha256"], f"{path.name}: not the release's bytes"
            out.write(json.dumps(fetch_cohort.compact(json.loads(path.read_text()))) + "\n")
    with got["deep-dive-preset-membership.tsv"].open() as fh:
        rows = [
            {"symbol": r["gene_symbol"], "db": int(r["low_lit_uniprot"] == "1")}
            for r in csv.DictReader(fh, delimiter="\t")
        ]
    (s1 / "catalog_rows.json").write_text(json.dumps(rows))
    caches = {
        n: sha256(s1 / GENE_LIST_INPUTS[n])
        for n in ("mygene_response.json", "mygene_response_full.json")
    }
    run([str(s1 / "build_expression_join.py")], s1)
    run([str(s1 / "build_expression_join.py"), "--full"], s1)
    for n, h in caches.items():
        assert sha256(s1 / GENE_LIST_INPUTS[n]) == h, f"{n}: resolver went to the network"
    run([str(s1 / "final_library.py")], s1)
    with (s1 / "library_membership.tsv").open() as fh:
        return {r["gene"] for r in csv.DictReader(fh, delimiter="\t")}


def verify_gene_list(work: Path) -> str:
    import csv

    passed = run_gene_filter(work)
    with EXCLUSIONS.open() as fh:
        excluded = {r["gene"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
    library = set(json.loads(gene_list(work / "genes").read_text()))
    unexplained = sorted((passed - set(excluded)) ^ library)
    assert not unexplained, f"filter output differs from the gene list for {unexplained}"
    return (
        f"S1E3 filter over surfaceome release 1.3.0 passes {len(passed)} genes; minus "
        f"{len(excluded)} recorded exclusion(s) ({', '.join(excluded)}) = the {len(library)} "
        "library genes exactly"
    )


def gene_list(dest: Path) -> Path:
    """The library's gene list: the frozen 2,368 root plus ``gene_list_additions.tsv``."""
    import csv

    root = stage_inputs(["library_members.json"], dest / "root")["library_members.json"]
    genes = json.loads(root.read_text())
    with ADDITIONS.open() as fh:
        added = [r["gene"] for r in csv.DictReader(fh, delimiter="\t")]
    clash = sorted(set(added) & set(genes))
    if clash:
        raise SystemExit(f"gene_list_additions.tsv repeats root genes: {clash}")
    out = dest / "library_members.json"
    out.write_text(json.dumps(sorted(genes + added)))
    return out


def run_reference(work: Path) -> Path:
    """Rerun the reference ladder offline over the gene list; returns the output folder."""
    names = [
        "appris_homo_sapiens_e115v50.txt",
        "Homo_sapiens.GRCh38.pep.all.fa.gz",
        "uniprot_human_reviewed_isoforms.tsv.gz",
        "peptideatlas_build607_pools.json",
        "Ribo-seq_ORFs.all.bed",
        "Ribo-seq_ORFs.all.faa",
        "Ribo-seq_ORFs.primary.bed",
        "Ribo-seq_ORFs.primary.faa",
        "Ribo-seq_ORFs.comprehensive.bed",
        "Ribo-seq_ORFs.comprehensive.faa",
        "human_canonical_non_hla.fasta",
        "ensg_to_uniprot.tsv",
        "deeptmhmm_human_canonical.tsv",
        "deeptmhmm_human_isoforms.tsv",
        "library_members.json",
        "reference_overrides.tsv",
        "library_reference_missing_topology.3line",
        "topology_shortlist.3line",
    ]
    got = stage_inputs(names, work / "in")
    got["library_members.json"] = gene_list(work / "genes")
    cache, ribo, repo = work / "cache", work / "riboseq", work / "repo"
    for n in names[:4]:
        (cache / n).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(got[n], cache / n)
    for n in names[4:10]:
        ribo.mkdir(exist_ok=True)
        shutil.copy(got[n], ribo / n)
    layout = {
        "human_canonical_non_hla.fasta": "data/external/deeptmhmm_surfaceome_inputs",
        "ensg_to_uniprot.tsv": "data/external/uniprot_ensembl_xrefs",
        "deeptmhmm_human_canonical.tsv": "data/processed/deeptmhmm",
        "deeptmhmm_human_isoforms.tsv": "data/processed/deeptmhmm",
    }
    for n, sub in layout.items():
        (repo / sub).mkdir(parents=True, exist_ok=True)
        shutil.copy(got[n], repo / sub / n)
    for n in names[15:]:
        shutil.copy(got[n], work / n)
    out = work / "out"
    out.mkdir()
    run(
        [
            str(AP / "select_reference.py"),
            "--offline",
            "--surfaceome-repo",
            str(repo),
            "--cache-dir",
            str(cache),
            "--riboseq-dir",
            str(ribo),
            "--library",
            str(got["library_members.json"]),
            "--out-dir",
            str(out),
        ],
        work,
    )
    return out


def verify_reference(work: Path) -> str:
    out = run_reference(work)
    for name in ("reference_selection.tsv", "reference_disagreements.tsv"):
        assert same_bytes(out / name, record(f"library-design/{name}", work / f"rec_{name}")), (
            f"{name} does not reproduce its record"
        )
    return "reference_selection.tsv and reference_disagreements.tsv reproduce byte for byte"


def verify_panels(work: Path) -> str:
    rej = HERE / "panel_rejections.tsv"
    for kind in ("isoform", "ortholog"):
        run([str(AP / f"sample_{kind}_panel.py"), "--snapshot", "--original"], work)
        spec = json.loads((AP / "panel_inputs.json").read_text())
        assert sha256(work / f"{kind}_panel.tsv") == spec["reproduces"][f"{kind}_panel.tsv"], (
            f"{kind}: original draw does not reproduce"
        )
        run([str(AP / f"sample_{kind}_panel.py"), "--snapshot", "--reject", str(rej)], work)
        assert same_bytes(work / f"{kind}_panel.tsv", HERE / f"{kind}_panel.tsv"), (
            f"{kind}: redraw does not reproduce the committed panel"
        )
    return "both panels: original draw and rule-based redraw reproduce byte for byte"


def run_candidate_cds(work: Path, reference: Path | None = None) -> dict:
    """Resolve every candidate's CDS from the panel snapshot and a reference table (default:
    the published record)."""
    snap = json.loads((AP / "panel_inputs.json").read_text())
    from deliverome_analysis.s3 import processed

    for name in ("isoform_cohorts.tsv", "ortholog_repick_mouse.tsv", "ortholog_repick_cyno.tsv"):
        processed.download(snap["prefix"] + name, work / name)
    if reference is None:
        record("library-design/reference_selection.tsv", work / "reference_selection.tsv")
    else:
        shutil.copy(reference, work / "reference_selection.tsv")
    run([str(AP / "resolve_candidate_cds.py")], work)
    return json.loads((work / "candidate_cds.json").read_text())


def verify_candidate_cds(work: Path) -> str:
    new = run_candidate_cds(work)
    rec = json.loads(record("library-candidates/candidate_cds.json", work / "rec.json").read_text())
    differ = sorted(k for k in rec if k not in new or new[k]["cds"] != rec[k]["cds"])
    assert not differ, f"designed CDS does not reproduce for {len(differ)}: {differ[:5]}"
    # An ORF without a stop reads through into the vector; Ensembl omits it from some
    # records, which the resolver now repairs (stop_added). None may slip through.
    no_stop = sorted(k for k, v in rec.items() if v["cds"][-3:] not in ("TAA", "TAG", "TGA"))
    assert not no_stop, f"{len(no_stop)} designed CDS end without a stop codon: {no_stop[:5]}"
    extra = sorted(set(new) - set(rec))
    return (
        f"all {len(rec)} recorded designed CDS reproduce and end in a stop codon "
        f"({sum('stop_added' in v for v in rec.values())} appended); "
        f"extra in rerun: {extra or 'none'}"
        " (iso:MICA:Q29983 is MICA's selected reference, so the dedupe stage retires it)"
    )


def verify_constructs(work: Path) -> str:
    """Re-domesticate every designed CDS and re-assemble every vector; compare with the
    records (about a minute)."""
    import csv

    sys.path.insert(0, str(AP))
    import assemble_candidates
    import candidate_inputs
    import domesticate_candidates
    from Bio import SeqIO

    from deliverome_analysis import seqtools

    cds = json.loads(record("library-candidates/candidate_cds.json", work / "c.json").read_text())
    con = json.loads(
        record("library-candidates/candidate_constructs.json", work / "k.json").read_text()
    )
    # The tranche manifest is the record; the pre-tranche candidate_manifest.tsv is a
    # historical intermediate (positional ids, and it predates MICA's repoint).
    with (HERE / "candidate_manifest_tranches.tsv").open() as fh:
        rows = {r["key"]: r for r in csv.DictReader(fh, delimiter="\t")}
    differ = [k for k in cds if domesticate_candidates.domesticate_one(cds[k]) != con.get(k)]
    assert not differ, f"{len(differ)} ORF(s) do not re-domesticate to the record: {differ[:5]}"
    bb = str(SeqIO.read(candidate_inputs.backbone(), "genbank").seq).upper()
    bb_sites = seqtools.find_sites(bb, assemble_candidates.ENZYMES)
    bad = []
    for key, row in rows.items():
        if key not in con:
            continue
        _, built, _ = assemble_candidates.assemble_one(
            bb, bb_sites, row["candidate_id"], key, con[key]
        )
        if any(str(v) != row[c] for c, v in built.items() if c in row and c != "candidate_id"):
            bad.append(key)
    assert not bad, f"{len(bad)} vector row(s) do not reproduce: {bad[:5]}"
    return (
        f"all {len(cds)} designed CDS re-domesticate to the record and "
        f"{sum(k in con for k in rows)} vectors re-assemble identically"
    )


def verify_cassettes(work: Path) -> str:
    run([str(HERE / "l1_swap_risk.py"), "--verify"], HERE)
    return "every issued linker-only cassette encoding replays from scratch, batch by batch"


def verify_giants(work: Path) -> str:
    """Rebuild the giant truncations and their ORFs from the pinned UniProt features."""
    import gzip

    got = stage_inputs(["giant_features.json"], work / "in")
    feats = json.loads(got["giant_features.json"].read_text())
    from deliverome_analysis.s3 import external

    table = external.download(
        "uniprot/2026-10-11/uniprot_human_reviewed.tsv.gz", work / "uniprot.tsv.gz"
    )
    with gzip.open(table, "rt") as fh:
        col = next(fh).rstrip("\n").split("\t").index("Sequence")
        seqs = {f[0]: f[col] for f in (line.rstrip("\n").split("\t") for line in fh)}
    off = [g for g, f in feats.items() if seqs.get(f["accession"]) != f["sequence"]]
    assert not off, f"features not on the pinned UniProt sequence: {off[:5]}"
    shutil.copy(AP / "giant_truncation_spec.tsv", work / "giant_truncation_spec.tsv")
    shutil.copy(got["giant_features.json"], work / "giant_features.json")
    run([str(AP / "build_giant_truncations.py")], work)
    assert same_bytes(work / "giant_truncations.tsv", AP / "giant_truncations.tsv"), (
        "giant_truncations.tsv does not reproduce"
    )
    run([str(AP / "build_giant_orfs.py")], work)
    rec = record("library-candidates/giant_constructs.json", work / "rec_giant.json")
    assert json.loads((work / "giant_constructs.json").read_text()) == json.loads(
        rec.read_text()
    ), "giant_constructs.json does not reproduce"
    return f"{len(feats)} genes' features pin to UniProt 2026_03; truncations and ORFs reproduce"


def verify_rescue(work: Path) -> str:
    """Re-select the expression-rescue arms, take their CDS from Ensembl 115 files, rebuild."""
    names = [
        "deeptmhmm_human_canonical.tsv",
        "deeptmhmm_human_isoforms.tsv",
        "library_reference_missing_topology.3line",
        "library_members.json",
    ]
    got = stage_inputs(names, work / "in")
    repo = work / "repo/data/processed/deeptmhmm"
    repo.mkdir(parents=True)
    for n in names[:2]:
        shutil.copy(got[n], repo / n)
    shutil.copy(got[names[2]], work / names[2])
    snap = json.loads((AP / "panel_inputs.json").read_text())
    from deliverome_analysis.s3 import processed

    processed.download(snap["prefix"] + "reference_selection.tsv", work / "reference_selection.tsv")
    run(
        [
            str(AP / "select_expression_rescue_arms.py"),
            "--surfaceome-repo",
            str(work / "repo"),
            "--library",
            str(got["library_members.json"]),
        ],
        work,
    )
    arms = record("library-design/expression_rescue_arms.tsv", work / "rec_arms.tsv")
    assert same_bytes(work / "expression_rescue_arms.tsv", arms), "rescue arms do not reproduce"
    run([str(AP / "fetch_rescue_cds.py"), "--offline"], work)
    assert json.loads((work / "rescue_cds.json").read_text()) == json.loads(
        (AP / "rescue_cds.json").read_text()
    ), "rescue CDS from Ensembl 115 files differ from the record"
    run([str(AP / "build_rescue_constructs.py")], work)
    assert same_bytes(work / "rescue_constructs.tsv", AP / "rescue_constructs.tsv"), (
        "rescue_constructs.tsv does not reproduce"
    )
    return "25 arms re-selected; CDS from Ensembl 115 GTF + cds files; constructs reproduce"


def verify_isoform_universe(work: Path) -> str:
    import reproduce_isoform_universe as riu

    return riu.verify_isoform_universe(work)


def verify_ortholog_repick(work: Path) -> str:
    import reproduce_isoform_universe as riu

    return riu.verify_ortholog_repick(work)


def verify_universe_coverage(work: Path) -> str:
    import check_universe_coverage as cuc

    out, summary = work / "universe_coverage.tsv", work / "universe_coverage_summary.json"
    s = cuc.run(out, summary)
    assert same_bytes(out, HERE / out.name), "universe_coverage.tsv differs"
    assert same_bytes(summary, HERE / summary.name), "universe_coverage_summary.json differs"
    return (
        f"all {sum(s['recorded_pairs_checked'].values())} recorded ortholog pairs are 1:1 "
        f"high-confidence in Compara 115; {s['library_genes_in_april_set']} of "
        f"{s['library_genes']} library genes were in the April query set "
        f"({s['mouse']['eligible_not_listed']} mouse / {s['cyno']['eligible_not_listed']} cyno "
        "eligible genes never queried; recorded, not repaired)"
    )


def verify_tranches(work: Path) -> str:
    run([str(HERE / "dedupe_library_members.py"), "--check"], HERE)
    return "no two active members share a designed protein"


def verify_topology(work: Path) -> str:
    fasta = work / "gap.fasta"
    run([str(HERE / "publish_topology_supplement.py"), "fasta", "--out", str(fasta)], HERE)
    assert fasta.read_text() == "", "active members without a DeepTMHMM v1 topology remain"
    return "every active L1 member has a v1 topology in the pinned snapshot"


def verify_tag_paths(work: Path) -> str:
    for lib in ("l1", "l2"):
        out = work / f"{lib}_tag_paths.tsv"
        run([str(HERE / "l1_tag_paths.py"), "--library", lib, "--out", str(out)], HERE)
        assert same_bytes(out, HERE / f"{lib}_tag_paths.tsv"), f"{lib}_tag_paths.tsv differs"
    return "l1_tag_paths.tsv and l2_tag_paths.tsv reproduce byte for byte"


# --------------------------------------------------------------------------- the stages


@dataclass
class Stage:
    id: str
    decides: str
    script: str
    inputs: str
    status: str  # verified | pinned
    why_pinned: str = ""
    verify: Callable[[Path], str] | None = None


STAGES = [
    Stage(
        "gene_list",
        "which genes are in the library: S1E3 filter over surfaceome release 1.3.0, minus "
        "recorded exclusions (gene_list_exclusions.tsv)",
        "analysis/sprint1/s1e3-surface-library-filters/final_library.py (offline)",
        "release 1.3.0 deep dives + preset table (Zenodo 10.5281/zenodo.23445521, on "
        "external://), expression TPM tables, pinned MyGene responses, curated S1E3 calls",
        "verified",
        verify=verify_gene_list,
    ),
    Stage(
        "gpi_set",
        "which genes are GPI-anchored, per form (the GPI set of record)",
        "analysis/sprint4/s4e9-orf-design-fg-panel/build_gpi_set.py",
        "NetGPI 1.1 predictions (committed) + curated gpi_anchored_genes.tsv",
        "verified",
        verify=verify_gpi_set,
    ),
    Stage(
        "reference",
        "each gene's reference proteoform: MANE+APPRIS agree -> the single "
        "source -> observed evidence (PeptideAtlas build 607, GENCODE Ribo-seq ORFs) -> MANE; "
        "named overrides in reference_overrides.tsv",
        "scratch/becca/appris-principal/select_reference.py --offline",
        "APPRIS e115v50, Ensembl 115, UniProt 2026_03 (S3); PeptideAtlas + Ribo-seq + "
        "accessible-surfaceome f9c73691 files (snapshot)",
        "verified",
        verify=verify_reference,
    ),
    Stage(
        "isoform_universe",
        "every alternate form of each gene, with MS / Ribo-seq / APPRIS "
        "evidence; tiered into cohorts",
        "isoform_universe.py -> openprot_evidence.py -> isoform_appris.py -> "
        "shortlist_for_topology.py + DeepTMHMM v1 -> apply_topology.py -> isoform_cohorts.py",
        "each script at its producing commit (reproduce_isoform_universe.py); OpenProt 2.2, "
        "APPRIS 2026_06.v50 full scores, PeptideAtlas build 607 replay, surfaceome f9c73691 "
        "(all on external://)",
        "verified",
        "",
        verify_isoform_universe,
    ),
    Stage(
        "ortholog_repick",
        "each gene's mouse and cyno ortholog proteoform",
        "scratch/becca/appris-principal/repick_orthologs.py",
        "rebuilt universe table, NetGPI ortholog calls, Ensembl 115 + UniProt 2026_03 "
        "ortholog files, surfaceome f9c73691 cohorts (repick_orthologs.py @c4dbe752)",
        "verified",
        "",
        verify_ortholog_repick,
    ),
    Stage(
        "universe_coverage",
        "which library genes the isoform and ortholog universes could draw from: the "
        "universes were built from the April 2026 surfaceome_expressed query set, not the "
        "library gene list, so 594 library genes were never offered to them",
        "analysis/sprint4/s4e9-orf-design-fg-panel/check_universe_coverage.py",
        "Ensembl Compara 115 homology dumps (external://), accessible-surfaceome f9c73691 "
        "query set and ortholog lists, Ensembl 115 peptides, UniProt 2026_03 isoforms",
        "verified",
        verify=verify_universe_coverage,
    ),
    Stage(
        "panels",
        "the 150-isoform and 200-ortholog panels: seeded stratified draws, "
        "rejected picks replaced by continuing the draw",
        "sample_isoform_panel.py / sample_ortholog_panel.py --snapshot; redraw_panels.py",
        "panel_inputs.json snapshot (13 files) + panel_rejections.tsv",
        "verified",
        verify=verify_panels,
    ),
    Stage(
        "rescue_arms",
        "25 expression-rescue arms (Igk leader / recoded), their CDS and constructs",
        "select_expression_rescue_arms.py, fetch_rescue_cds.py --offline, "
        "build_rescue_constructs.py (seed 20261011)",
        "Tedman 2025 screen (S3), Ensembl 115 GTF canonical tag + cds FASTA, topology tables",
        "verified",
        verify=verify_rescue,
    ),
    Stage(
        "giants",
        "truncations of the oversized genes and their ORFs",
        "build_giant_truncations.py (curated giant_truncation_spec.tsv), build_giant_orfs.py",
        "UniProt 2026_03 features (REST fetch of 2026-10-11, pinned on external://)",
        "verified",
        verify=verify_giants,
    ),
    Stage(
        "candidate_cds",
        "each candidate's designed CDS, verified to encode its protein",
        "scratch/becca/appris-principal/resolve_candidate_cds.py (candidate_inputs.py)",
        "Ensembl 115 CDS/pep, UniProt 2026-10-11 (S3); reference record; panel snapshot",
        "verified",
        verify=verify_candidate_cds,
    ),
    Stage(
        "constructs",
        "domesticated ORFs (BbsI/Esp3I/BsaI) and p0039 vectors",
        "domesticate_candidates.py, assemble_candidates.py (seed 20261011)",
        "candidate_cds.json, p0039 backbone (S3)",
        "verified",
        verify=verify_constructs,
    ),
    Stage(
        "tranches",
        "L1/L2 membership; duplicates retired on the designed protein",
        "split_library_tranches.py (once), dedupe_library_members.py, redraw_panels.py",
        "candidate_manifest_tranches.tsv (committed record)",
        "verified",
        verify=verify_tranches,
    ),
    Stage(
        "topology",
        "DeepTMHMM v1 topology for every active member's designed protein",
        "publish_topology_supplement.py + modal deeptmhmm_v1_app.py",
        "external://accessible-surfaceome/topology/2026-11-01/topology.parquet",
        "verified",
        verify=verify_topology,
    ),
    Stage(
        "tag_paths",
        "each member's tag placement rung",
        "l1_tag_paths.py",
        "topology snapshot, gpi_set.tsv, shed_form_risk.tsv",
        "verified",
        verify=verify_tag_paths,
    ),
    Stage(
        "cassettes",
        "each C-terminal L1 member's tag cassette: SpyTag003 and P-ALFA-P on p003's codons, "
        "GGGGS linkers recoded so the terminal 15-nt block is unique; issued in batches, pinned",
        "l1_swap_risk.py (seed 20261029); l1_cassette_encodings.tsv records each batch",
        "member CDS, tag paths, p003 cassette codons",
        "verified",
        verify=verify_cassettes,
    ),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    v = sub.add_parser("verify")
    v.add_argument("--stage", choices=[s.id for s in STAGES])
    args = ap.parse_args(argv)

    if args.cmd == "plan":
        for n, s in enumerate(STAGES, 1):
            print(
                f"{n:>2}. {s.id:<18} [{s.status}] {s.decides}\n      {s.script}\n"
                f"      inputs: {s.inputs}"
                + (f"\n      pinned: {s.why_pinned}" if s.why_pinned else "")
            )
        return 0

    failed = 0
    for s in STAGES:
        if args.stage and s.id != args.stage:
            continue
        if s.verify is None:
            print(f"-- {s.id:<18} pinned: {s.why_pinned}")
            continue
        work = Path(tempfile.mkdtemp(prefix=f"pipeline-{s.id}-"))
        try:
            print(f"OK {s.id:<18} {s.verify(work)}")
        except (AssertionError, SystemExit) as exc:
            failed += 1
            print(f"FAIL {s.id:<16} {exc}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
