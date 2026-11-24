"""The surface-protein ORF library, from gene list to ordered proteoforms, as one pipeline.

Every stage is listed here in run order with what it decides, the script that decides it,
where its inputs come from, and whether it is **verified** (rerun from pinned inputs and
compared with its record) or **pinned** (its inputs are lost or live, so its published
output is the frozen record). Freeze-then-verify (Becca, 2026-11-01): nothing is rebuilt
from today's sources; each rerunnable stage must reproduce what the library actually used.

Inputs that are not already immutable S3 keys live in one dated snapshot,
``processed://sprint4/s4e9-orf-design-fg-panel/library-design/pipeline-inputs/2026-11-01/``,
with a SHA-256 per file in ``library_pipeline_inputs.json``. The panel draws keep their own
snapshot (``library/l1/proteoform-selection/panel_inputs.json``).

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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[2]
AP = REPO / "library/l1/proteoform-selection"
INPUTS = HERE / "library_pipeline_inputs.json"
ADDITIONS = HERE / "gene_list_additions.tsv"
EXCLUSIONS = HERE / "gene_list_exclusions.tsv"
CORECEPTOR_ADDITIONS = HERE / "coreceptor_rule_additions_proposed.tsv"
GAP_ADDITIONS = HERE / "gap_gene_additions_proposed.tsv"
S1E3 = REPO / "library/l1/surface-filters"
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
#: The reference ladder's current record (under PREFIX). 2026-11-12: the five named
#: references of apply_reference_fixes.py; the earlier record stays at library-design/.
#: Second round (LMLN, DISP2, TSPAN10 named) at its own key; the first round's is kept.
#: 2026-11-21: plus the co-receptor additions (add_library_genes.py --publish); every earlier row
#: unchanged, and the 2026-11-12 record stays at its key.
REFERENCE_RECORD = "library-design/reference-selection/2026-11-23-pre-freeze/"


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
    # The co-receptor gate restated on 2026-11-21 (Becca: every required partner >= 10 TPM in
    # >= 2 of the 42 S2E18 shortlist lines, replacing S1E3's median >= 5 TPM). It admits genes
    # that pass every other S1E3 gate; coreceptor_host_line.py derives them from the same
    # filter rerun and writes them to coreceptor_rule_additions_proposed.tsv. It drops none.
    with CORECEPTOR_ADDITIONS.open() as fh:
        restated = {r["gene"] for r in csv.DictReader(fh, delimiter="\t")}
    assert not restated & passed, f"restated additions already pass S1E3: {restated & passed}"
    # Genes Becca added from the surface-proteomics challenge set (2026-11-20 and 2026-11-22),
    # taking survey detection as surface evidence in place of the catalog's; they pass the
    # construct gates but not S1E3's surface-evidence gate, so they are listed, not derived.
    with GAP_ADDITIONS.open() as fh:
        gap = {r["gene"] for r in csv.DictReader(fh, delimiter="\t")}
    library = set(json.loads(gene_list(work / "genes").read_text()))
    unexplained = sorted(((passed | restated | gap) - set(excluded)) ^ library)
    assert not unexplained, f"filter output differs from the gene list for {unexplained}"
    return (
        f"S1E3 filter over surfaceome release 1.3.0 passes {len(passed)} genes; plus "
        f"{len(restated)} the restated co-receptor gate admits; plus {len(gap)} surface-proteomics "
        f"additions; minus {len(excluded)} recorded "
        f"exclusion(s) ({', '.join(excluded)}) = the {len(library)} library genes exactly"
    )


def additions(withdrawn: bool = False) -> list[str]:
    """Genes in ``gene_list_additions.tsv``: the active ones, or with ``withdrawn=True`` the
    ones a later decision took off the gene list (dated in the ``withdrawn`` column). A
    withdrawn gene keeps its row, and its member keeps its records, retired in the manifest."""
    import csv

    with ADDITIONS.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [r["gene"] for r in rows if bool((r.get("withdrawn") or "").strip()) == withdrawn]


def gene_list(dest: Path) -> Path:
    """The library's gene list: the frozen 2,368 root plus the active ``gene_list_additions.tsv`` rows."""

    root = stage_inputs(["library_members.json"], dest / "root")["library_members.json"]
    genes = json.loads(root.read_text())
    added = additions()
    clash = sorted(set(added) & set(genes))
    if clash:
        raise SystemExit(f"gene_list_additions.tsv repeats root genes: {clash}")
    out = dest / "library_members.json"
    out.write_text(json.dumps(sorted(genes + added)))
    return out


def run_reference(work: Path, offline: bool = True) -> Path:
    """Rerun the reference ladder over the gene list; returns the output folder.

    Offline (the default) reads only the pinned PeptideAtlas cache and fails on an accession it
    lacks. ``offline=False`` is for extending that cache for genes added to the list
    (``add_library_genes.py --fetch-peptides``): missing accessions are queried and appended
    to ``work/cache``, which is then published as a new dated snapshot and repinned."""
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
        # MANE Select read from the Ensembl 115 GTF (derive_mane_select.py), 2026-11-22.
        "mane_select_ensembl115.tsv",
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
            *(["--offline"] if offline else []),
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
        assert same_bytes(out / name, record(REFERENCE_RECORD + name, work / f"rec_{name}")), (
            f"{name} does not reproduce its record"
        )
    # The chosen protein per gene, by hash: the selection table carries source and length but
    # not the sequence. Pinned in the repo beside select_reference.py (2026-11-14).
    assert same_bytes(
        out / "reference_selection_sequences.tsv", AP / "reference_selection_sequences.tsv"
    ), "a gene's chosen reference protein changed (reference_selection_sequences.tsv)"
    return (
        "reference_selection.tsv and reference_disagreements.tsv reproduce byte for byte; "
        "every gene's chosen protein matches reference_selection_sequences.tsv"
    )


def verify_isoform_gate(work: Path) -> str:
    run([str(HERE / "gate_isoform_cohort.py"), "--check"], HERE)
    spec = json.loads((AP / "panel_inputs.json").read_text())["isoform_gate"]
    assert sha256(HERE / "isoform_gate.tsv") == spec["sha256"], "gate table differs from its pin"
    from deliverome_analysis.s3 import processed

    pinned = Path(processed.download(spec["key"], work / "isoform_gate.tsv"))
    assert sha256(pinned) == spec["sha256"], "pinned gate key differs from panel_inputs.json"
    s = json.loads((HERE / "isoform_gate_summary.json").read_text())
    return (
        f"isoform_gate.tsv reproduces from pinned inputs and matches its pinned key: "
        f"{s['pass']} of {s['forms']} forms pass, {s['fail']} fail "
        f"({s['unverifiable_only']} only for want of a DeepTMHMM v1 row)"
    )


def verify_panels(work: Path) -> str:
    rej = HERE / "panel_rejections.tsv"
    spec = json.loads((AP / "panel_inputs.json").read_text())
    # The isoform redraw runs through the cohort gate (2026-11-12); the original ortholog
    # panel keeps its recorded redraw (retired from L1, kept as the record).
    extra = {"isoform": ["--gate", str(HERE / "isoform_gate.tsv")], "ortholog": []}
    for kind in ("isoform", "ortholog"):
        run([str(AP / f"sample_{kind}_panel.py"), "--snapshot", "--original"], work)
        assert sha256(work / f"{kind}_panel.tsv") == spec["reproduces"][f"{kind}_panel.tsv"], (
            f"{kind}: original draw does not reproduce"
        )
        run(
            [str(AP / f"sample_{kind}_panel.py"), "--snapshot", "--reject", str(rej), *extra[kind]],
            work,
        )
        assert same_bytes(work / f"{kind}_panel.tsv", HERE / f"{kind}_panel.tsv"), (
            f"{kind}: redraw does not reproduce the committed panel"
        )
    return (
        "both panels: original draw and rule-based redraw reproduce byte for byte (isoform "
        "redraw through the cohort gate)"
    )


def run_candidate_cds(work: Path, reference: Path | None = None) -> dict:
    """Resolve every candidate's CDS from the panel snapshot and a reference table (default:
    the published record)."""
    snap = json.loads((AP / "panel_inputs.json").read_text())
    from deliverome_analysis.s3 import processed

    for name in ("isoform_cohorts.tsv", "ortholog_repick_mouse.tsv", "ortholog_repick_cyno.tsv"):
        processed.download(snap["prefix"] + name, work / name)
    if reference is None:
        record(REFERENCE_RECORD + "reference_selection.tsv", work / "reference_selection.tsv")
    else:
        shutil.copy(reference, work / "reference_selection.tsv")
    run([str(AP / "resolve_candidate_cds.py")], work)
    return json.loads((work / "candidate_cds.json").read_text())


def verify_candidate_cds(work: Path) -> str:
    import apply_reference_fixes as fixes

    from deliverome_analysis import library_design as ld
    from deliverome_analysis.s3 import processed

    # Rerun against the current reference record (named references of 2026-11-12).
    new = run_candidate_cds(work)
    # The resolver's own record is the pre-adoption one; the current record extends it.
    rec = json.loads(Path(processed.download(ld.BASE_CDS_KEY, work / "rec.json")).read_text())
    named = {f"lib:{g}" for g in fixes.NAMED}
    # References moved by the GTF MANE rerun (apply_mane_gtf.py, 2026-11-22): like the named
    # ones, each must equal the rerun's entry.
    mane_gtf = HERE / "mane_gtf_changes_2026-11-22.tsv"
    if mane_gtf.exists():
        import csv

        with mane_gtf.open() as fh:
            named |= {r["key"] for r in csv.DictReader(fh, delimiter="\t")}
    differ = sorted(k for k in rec if k not in new or new[k]["cds"] != rec[k]["cds"])
    # (a gene added after the resolver's record is checked with the other additions below)
    assert set(differ) == named & set(rec), (
        f"designed CDS does not reproduce for {sorted(set(differ) ^ (named & set(rec)))}"
    )
    cur = json.loads(Path(processed.download(ld.CDS_KEY, work / "cur.json")).read_text())
    starts = {k: d for k, d in fixes.approved_starts().items() if k in cur and k.startswith("lib:")}
    trims = {k: s for k, s in fixes.TRIMS.items() if k in cur}
    # Cyno 5' restarts (apply_ortholog_start_trims.py) of orthologs already in the base record;
    # that stage checks each against its recorded CDS and protein hash.
    import apply_ortholog_start_trims as ortho_trims

    restarts = {
        r["key"]: r for r in ortho_trims.ledger() if r["status"] == "applied" and r["key"] in rec
    }
    assert all(
        ortho_trims.sha(cur[k]["cds"]) == r["restarted_cds_sha256"] for k, r in restarts.items()
    ), "an ortholog restart is not the recorded one"
    moved = sorted(k for k in rec if cur.get(k) != rec[k])
    allowed = (named & set(rec)) | set(starts) | set(trims) | set(restarts)
    assert set(moved) == allowed, f"unexpected changed entries {sorted(set(moved) ^ allowed)}"
    assert all(cur[k] == new[k] for k in named), "a named reference's CDS is not the rerun's"
    assert all(cur[k] == fixes.start_replaced(rec[k], d) for k, d in starts.items()), (
        "an ATG start replacement is not the approved one"
    )
    assert all(cur[k] == fixes.trimmed(rec[k], s) for k, s in trims.items()), (
        "an N-terminal trim is not the approved one"
    )
    appended = sorted(set(cur) - set(rec))
    # An ORF without a stop reads through into the vector; Ensembl omits it from some
    # records, which the resolver now repairs (stop_added). None may slip through.
    no_stop = sorted(k for k, v in cur.items() if v["cds"][-3:] not in ("TAA", "TAG", "TGA"))
    assert not no_stop, f"{len(no_stop)} designed CDS end without a stop codon: {no_stop[:5]}"
    # Ortholog-policy forms (2026-11-10) and the dated gene additions (gene_list_additions.tsv,
    # appended by add_library_genes.py); an addition must equal the rerun's CDS for it.
    import csv

    with ADDITIONS.open() as fh:
        added_genes = {f"lib:{r['gene']}" for r in csv.DictReader(fh, delimiter="\t")}
    # Root genes built from a recorded working copy (allele_restorations.py, 2026-11-22): no
    # Ensembl 115 transcript encodes them, so each must equal its ledger-checked CDS.
    restored = json.loads((HERE / "allele_restorations_cds.json").read_text())
    assert all(k.startswith("ort:") or k in added_genes or k in restored for k in appended), (
        "only ortholog forms, dated additions and allele restorations were appended"
    )
    assert all(cur[k]["cds"] == restored[k]["cds"] for k in appended if k in restored), (
        "an allele restoration's CDS is not the ledger's"
    )
    # A withdrawn addition (SCNN1B) keeps its record entry but is no longer on the gene list,
    # so the rerun no longer produces it.
    active_added = {f"lib:{g}" for g in additions()}
    assert all(cur[k] == new[k] for k in appended if k in active_added), (
        "an added gene's CDS is not the rerun's"
    )
    extra = sorted(set(new) - set(rec))
    return (
        f"all {len(rec)} resolver-recorded designed CDS reproduce except the {len(named)} named "
        f"references, which equal the current record; {len(starts)} approved ATG start "
        f"replacements; {len(trims)} in-frame N-terminal trim; the current record adds {len(appended)} ortholog-policy forms "
        f"(verified in ortholog_adoption), {sum(k in added_genes for k in appended)} dated "
        f"gene additions and {sum(k in restored for k in appended)} allele restorations, and "
        "every entry ends in a stop codon "
        f"({sum('stop_added' in v for v in cur.values())} appended stops); extra in rerun: "
        f"{extra or 'none'} (iso:MICA:Q29983 is MICA's selected reference, so the dedupe stage "
        "retires it)"
    )


def verify_reference_fixes(work: Path) -> str:
    import apply_reference_fixes as fixes

    return fixes.verify(work)


def verify_ortholog_start_trims(work: Path) -> str:
    import apply_ortholog_start_trims as trims

    return trims.verify(work)


def verify_constructs(work: Path) -> str:
    """Re-domesticate every designed CDS and re-assemble every vector; compare with the
    records (about a minute)."""
    import csv

    sys.path.insert(0, str(AP))
    import assemble_candidates
    import candidate_inputs
    import domesticate_candidates
    from Bio import SeqIO

    from deliverome_analysis import library_design as ld
    from deliverome_analysis import seqtools
    from deliverome_analysis.s3 import processed

    cds = json.loads(Path(processed.download(ld.CDS_KEY, work / "c.json")).read_text())
    con = json.loads(Path(processed.download(ld.CONSTRUCTS_KEY, work / "k.json")).read_text())
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


def verify_cterm_review(work: Path) -> str:
    out = work / "l1_cterm_review.tsv"
    run([str(HERE / "review_l1_cterm_members.py"), "--out", str(out)], HERE)
    assert same_bytes(out, HERE / "l1_cterm_review.tsv"), "l1_cterm_review.tsv differs"
    run([str(HERE / "l1_cterm_decisions.py")], work)  # writes next to the script
    proc = subprocess.run(["git", "diff", "--quiet", "--", "l1_cterm_decisions.tsv"], cwd=HERE)
    assert proc.returncode == 0, "l1_cterm_decisions.tsv does not reproduce from the review"
    return "l1_cterm_review.tsv reproduces byte for byte from the pinned UniProt snapshot"


def verify_cterm_arms(work: Path) -> str:
    run([str(HERE / "diversify_snorkel_modules.py"), "--check"], HERE)
    run([str(HERE / "build_l1_cterm_arms.py"), "--check"], HERE)
    return (
        "every C-terminal-rung arm rebuilds byte for byte: protein = native[:n] + cassette + "
        "native[n:] (n = the last residue, or the decided site for ITM2B), one stop, "
        "no BbsI/Esp3I/BsaI beyond the wrapper (also in p0039), swap tract re-measured"
    )


#: The four tag-rung builders whose arm tables the synthesis stage reads (the C-terminal one is
#: its own stage above). Each rebuilds its table from its decision table and checks it.
RUNG_BUILDERS = (
    "build_l1_nterm_arms.py",
    "build_l1_multipass_arms.py",
    "build_l1_type3_arms.py",
    "build_l1_internal_arms.py",
)


def verify_rung_arms(work: Path) -> str:
    scripts = ("diversify_arm_cassettes.py", *RUNG_BUILDERS)
    with ThreadPoolExecutor(max_workers=len(scripts)) as pool:
        list(pool.map(lambda s: run([str(HERE / s), "--check"], HERE), scripts))
    return (
        "every upstream arm's linker encoding replays (distinct per arm, p003 tag codons), and "
        "the N-after-SP, multipass, type III and internal arm tables rebuild from their decision "
        "tables and match the committed records"
    )


def verify_consistency(work: Path) -> str:
    run([str(HERE / "l1_consistency.py"), "--check"], HERE)
    run([str(HERE / "l1_completeness.py"), "--check"], HERE)
    return (
        "every built arm matches its protein, layout and fold; every failing refold has a "
        "recorded disposition; every active member is placed on its current ids"
    )


def verify_saved_inputs(work: Path) -> str:
    import pipeline_input_manifest

    return pipeline_input_manifest.check()


def verify_pinned_refs(work: Path) -> str:
    import pinned_refs

    return pinned_refs.check()


def verify_structure_models(work: Path) -> str:
    import pin_structure_models

    return pin_structure_models.check()


def verify_rung_builds(work: Path) -> str:
    import save_rung_builds

    return save_rung_builds.check()


def verify_external_inputs(work: Path) -> str:
    import external_inputs

    return external_inputs.check()


def verify_synthesis_compliance(work: Path) -> str:
    run([str(HERE / "synthesis_compliance.py"), "--check"], HERE)
    return (
        "every tagged insert at synthesis_compliance.INPUT_REF re-checks and recodes to the "
        "committed tables: synonymous codons only, protein, cassette DNA, wrapper and readout "
        "block unchanged, no new Type IIS site, splice donor or Tn5-ME match"
    )


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


def verify_ortholog_adoption(work: Path) -> str:
    """The L1 ortholog panel: the final-v2 draw, the correspondence swaps, the members."""
    import adopt_ortholog_cds_policy as adopt
    import review_ortholog_correspondence as review

    rows, _, _ = review.replay_draw(work / "draw")  # asserts the committed panel back
    detail = adopt.verify(work / "adopt")
    return f"final-v2 draw over {len(rows)} sampleable forms reproduces the panel; {detail}"


def verify_tranches(work: Path) -> str:
    run([str(HERE / "dedupe_library_members.py"), "--check"], HERE)
    return "no two active members share a designed protein"


def verify_topology(work: Path) -> str:
    from deliverome_analysis import topology
    from deliverome_analysis.s3 import external

    snap = Path(external.download(topology.SNAPSHOT_KEY, work / "topology.parquet"))
    assert sha256(snap) == topology.SNAPSHOT_SHA256, "topology snapshot differs from its pin"
    fasta = work / "gap.fasta"
    run([str(HERE / "publish_topology_supplement.py"), "fasta", "--out", str(fasta)], HERE)
    assert fasta.read_text() == "", "active members without a DeepTMHMM v1 topology remain"
    return (
        f"every active L1 member has a v1 topology in the pinned snapshot "
        f"({topology.SNAPSHOT_KEY.split('/')[-2]}, SHA-256 checked)"
    )


def verify_tag_paths(work: Path) -> str:
    for lib in ("l1", "l2"):
        out = work / f"{lib}_tag_paths.tsv"
        run([str(HERE / "l1_tag_paths.py"), "--library", lib, "--out", str(out)], HERE)
        assert same_bytes(out, HERE / f"{lib}_tag_paths.tsv"), f"{lib}_tag_paths.tsv differs"
    return "l1_tag_paths.tsv and l2_tag_paths.tsv reproduce byte for byte"


def verify_renumber(work: Path) -> str:
    run([str(HERE / "renumber_library_members.py"), "--check"], HERE)
    summary = json.loads((HERE / "l1_member_renumbering_summary.json").read_text())
    return (
        f"{summary['inserts']} inserts of {summary['members_built']} built L1 members number "
        f"to the committed l1-i ids on {summary['plates']} plates; members are keyed by design "
        "id and every re-keyed record is current"
    )


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
        "library/l1/surface-filters/final_library.py (offline)",
        "release 1.3.0 deep dives + preset table (Zenodo 10.5281/zenodo.23445521, on "
        "external://), expression TPM tables, pinned MyGene responses, curated S1E3 calls",
        "verified",
        verify=verify_gene_list,
    ),
    Stage(
        "gpi_set",
        "which genes are GPI-anchored, per form (the GPI set of record)",
        "library/l1/orf-design/build_gpi_set.py",
        "NetGPI 1.1 predictions (committed) + curated gpi_anchored_genes.tsv",
        "verified",
        verify=verify_gpi_set,
    ),
    Stage(
        "reference",
        "each gene's reference proteoform: MANE+APPRIS agree -> the single "
        "source -> observed evidence (PeptideAtlas build 607, GENCODE Ribo-seq ORFs) -> MANE; "
        "named overrides in reference_overrides.tsv",
        "library/l1/proteoform-selection/select_reference.py --offline",
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
        "each gene's mouse and cyno ortholog proteoform in the April (Compara 115) pool: the "
        "source of the L2 ortholog members and of the L1 ortholog panel retired on 2026-11-10",
        "library/l1/proteoform-selection/repick_orthologs.py",
        "rebuilt universe table, NetGPI ortholog calls, Ensembl 115 + UniProt 2026_03 "
        "ortholog files, surfaceome f9c73691 cohorts (repick_orthologs.py @c4dbe752)",
        "verified",
        "",
        verify_ortholog_repick,
    ),
    Stage(
        "universe_coverage",
        "recorded ortholog coverage of that April pool relative to the library gene list, and "
        "the selected isoform cohort relative to the current library gene list; "
        "isoform cohort absence is not universe absence",
        "library/l1/orf-design/check_universe_coverage.py",
        "Ensembl Compara 115 homology dumps (external://), accessible-surfaceome f9c73691 "
        "query set and ortholog lists, Ensembl 115 peptides, UniProt 2026_03 isoforms",
        "verified",
        verify=verify_universe_coverage,
    ),
    Stage(
        "isoform_gate",
        "which isoform-cohort forms may be drawn: native ATG start supported by an Ensembl "
        "115 CDS, and the reference's signal peptide, N-terminal orientation and first TM "
        "helix kept (DeepTMHMM v1); Becca, 2026-11-12",
        "gate_isoform_cohort.py",
        "panel snapshot cohort, Ensembl 115 pep/cds, UniProt 2026-10-11, L1 reference "
        "designed proteins, pinned topology snapshot",
        "verified",
        verify=verify_isoform_gate,
    ),
    Stage(
        "panels",
        "the 150-isoform panel (L1; redrawn through the cohort gate on 2026-11-12) and the "
        "original 200-ortholog panel (retired from L1 on 2026-11-10, kept as the record): "
        "seeded stratified draws, rejected picks replaced by continuing the draw",
        "sample_isoform_panel.py / sample_ortholog_panel.py --snapshot; redraw_panels.py",
        "panel_inputs.json snapshot (13 files) + panel_rejections.tsv + isoform_gate.tsv "
        "(pinned key)",
        "verified",
        verify=verify_panels,
    ),
    Stage(
        "ortholog_adoption",
        "the L1 ortholog panel: the Compara 116 CDS-eligible draw (final-v2, seeds 20261011/"
        "20261012), with the 11 correspondence-review swaps and 6 asks kept "
        "(ortholog_panel_adoption_judgements.tsv); members kept, moved from L2 or added",
        "draw_eligible_ortholog_panel.py (replayed by review_ortholog_correspondence.py); "
        "adopt_ortholog_cds_policy.py",
        "final-v2 archive (processed://.../ortholog-cds-eligibility/2026-11-05/final-v2/), "
        "ortholog-rebuild-inputs/2026-11-05 manifest + CDS, committed review tables",
        "verified",
        verify=verify_ortholog_adoption,
    ),
    Stage(
        "rescue_arms",
        "25 expression-rescue arms, their CDS and constructs (retired 2026-11-13, one design per gene; "
        "kept and verified as the record)",
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
        "library/l1/proteoform-selection/resolve_candidate_cds.py (candidate_inputs.py)",
        "Ensembl 115 CDS/pep, UniProt 2026-10-11 (S3); reference record; panel snapshot; "
        "the 56 ortholog-policy forms appended from final-v2 (stage ortholog_adoption)",
        "verified",
        verify=verify_candidate_cds,
    ),
    Stage(
        "reference_fixes",
        "the approved reference-form fixes: eight named references (CD200R1L, LRRC37B, "
        "CACHD1, ADGRF3, RGMB; LMLN, DISP2, TSPAN10), three ATG start replacements "
        "(SLC22A17, TRPV6, OR8G3) and EPHA6's in-frame trim to Q9UF33",
        "apply_reference_fixes.py (ledger reference_changes.tsv)",
        "reference_overrides.tsv (pipeline-inputs/2026-11-12-round2b), "
        "proteoform_decisions_applied.tsv, "
        "library_cds_boundary_decisions.tsv, the 2026-11-10 CDS/construct records",
        "verified",
        verify=verify_reference_fixes,
    ),
    Stage(
        "ortholog_start_trims",
        "cyno L1 orthologs restarted at the native Met aligned to human Met1 (5' trim, no codon "
        "edited; a recorded exception to the ortholog CDS policy), applied once the restarted "
        "protein has a DeepTMHMM v1 row: all four applied (CHRM3, NIPAL4, OR3A2 after their "
        "v1 rows were published on 2026-11-12)",
        "apply_ortholog_start_trims.py (ledger ortholog_start_trims.tsv)",
        "the reference-fixes-2 CDS/construct records, the pinned topology snapshot",
        "verified",
        verify=verify_ortholog_start_trims,
    ),
    Stage(
        "constructs",
        "domesticated ORFs (BbsI/Esp3I/BsaI) and p0039 vectors",
        "domesticate_candidates.py, assemble_candidates.py (seed 20261011)",
        "library_design.CDS_KEY / CONSTRUCTS_KEY records, p0039 backbone (S3)",
        "verified",
        verify=verify_constructs,
    ),
    Stage(
        "tranches",
        "L1/L2 membership; duplicates retired on the designed protein",
        "split_library_tranches.py (once), dedupe_library_members.py, redraw_panels.py, "
        "adopt_ortholog_cds_policy.py",
        "candidate_manifest_tranches.tsv (committed record)",
        "verified",
        verify=verify_tranches,
    ),
    Stage(
        "topology",
        "DeepTMHMM v1 topology for every active member's designed protein",
        "publish_topology_supplement.py + modal deeptmhmm_v1_app.py",
        "external://accessible-surfaceome/topology/2026-11-12-cyno-restarts/"
        "topology.parquet (topology.SNAPSHOT_KEY, SHA-256 pinned)",
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
        "l1_swap_risk.py (seed 20261029); l1_cassette_encodings.tsv records each batch, "
        "including retired members' encodings, which stay issued",
        "member CDS, tag paths, p003 cassette codons",
        "verified",
        verify=verify_cassettes,
    ),
    Stage(
        "renumber",
        "each active L1 member's id (l1-mNNNN) in plate order, by (insert_nt, design id), 96 "
        "per plate column-major; the design id (l1-pNNNN) is kept; current L1 records re-keyed",
        "renumber_library_members.py; l1_member_renumbering.tsv records the mapping",
        "candidate_manifest_tranches.tsv, designed CDS (candidate_cds.json, rescue constructs)",
        "verified",
        verify=verify_renumber,
    ),
    Stage(
        "cterm_review",
        "which C-terminal members keep a C-terminal tag on the cell: held when a C-terminal "
        "peptide is cut off (UniProt) or the gene is GPI; every other risk recorded as a flag",
        "review_l1_cterm_members.py",
        "tag paths, v1/v2 topology snapshot, UniProt entries (snapshot of 2026-11-10 on S3), "
        "member structure models (AlphaFold DB v6 / ColabFold on S3), gpi_set.tsv",
        "verified",
        verify=verify_cterm_review,
    ),
    Stage(
        "cterm_arms",
        "each C-terminal member's insert: domesticated ORF + its issued cassette + TAA, in "
        "the p003 layout and Esp3I wrapper (barcode slot a placeholder until issue #1)",
        "build_l1_cterm_arms.py",
        "candidate_constructs.json, l1_cassette_encodings.tsv, l1_snorkel_module_encodings.tsv "
        "(each snorkel arm's own module DNA, diversify_snorkel_modules.py), l1_cterm_review.tsv, "
        "p0039",
        "verified",
        verify=verify_cterm_arms,
    ),
    Stage(
        "rung_arms",
        "the tag cassette in every N-after-SP, residue-1 multipass, type III and internal arm, "
        "at the decided site (the decisions are the rung sessions' tables), each on its own "
        "linker encoding with p003's tag codons (Becca, 2026-11-21: no shared linkers)",
        "diversify_arm_cassettes.py (seed 20261121; l1_arm_cassette_encodings.tsv), "
        + ", ".join(RUNG_BUILDERS),
        "reviewer-package/nterm/nterm_decisions.tsv, nterm_res1_multipass.tsv, "
        "type3_nterm_calls.tsv, S5E13 l1_internal_decisions.tsv; candidate_constructs.json",
        "verified",
        verify=verify_rung_arms,
    ),
    Stage(
        "consistency",
        "the built arms are internally consistent and complete: protein, cassette and layout per "
        "arm, the fold of the exact protein, a disposition for every refold failure",
        "l1_consistency.py, l1_completeness.py",
        "the five constructs_l1_*.tsv/.fasta, tag_refold_structures.tsv, l1_member_renumbering.tsv",
        "verified",
        verify=verify_consistency,
    ),
    Stage(
        "synthesis_compliance",
        "vendor synthesis check on every tagged insert (rules calibrated on Twist's report of "
        "the first draft; an approximation) and the minimal synonymous recodes that clear it",
        "synthesis_compliance.py",
        "each rung's constructs_l1_<rung>.fasta/.tsv at INPUT_REF; l1_tag_paths.tsv",
        "verified",
        verify=verify_synthesis_compliance,
    ),
    Stage(
        "saved_inputs",
        "every S3 object the stages read is still saved, unchanged: the 99 objects a clean-clone "
        "verify with an empty cache downloaded, each with its SHA-256",
        "pipeline_input_manifest.py",
        "library_pipeline_inputs_saved.tsv (+ .json: tool versions, lockfiles, git refs)",
        "verified",
        verify=verify_saved_inputs,
    ),
    Stage(
        "pinned_refs",
        "every git commit a stage replays at (synthesis INPUT_REF, the fold check's refold "
        "table, ...) is an ancestor of HEAD, so main alone rebuilds the order",
        "pinned_refs.py",
        "the *_REF / *_COMMIT constants in the s4e9 scripts",
        "verified",
        verify=verify_pinned_refs,
    ),
    Stage(
        "structure_models",
        "every AlphaFold DB model the design read (structure payload, refold natives) is in the "
        "write-once external://alphafold-db/v6/ archive with its SHA-256",
        "pin_structure_models.py",
        "tag_refold_afdb_archive.tsv (l1-tag-refold's archive_afdb_models.py)",
        "verified",
        verify=verify_structure_models,
    ),
    Stage(
        "rung_builds",
        "the five rung builds every order of record was made from are saved on S3, so no order "
        "depends on a git branch; the final order's inputs are the saved ones",
        "save_rung_builds.py",
        "rung_builds_saved.json",
        "verified",
        verify=verify_rung_builds,
    ),
    Stage(
        "external_inputs",
        "every external source a script reads (UniProt, Ensembl, MyGene, AlphaFold DB, "
        "TmAlphaFold, NetGPI, D1, licensed predictors, Modal images) has a saved copy or a "
        "recorded tool identity, and every saved copy is unchanged",
        "external_inputs.py",
        "external_inputs.tsv",
        "verified",
        verify=verify_external_inputs,
    ),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    v = sub.add_parser("verify")
    v.add_argument("--stage", choices=[s.id for s in STAGES])
    v.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 2),
        help="max stages to verify in parallel (default: min(4, cpu_count))",
    )
    args = ap.parse_args(argv)

    if args.cmd == "plan":
        for n, s in enumerate(STAGES, 1):
            print(
                f"{n:>2}. {s.id:<18} [{s.status}] {s.decides}\n      {s.script}\n"
                f"      inputs: {s.inputs}"
                + (f"\n      pinned: {s.why_pinned}" if s.why_pinned else "")
            )
        return 0

    selected = [s for s in STAGES if not args.stage or s.id == args.stage]

    def _run_one(s: Stage) -> tuple[Stage, str, str]:
        if s.verify is None:
            return s, "PINNED", s.why_pinned
        work = Path(tempfile.mkdtemp(prefix=f"pipeline-{s.id}-"))
        try:
            return s, "OK", s.verify(work)
        except (AssertionError, SystemExit) as exc:
            return s, "FAIL", str(exc)

    jobs = max(1, min(args.jobs, len(selected)))
    if jobs == 1:
        outcomes = [_run_one(s) for s in selected]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            outcomes = list(pool.map(_run_one, selected))

    failed = 0
    for s, status, msg in outcomes:
        if status == "PINNED":
            print(f"-- {s.id:<18} pinned: {msg}")
        elif status == "OK":
            print(f"OK {s.id:<18} {msg}")
        else:
            failed += 1
            print(f"FAIL {s.id:<16} {msg}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
