"""Every copy of an L1 definition agrees with its one home in the shared package.

The insert wrapper, the cassette, the tag peptides, the p003 tag codons, the Igk leader and the
insert/FASTA readers were each written out in several scripts (CODE_REVIEW_2026-11-22.md). Until
each copy is replaced by an import, these tests pin it to the definition in
``deliverome_analysis.plasmid_qc`` / ``seqtools``, so a copy that changes on its own fails here
instead of building or checking a different library. They read only committed files and import
only scripts with no work at import time.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import random
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
S5E13 = REPO / "library/l1/tag-site-benchmarking"
AP = REPO / "library/l1/proteoform-selection"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(AP))

from deliverome_analysis.cloning import plasmid_qc as qc  # noqa: E402
from deliverome_analysis.cloning import seqtools  # noqa: E402

#: The six rung tables; every built L1 insert is in exactly one.
RUNG_FASTAS = sorted(HERE.glob("constructs_l1_*.fasta"))
RUNG_FASTAS = [p for p in RUNG_FASTAS if not p.name.endswith("_proteins.fasta")]


def _load(path: Path) -> ModuleType:
    """Import a script by path under a unique name (two folders both have e.g. refold_queue)."""
    name = "l1defs_" + path.relative_to(REPO).as_posix().replace("/", "_").replace("-", "_")
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        del sys.modules[name]
        raise
    return mod


def _functions(path: Path, *names: str) -> SimpleNamespace:
    """Named top-level functions and constants of a script, executed from its source alone.

    For scripts whose own imports resolve by folder name (``tag_refold/`` and s5e13 both have a
    ``refold_queue``), so importing them depends on which test ran first. Only the requested
    definitions run, in file order, with the standard modules they use.
    """
    nodes: list[ast.stmt] = []
    found: set[str] = set()
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            nodes.append(node)
            found.add(node.name)
        elif isinstance(node, ast.Assign):
            hit = {t.id for t in node.targets if isinstance(t, ast.Name) and t.id in names}
            if hit:
                nodes.append(node)
                found |= hit
    assert found == set(names), f"{path.name} lacks {set(names) - found}"
    namespace: dict[str, object] = {
        "Path": Path,
        "hashlib": hashlib,
        "qc": qc,
        "seqtools": seqtools,
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    return SimpleNamespace(**{n: namespace[n] for n in names})


def _literal(path: Path, name: str) -> object:
    """A module-level constant read from source without importing the script."""
    for node in ast.parse(path.read_text()).body:
        targets = node.targets if isinstance(node, ast.Assign) else []
        if isinstance(node, ast.AnnAssign) and node.value is not None:
            targets = [node.target]
        if any(isinstance(t, ast.Name) and t.id == name for t in targets):
            value = node.value
            assert value is not None
            for sub in ast.walk(value):
                assert isinstance(sub, ast.Constant | ast.BinOp | ast.Add | ast.Tuple | ast.Load), (
                    f"{path.name}:{name} is not a plain literal"
                )
            return eval(compile(ast.Expression(value), str(path), "eval"), {"__builtins__": {}})
    raise AssertionError(f"{path.name} defines no {name}")


@pytest.fixture(scope="module")
def inserts() -> dict[str, str]:
    out: dict[str, str] = {}
    for path in RUNG_FASTAS:
        for key, seq in seqtools.fasta_by_id(path).items():
            assert key not in out, f"{key} in two rung tables"
            out[key] = seq
    return out


#: Rung tables not rebuilt with the L1 readout region on 2026-11-22: they keep p003's layout.
P003_LAYOUT_TABLES = {"constructs_l1_rescue.fasta"}


# --- 1. The insert wrapper ---------------------------------------------------------------- #


def test_wrapper_constants_agree_across_the_three_modules_that_define_them():
    import assemble_candidates as ac
    import build_constructs as bc

    for name in ("KOZAK", "STOP", "SPACER", "FWD_PRIMER", "REV_PRIMER", "CS1"):
        assert getattr(bc, name) == getattr(qc, name), f"build_constructs.{name}"
    for name in ("KOZAK", "REV_PRIMER", "CS1", "ESP3I_5", "ESP3I_3"):
        assert getattr(ac, name) == getattr(qc, name), f"assemble_candidates.{name}"
    assert (ac.SLOT_START, ac.SLOT_LEN) == (qc.BACKBONE_SLOT_START, qc.BACKBONE_SLOT_LEN)
    assert tuple(ac.ENZYMES) == tuple(bc.ENZYMES) == qc.DOMESTICATION_ENZYMES


def test_every_builder_and_the_pool_qc_use_the_same_wrapper():
    import l1_pool_qc as pool

    for name in ("nterm", "cterm", "internal", "multipass", "type3"):
        mod = _load(HERE / f"build_l1_{name}_arms.py")
        assert mod.BARCODE_SLOT == qc.BARCODE_SLOT, name
        assert tuple(mod.ENZYMES) == qc.DOMESTICATION_ENZYMES, name
    assert pool.BARCODE_SLOT == qc.BARCODE_SLOT
    assert pool.OPENING == qc.ESP3I_5 + qc.KOZAK + "ATG"
    assert pool.READOUT_BLOCK == qc.wrap_insert("")[qc.CODING_START :]
    assert _load(HERE / "build_l1_cterm_arms.py").TAIL == qc.STOP_TO_BARCODE
    assert _load(HERE / "l1_swap_risk.py").TAIL == qc.STOP_TO_BARCODE


def test_every_built_insert_is_its_coding_sequence_in_the_wrapper(inserts):
    assert len(inserts) > 3000
    p003 = {k for p in RUNG_FASTAS if p.name in P003_LAYOUT_TABLES for k in seqtools.fasta_by_id(p)}
    assert p003, "the rescue arms keep p003's layout"
    for key, insert in inserts.items():
        wrap = qc.wrap_p003_insert if key in p003 else qc.wrap_insert
        tail = len(wrap("")) - qc.CODING_START
        coding = insert[qc.CODING_START : len(insert) - tail]
        assert wrap(coding) == insert, key
        assert coding.startswith("ATG") and seqtools.translate(coding[-3:]) == "*", key
        assert qc.wrapper_sites_only(insert), key


def test_the_pool_qc_places_an_insert_where_the_shared_function_does():
    import l1_pool_qc as pool

    rng = random.Random(1)
    backbone = "".join(rng.choice("ACGT") for _ in range(qc.BACKBONE_SLOT_START + 400))
    insert = qc.wrap_insert("ATG" + "GCC" * 10 + "TAA")
    assert pool.assemble(backbone, insert) == qc.place_in_backbone(backbone, insert)


# --- 2 and 3. Cassette, tag peptides, p003 tag codons ------------------------------------- #


def test_every_cassette_peptide_copy_is_the_canonical_one():
    import build_constructs as bc

    assert bc.CASSETTE_AA == qc.P003_CASSETTE_PEPTIDE
    assert bc.CASSETTE == qc.P003_CASSETTE_NT
    assert _load(HERE / "tag_refold/build_refold_inputs.py").CASSETTE == qc.P003_CASSETTE_PEPTIDE
    assert _load(S5E13 / "build_internal_tagged.py").CASSETTE_AA == qc.P003_CASSETTE_PEPTIDE
    assert _load(S5E13 / "l1_internal_additions.py").CASSETTE_AA == qc.P003_CASSETTE_PEPTIDE


def test_every_tag_peptide_copy_is_the_canonical_one():
    import diversify_cassettes as dc

    assert tuple(dc.TAG_SEGMENTS) == tuple(qc.TAG_PEPTIDES.values())
    refold = _load(HERE / "tag_refold/build_refold_inputs.py")
    assert refold.SPYTAG003 == qc.TAG_PEPTIDES["SpyTag003"]
    # build_refold_inputs keeps ALFA without its flanking prolines.
    assert "P" + refold.ALFA + "P" == qc.TAG_PEPTIDES["P-ALFA-P"]


def test_the_p003_tag_codons_are_the_same_wherever_they_are_computed():
    import build_constructs as bc
    import diversify_cassettes as dc

    assert dc.tag_codons() == qc.p003_tag_codons()
    assert seqtools.pin_codons(bc.CASSETTE, bc.CASSETTE_AA, dc.TAG_SEGMENTS) == qc.p003_tag_codons()


# --- 8. The Igk leader -------------------------------------------------------------------- #


def test_every_igk_leader_copy_is_build_gpcr_rescue_arms_leader():
    import build_gpcr_rescue_arms as g

    assert seqtools.translate(g.LEADER) == g.LEADER_AA
    assert seqtools.translate(g.PLUS1) == g.PLUS1_AA
    # synthesis_compliance protects the leader and the +1 Asp from recoding.
    assert _load(HERE / "synthesis_compliance.py").IGK_LEADER == g.LEADER_AA + g.PLUS1_AA
    # verify_ordered_constructs and rung_numbers strip the leader alone to reach the mature chain.
    assert _load(HERE / "tag_refold/verify_ordered_constructs.py").IGK_LEADER == g.LEADER_AA
    assert _load(HERE / "reviewer-package/rung_numbers.py").IGK_MATURE_LEADER == g.LEADER_AA


# --- 4. Reading an insert's protein ------------------------------------------------------- #


def test_every_insert_translator_agrees_with_the_shared_one(inserts):
    import l1_consistency

    refold = _functions(HERE / "tag_refold/build_refold_inputs.py", "CODON", "translate_insert")
    fold = _functions(
        HERE / "reviewer-package/fold_order_check.py", "BASES", "AMINO", "CODON", "translate"
    )
    for key, insert in inserts.items():
        protein = qc.translate_insert(insert)
        assert l1_consistency.translate_insert(insert) == protein, key
        assert refold.translate_insert(insert) == protein, key
        assert fold.translate(insert) == protein, key


# --- 5. FASTA readers --------------------------------------------------------------------- #


@pytest.mark.parametrize("path", RUNG_FASTAS, ids=lambda p: p.name)
def test_every_fasta_reader_agrees_with_the_shared_one_on_the_rung_tables(path):
    import diversify_arm_cassettes as dac
    import l1_consistency
    import l1_pool_qc as pool

    by_id = seqtools.fasta_by_id(path)
    pairs = seqtools.read_fasta(path)
    upper = [(h, s.upper()) for h, s in pairs]
    assert dac.read_fasta(path) == by_id
    assert l1_consistency.fasta(path) == by_id
    assert pool.read_fasta(path) == upper
    assert _functions(HERE / "tag_refold/build_refold_inputs.py", "fasta").fasta(path) == upper
    verify = _functions(HERE / "tag_refold/verify_ordered_constructs.py", "fasta")
    assert verify.fasta(path.read_text()) == {k: v.upper() for k, v in by_id.items()}


def test_build_constructs_reader_agrees_on_the_files_annotate_insert_reads():
    import annotate_insert
    import build_constructs as bc

    for name in annotate_insert.SOURCES:
        path = HERE / name
        # build_constructs keys on the header up to its first "|" only; these headers have no
        # space before it, so both keyings agree here (they would not on a multipass table).
        assert bc.read_fasta(path) == seqtools.fasta_by_id(path), name


# --- 9. The swap-tract measure ------------------------------------------------------------ #


def test_common_suffix_copies_agree():
    swap = _load(HERE / "l1_swap_risk.py")
    pool_model = _load(HERE / "l1_swap_pool_model.py")
    rng = random.Random(7)
    for _ in range(500):
        tail = "".join(rng.choice("AC") for _ in range(rng.randint(0, 30)))
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 20))) + tail
        b = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 20))) + tail
        expect = seqtools.common_suffix(a, b)
        assert swap.common_suffix(a, b) == expect
        assert pool_model.common_suffix(a, b) == expect


def test_plot_rescue_selection_swap_rate_matches_seqtools():
    plot = _load(HERE / "plot_rescue_selection.py")
    for bp in (0, 50, 100, 432, 762, 972):
        assert plot.swap_rate(bp) == pytest.approx(seqtools.barcode_swap_rate(bp))


def test_pool_qc_longest_shared_matches_seqtools():
    import l1_pool_qc as pool

    rng = random.Random(19)
    windows = ["".join(rng.choice("ACGT") for _ in range(80)) for _ in range(10)]
    windows[1] = windows[1][:20] + windows[0][10:42] + windows[1][52:]
    owners = [0, 1, 1, 2, 3, 4, 5, 6, 7, 8]
    assert pool.longest_shared(windows, owners) == seqtools.longest_shared_across_pool(
        windows, owners
    )


# --- 7. Pinned keys shared across scripts ------------------------------------------------- #


def test_order_of_record_fasta_key_agrees_across_synthesis_and_fold_check():
    syn = _load(HERE / "synthesis_compliance.py")
    fold_key = _literal(HERE / "reviewer-package/fold_order_check.py", "FASTA_KEY")
    assert fold_key == syn.FASTA_KEY


def test_compara116_archive_pins_agree_between_review_and_adoption():
    adopt = _load(HERE / "adopt_ortholog_cds_policy.py")
    review = _load(HERE / "review_ortholog_correspondence.py")
    pairs = [
        ("selected_cds", "selected_cds"),
        ("eligible_proteins", "draw_eligible"),
        ("provisional_picks", "draw_picks"),
        ("base_constructs", "draw_constructs"),
        ("pre_manifest", "draw_manifest"),
    ]
    for adopt_name, review_name in pairs:
        assert adopt.INPUTS[adopt_name] == review.INPUTS[review_name], (
            f"Compara-116 pin drifted: adopt.{adopt_name} vs review.{review_name}"
        )
    # base_cds / human_cds point to two keys with the same pinned SHA-256.
    assert adopt.INPUTS["base_cds"][2] == review.INPUTS["human_cds"][2]
