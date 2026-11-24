"""Re-pick each ortholog's proteoform by correspondence to the selected human reference.

The Compara one-to-one high-confidence mapping is gene level and stands. What is replaced
is the proteoform taken inside each confirmed ortholog gene, which was the **longest
symbol-matched entry, reviewed or not** -- a length tiebreak with no reviewed-entry floor,
leaving 654 of 1,971 mouse picks and 1,855 of 1,927 cyno picks on unreviewed TrEMBL.

An ortholog construct exists to be the counterpart of a human one, so the property that
matters is resemblance to the human form we actually selected, not which proteoform is
principal in the ortholog's own species. Selecting on correspondence also makes "the human
principal and the mouse principal are different splice forms" impossible rather than merely
detectable.

**Candidate pool.** UniProt proteoforms for the ortholog gene, the current cohort pick, and
**every Ensembl translation of the ortholog gene**. That last source matters: scored without
it, APPRIS beat correspondence on 58 mouse genes purely because its candidate was not in the
pool.

**Score.** Residues in exactly-matching contiguous blocks over the longer of the two
lengths. Dividing by the longer length is what stops a fragment winning -- a 200-residue
candidate against a 600-residue reference cannot exceed 0.33. Not Needleman-Wunsch identity;
see METHODS_library_design.md.

**Ties** break to a reviewed UniProt entry, then to the APPRIS principal where the species
has one (mouse does, cynomolgus does not), then to the current pick so the result is stable.

    uv run --locked python repick_orthologs.py --species mouse
    uv run --locked python repick_orthologs.py --species cyno
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from pathlib import Path

from gpi import gpi_genes as gpi_gene_set
from select_principal import (
    APPRIS_RELEASE,
    load_uniprot_proteoforms,
    parse_appris,
    parse_ensembl_pep,
    parse_uniprot_fasta,
    select_principal,
)
from select_reference import translations_by_gene

TIE = 0.001  # scores within this are a tie, broken on provenance rather than noise

# An ortholog reference must not be a soluble form, for the same reason a human one must
# not. The test is keyed on the **human** gene's GPI status, not the ortholog's: a
# GPI-anchored protein has a signal peptide and no transmembrane helix, so DeepTMHMM calls
# it SP, and rejecting SP outright would discard 57 legitimate GPI orthologs to catch 3 real
# soluble picks. Orthologs of GPI proteins are GPI, so the human keyword is the right source.
SURFACE_LABELS = {"TM", "SP+TM", "BETA"}
NON_SURFACE_LABELS = {"SP", "GLOB"}

SPECIES = {
    "mouse": {
        "cohort": "mouse_ortholog_one2one_highconf_non_hla",
        "pep": "Mus_musculus.GRCm39.pep.all.fa.gz",
        "uniprot": "uniprot_mouse_reviewed_isoforms.tsv.gz",
        "appris": f"appris_mus_musculus_{APPRIS_RELEASE}.txt",
    },
    "cyno": {
        "cohort": "cyno_ortholog_one2one_highconf_non_hla",
        "pep": "Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz",
        "uniprot": "uniprot_cyno_all_isoforms.tsv.gz",
        "appris": None,  # APPRIS covers M. mulatta, not M. fascicularis
    },
}


def score(ref: str, cand: str) -> float:
    sm = SequenceMatcher(None, ref, cand, autojunk=False)
    return sum(b.size for b in sm.get_matching_blocks()) / max(len(ref), len(cand))


def load_cyno_uniprot(path: Path) -> tuple[dict[str, list[dict]], set[str]]:
    """{gene symbol: [proteoform]} plus the set of reviewed accessions.

    Cynomolgus is 96% TrEMBL, so unlike mouse the unreviewed entries cannot be discarded --
    they are most of what exists. Reviewed status is kept as a tiebreak instead.
    """
    by_symbol: dict[str, list[dict]] = defaultdict(list)
    reviewed: set[str] = set()
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym = (row.get("Gene Names (primary)") or "").strip()
            acc = row["Entry"]
            if not sym:
                continue
            by_symbol[sym].append({"accession": acc, "seq": row["Sequence"]})
            if (row.get("Reviewed") or "").lower().startswith("reviewed"):
                reviewed.add(acc)
    return dict(by_symbol), reviewed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", choices=sorted(SPECIES), required=True)
    ap.add_argument("--surfaceome-repo", default="~/Git/accessible-surfaceome")
    ap.add_argument("--cache-dir", default="~/.cache/deliverome/appris")
    ap.add_argument("--universe", default="isoform_universe_scored.tsv")
    ap.add_argument(
        "--dump-picked", default="",
        help="write every pick's sequence to this FASTA. repick holds the candidate "
             "sequences already, including the unreviewed TrEMBL entries that are not in "
             "any cached reviewed table.",
    )
    ap.add_argument(
        "--netgpi", default="netgpi_ortho_candidate_predictions.txt",
        help="NetGPI calls for the GPI-gene candidates, from --dump-gpi-candidates",
    )
    ap.add_argument(
        "--dump-gpi-candidates", default="",
        help="write every candidate ortholog of a GPI gene to this FASTA and stop. The "
             "anchor guard needs a NetGPI call per candidate, not per picked form.",
    )
    args = ap.parse_args()

    cfg = SPECIES[args.species]
    repo = Path(args.surfaceome_repo).expanduser()
    cache = Path(args.cache_dir).expanduser()
    seqdir = repo / "data/external/deeptmhmm_surfaceome_inputs"

    # ---- the human reference sequence, pooled the same way the universe was built ------
    human_up = load_uniprot_proteoforms(cache / "uniprot_human_reviewed_isoforms.tsv.gz")
    human_cohort = parse_uniprot_fasta(seqdir / "human_canonical_non_hla.fasta")
    human_pep = parse_ensembl_pep(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    by_len_gene: dict[tuple[str, int], str] = {}
    for gene, forms in human_up.items():
        for p in forms:
            by_len_gene.setdefault((gene, len(p["seq"])), p["seq"])
    for r in human_cohort.values():
        if r["symbol"]:
            by_len_gene.setdefault((r["symbol"], len(r["seq"])), r["seq"])

    universe = [r for r in csv.DictReader(open(args.universe), delimiter="\t")]
    ens_by_gene_h = translations_by_gene(cache / "Homo_sapiens.GRCh38.pep.all.fa.gz")
    human_seq: dict[str, str] = {}
    for r in universe:
        if r["is_reference"] != "True":
            continue
        gene, want = r["gene"], int(r["length"])
        seq = by_len_gene.get((gene, want))
        if seq is None:
            # Ensembl-only reference: recover it from the gene's translations by length.
            for ensg in (r.get("ensembl_gene_id") or "").split(";"):
                for s in sorted(ens_by_gene_h.get(ensg.strip(), ())):
                    if len(s) == want:
                        seq = s
                        break
                if seq:
                    break
        if seq is None:
            for enst in (r.get("ensembl_transcripts") or "").split(";"):
                s = human_pep.get(enst.strip())
                if s and len(s) == want:
                    seq = s
                    break
        if seq:
            human_seq[gene] = seq
    print(f"human references with a sequence: {len(human_seq)}")

    # The ortholog metadata records the human gene it was found from as an ENSG. Joining on
    # that rather than on query_input_gene_symbols recovers the orthologs whose stored
    # symbol differs from the library's -- 565 of 1,971 mouse rows on the symbol join.
    ensg_for_acc: dict[str, set[str]] = defaultdict(set)
    with (repo / "data/external/uniprot_ensembl_xrefs/ensg_to_uniprot.tsv").open() as fh:
        next(fh)
        for line in fh:
            g, a = line.split("\t")[:2]
            ensg_for_acc[a].add(g)
    seq_by_human_ensg: dict[str, str] = {}
    symbol_by_human_ensg: dict[str, str] = {}
    for r in universe:
        if r["is_reference"] != "True" or r["gene"] not in human_seq:
            continue
        base = (r.get("uniprot") or "").split("-")[0]
        for ensg in ensg_for_acc.get(base, ()):
            seq_by_human_ensg.setdefault(ensg, human_seq[r["gene"]])
            symbol_by_human_ensg.setdefault(ensg, r["gene"])
    print(f"  reachable by human ENSG: {len(seq_by_human_ensg)}")

    # ---- ortholog candidates ------------------------------------------------------
    cohort = parse_uniprot_fasta(seqdir / f"{cfg['cohort']}.fasta")
    meta = [
        m for m in csv.DictReader((seqdir / f"{cfg['cohort']}_metadata.csv").open())
        if m.get("status") == "ok"
    ]
    ens_by_gene = translations_by_gene(cache / cfg["pep"])
    # Real protein ids for the Ensembl candidates: a synthetic key cannot be looked up
    # later, and the D1 table for Ensembl-only forms is keyed on the ENSP.
    ensp_for_seq: dict[str, str] = {}
    pep_all = parse_ensembl_pep(cache / cfg["pep"])
    import gzip as _gzip, re as _re

    with _gzip.open(cache / cfg["pep"], "rt") as fh:
        cur_id, buf = None, []
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    ensp_for_seq.setdefault("".join(buf).rstrip("*"), cur_id)
                cur_id = line[1:].split()[0].split(".")[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id:
            ensp_for_seq.setdefault("".join(buf).rstrip("*"), cur_id)
    if args.species == "cyno":
        up_by_symbol, reviewed_accs = load_cyno_uniprot(cache / cfg["uniprot"])
    else:
        up = load_uniprot_proteoforms(cache / cfg["uniprot"])
        up_by_symbol = up
        reviewed_accs = {p["accession"] for v in up.values() for p in v}

    appris_seq_for: dict[str, str] = {}
    if cfg["appris"]:
        appris = parse_appris(cache / cfg["appris"])
        pep = parse_ensembl_pep(cache / cfg["pep"])
        by_sym: dict[str, list[dict]] = defaultdict(list)
        for rows in appris.values():
            for r in rows:
                by_sym[r["symbol"]].append(r)
        for sym, rows in by_sym.items():
            sel = select_principal(rows, pep)
            if sel and sel["seq"]:
                appris_seq_for[sym] = sel["seq"]

    # Topology for whatever ortholog forms have a call: the cohort tables, plus any run of
    # re-picked forms. A candidate with no call is never rejected -- absence of evidence is
    # not evidence, the same rule the human guard uses.
    label_for_seq: dict[str, str] = {}
    tsv = repo / "data/processed/deeptmhmm" / f"deeptmhmm_{args.species}_ortholog.tsv"
    if tsv.exists():
        seq_for_acc = {a: r["seq"] for a, r in cohort.items()}
        with tsv.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = row.get("uniprot_accession_full") or row.get("uniprot_accession")
                seq = seq_for_acc.get(acc) or seq_for_acc.get((acc or "").split("-")[0])
                if seq:
                    label_for_seq.setdefault(seq, row["deeptmhmm_label"])
    for name in ("ortholog_repick_changed.3line", "ortholog_missing_topology.3line"):
        prior = Path(name)
        if not prior.exists():
            continue
        lines = prior.read_text().splitlines()
        for i in range(len(lines) - 2):
            if lines[i].startswith(">"):
                _n, lab = lines[i][1:].rsplit("|", 1)
                label_for_seq.setdefault(lines[i + 1].strip(), lab.strip())
    gpi_human = gpi_gene_set()
    # NetGPI per candidate. The TM guard refuses to move a reference onto a form that lost
    # every helix; for a GPI gene the helix count is zero either way, so the anchor is the
    # equivalent test and it needs a call on each candidate's own sequence.
    anchor_call: dict[str, bool] = {}
    netgpi_path = Path(args.netgpi)
    if netgpi_path.exists():
        for line in netgpi_path.read_text().splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 6:
                continue
            anchored = f[2] == "GPI-Anchored"
            # A low-confidence negative is no opinion, not a rejection: every genuine false
            # negative on the control set sat below this, and NetGPI misses whole families
            # (every CRIPTO form at 0.745, PMID 21055389).
            if not anchored and float(f[4]) < 0.55:
                continue
            anchor_call[f[0]] = anchored
        print(f"NetGPI candidate calls loaded: {len(anchor_call)}")
    print(f"ortholog topology labels available: {len(label_for_seq)}")

    def surface_ok(seq: str, human_symbol: str) -> bool | None:
        lab = label_for_seq.get(seq)
        if lab is None:
            return None
        # Keyed on the human gene: the architecture attribute is annotated per human gene,
        # and an ortholog of an anchored protein is normally anchored too. Normally, not
        # always — human nyctalopin is GPI-anchored where the mouse ortholog is not
        # (PMID 15905181) — so this is a guard against moving onto a soluble form, not a
        # claim about the ortholog's own anchor.
        return lab in SURFACE_LABELS or human_symbol in gpi_human

    def anchor_ok(acc: str, human_symbol: str) -> bool | None:
        """For a GPI gene, does this candidate still carry the anchor? None if unknown.

        Without this the guard is a no-op on exactly the genes it matters for: a GPI gene
        satisfies ``human_symbol in gpi_human``, so every candidate passed whatever its own
        sequence said.
        """
        if human_symbol not in gpi_human:
            return None
        return anchor_call.get(f"{args.species}_{human_symbol}_{acc}")

    def candidate_ok(acc: str, seq: str, human_symbol: str) -> bool | None:
        anchored = anchor_ok(acc, human_symbol)
        if anchored is not None:
            return anchored
        return surface_ok(seq, human_symbol)

    out, skipped, guarded = [], Counter(), []
    gpi_candidates: list[tuple[str, str]] = []
    for m in meta:
        q_ensg = (m.get("query_ensembl_gene_id") or "").strip()
        ref = seq_by_human_ensg.get(q_ensg)
        human_symbol = symbol_by_human_ensg.get(q_ensg) or (
            m["query_input_gene_symbols"] or ""
        ).split(";")[0].strip()
        if not ref:
            ref = human_seq.get(human_symbol)
        if not ref:
            skipped["no human reference sequence"] += 1
            continue
        cur_acc = m["selected_uniprot_accession"]
        cur = cohort.get(cur_acc)
        if not cur:
            skipped["current pick missing from the cohort FASTA"] += 1
            continue
        sym = m["target_gene_symbol"]

        cands: dict[str, dict] = {
            cur_acc: {"seq": cur["seq"], "source": "current", "reviewed": cur_acc in reviewed_accs}
        }
        for p in up_by_symbol.get(sym, []):
            cands.setdefault(
                p["accession"],
                {"seq": p["seq"], "source": "uniprot", "reviewed": p["accession"] in reviewed_accs},
            )
        for ensg in {m["target_ensembl_gene_id"]} - {""}:
            # Sorted: translations_by_gene hands back a set, and iterating a set of
            # strings follows this process's hash order, which made the tie-break between
            # equally-scoring candidates differ from run to run (TMIGD2 flipped between two
            # candidates that both score 0.8489).
            for s in sorted(ens_by_gene.get(ensg, ())):
                cands.setdefault(
                    ensp_for_seq.get(s, f"ENSP_UNKNOWN:{len(s)}"),
                    {"seq": s, "source": "ensembl", "reviewed": False},
                )
        appris_s = appris_seq_for.get(sym)

        if args.dump_gpi_candidates and human_symbol in gpi_human:
            for a, c in cands.items():
                gpi_candidates.append((f"{args.species}|{human_symbol}|{a}", c["seq"]))
            continue
        scores = {a: score(ref, c["seq"]) for a, c in cands.items()}
        # Drop candidates that are explicitly non-surface, unless nothing else survives.
        allowed = {
            a: v for a, v in scores.items()
            if candidate_ok(a, cands[a]["seq"], human_symbol) is not False
        }
        rejected = sorted(set(scores) - set(allowed))
        if allowed:
            scores_eff = allowed
        else:
            scores_eff = scores
            rejected = []
        best = max(scores_eff.values())
        tied = [a for a, v in scores_eff.items() if best - v <= TIE]
        if len(tied) > 1:
            tied.sort(
                key=lambda a: (
                    not cands[a]["reviewed"],
                    appris_s is None or cands[a]["seq"] != appris_s,
                    a != cur_acc,
                )
            )
        pick = tied[0]
        if rejected and cur_acc in rejected:
            guarded.append((human_symbol, sym, cur_acc, pick))
        out.append(
            {
                "human_symbol": human_symbol,
                "ortholog_symbol": sym,
                "human_ref_length": len(ref),
                "current_acc": cur_acc,
                "current_length": len(cur["seq"]),
                "current_score": round(scores[cur_acc], 4),
                "current_reviewed": m["selected_reviewed"],
                "picked_acc": pick,
                "picked_length": len(cands[pick]["seq"]),
                "picked_score": round(scores[pick], 4),
                "picked_source": cands[pick]["source"],
                "picked_reviewed": cands[pick]["reviewed"],
                "n_candidates": len(cands),
                "n_tied": len(tied),
                "changed": pick != cur_acc,
                "score_gain": round(scores[pick] - scores[cur_acc], 4),
                "soluble_guard_rejected": ";".join(rejected),
                "picked_anchor_call": {True: "anchored", False: "not_anchored"}.get(
                    anchor_ok(pick, human_symbol), ""
                ),
                "picked_topology_label": label_for_seq.get(cands[pick]["seq"], ""),
                # The guard only declines to MOVE onto a soluble form. This says whether the
                # form we end up with is itself soluble and un-anchored, which is the thing
                # the library would have to drop.
                "restrict_soluble_not_gpi": (
                    label_for_seq.get(cands[pick]["seq"]) in NON_SURFACE_LABELS
                    and human_symbol not in gpi_human
                ),
                "picked_sequence": cands[pick]["seq"],
            }
        )

    if args.dump_gpi_candidates:
        seen: set[str] = set()
        with open(args.dump_gpi_candidates, "w") as fh:
            for name_, seq_ in gpi_candidates:
                if name_ in seen:
                    continue
                seen.add(name_)
                fh.write(f">{name_}\n{seq_.rstrip('*')}\n")
        print(f"wrote {args.dump_gpi_candidates} ({len(seen)} candidates)")
        return

    if args.dump_picked:
        with open(args.dump_picked, "w") as fh:
            for r in out:
                fh.write(f">{args.species}|{r['human_symbol']}|{r['picked_acc']}\n"
                         f"{r['picked_sequence'].rstrip('*')}\n")
        print(f"wrote {args.dump_picked} ({len(out)} picks)")

    name = f"ortholog_repick_{args.species}"
    cols = [c for c in out[0] if c != "picked_sequence"]
    with open(f"{name}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    # Only the changed picks need a topology call; an unchanged pick keeps the one it has.
    changed_rows = [r for r in out if r["changed"]]
    Path(f"{name}_changed.fasta").write_text(
        "".join(
            f">{args.species}|{r['ortholog_symbol']}|{r['picked_acc']}\n{r['picked_sequence']}\n"
            for r in changed_rows
        )
    )
    print(f"wrote {name}_changed.fasta ({len(changed_rows)} sequences needing topology)")
    print(f"wrote {name}.tsv ({len(out)} orthologs)")
    for k, n in skipped.most_common():
        print(f"  skipped — {k}: {n}")

    print(f"  soluble guard rejected a candidate on {sum(1 for r in out if r['soluble_guard_rejected'])} orthologs")
    ch = [r for r in out if r["changed"]]
    lr = lambda r, k: min(r[k], r["human_ref_length"]) / max(r[k], r["human_ref_length"])
    summary = {
        "species": args.species,
        "orthologs": len(out),
        "changed": len(ch),
        "median_score_current": round(statistics.median(r["current_score"] for r in out), 4),
        "median_score_picked": round(statistics.median(r["picked_score"] for r in out), 4),
        "picked_source": dict(Counter(r["picked_source"] for r in out).most_common()),
        "picked_reviewed": sum(1 for r in out if r["picked_reviewed"]),
        "current_reviewed": sum(1 for r in out if r["current_reviewed"] == "reviewed"),
        "median_length_ratio_current": round(
            statistics.median(lr(r, "current_length") for r in out), 4
        ),
        "median_length_ratio_picked": round(
            statistics.median(lr(r, "picked_length") for r in out), 4
        ),
        "gain_over_5pct": sum(1 for r in out if r["score_gain"] > 0.05),
        "soluble_guard_rejected_a_candidate": sum(
            1 for r in out if r["soluble_guard_rejected"]
        ),
        "never_worse": all(r["score_gain"] >= -TIE for r in out),
    }
    Path(f"{name}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
