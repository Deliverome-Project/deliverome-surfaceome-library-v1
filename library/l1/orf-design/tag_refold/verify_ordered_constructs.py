"""Prove every refold is the protein that will be ordered.

Becca, 2026-11-12: "alphafold repredict and save the structure with the diversified gggs but
constant alfa and spy". Each L1 member's cassette DNA is diversified (member-specific
linker codons, ``l1_cassette_encodings.tsv``), so the amino acids must still be the p003
51-aa tandem. This translates every built insert on the draft-L1 integration branch, checks

1. the translated protein's SHA-256 equals the SHA-256 its rung table states;
2. the cassette amino acids are exactly the p003 tandem, once (asserted, not assumed);
3. the insert carries the member's own cassette encoding where one is listed;
4. the protein (or, for Igk designs folded without the cleaved leader, its mature chain) is
   a sequence this folder folded -- by SHA-256, so any codon change that altered a
   residue would show as unfolded,

and writes one row per built construct to ``tag_refold_construct_check.tsv``. Constructs
whose protein was never folded are written to ``tag_refold_ordered_unfolded.tsv``, which
``build_refold_inputs.py`` queues at tier 1.

Implemented by Claude Code session l1-tag-refold (coordinator request, 2026-11-12).

    uv run --locked python library/l1/orf-design/tag_refold/verify_ordered_constructs.py \\
        --folds ~/.cache/deliverome-analysis/tagged-folds/tagged-2026-11-12
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
REPO = EXP.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EXP))

import build_gpcr_rescue_arms as gpcr_arms  # noqa: E402
from build_refold_inputs import CASSETTE, CODON, translate_insert  # noqa: E402


def _sibling(name: str):
    """Load this folder's own ``name``.py: S5E13 has a module of the same name, and whichever is
    imported first otherwise wins ``sys.modules`` in a full test run."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(f"tag_refold_{name}", HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


read_queue = _sibling("refold_queue").read_queue

from deliverome_analysis.cloning import seqtools  # noqa: E402

#: Draft-L1 integration branch (PR #299) at the commit verified, and l1-reviewer-package's
#: encodings table (PR #283) with the member-specific diversified cassettes.
INTEGRATION_COMMIT = "5cfd91c02d2d7e889097e55b48d874277ff1306e"
ENCODINGS_COMMIT = "59050daa"
#: Overnight review (PR #303): Twist synthesis-compliance recodes, codon by codon.
RECODES_COMMIT = "8d4b935e22b10357b6cb45cacd615fc4c8a3c71e"
RECODES_OUT = EXP / "tag_refold_recode_check.tsv"
#: Historical 2026-11-13 order verified by this refold check (superseded for current orders by
#: ``synthesis_compliance.FASTA_KEY``).
ORDER_KEY = "sprint4/s4e9-orf-design-fg-panel/synthesis/2026-11-13-one-arm-renamed/l1_order_inserts_recoded.fasta"
ORDER_OUT = EXP / "tag_refold_order_check.tsv"
IGK_LEADER = gpcr_arms.LEADER_AA
EXP_REL = "library/l1/orf-design"
RUNGS = ("cterm", "nterm", "internal", "multipass", "type3")
CHECK_OUT = EXP / "tag_refold_construct_check.tsv"
UNFOLDED_OUT = EXP / "tag_refold_ordered_unfolded.tsv"


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def show(commit: str, rel: str) -> str:
    from deliverome_analysis.core.provenance import read_pinned_text

    return read_pinned_text(commit, f"{EXP_REL}/{rel}", repo=REPO)


def rows(text: str) -> list[dict[str, str]]:
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def fasta(text: str) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name: str | None = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = seqtools.fasta_id(line[1:])
            out[name] = []
        elif name is not None:
            out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", required=True)
    args = ap.parse_args()
    folds = Path(args.folds).expanduser()

    queue = read_queue(EXP / "tag_refold_queue.tsv")
    folded = {
        q["tagged_sha256"]: q["fold_id"] for q in queue if (folds / f"{q['fold_id']}.pdb").exists()
    }
    enc = {
        r["candidate_id"]: r["cassette_nt"]
        for r in rows(show(ENCODINGS_COMMIT, "l1_cassette_encodings.tsv"))
    }

    out, unfolded = [], []
    for rung in RUNGS:
        table = rows(show(INTEGRATION_COMMIT, f"constructs_l1_{rung}.tsv"))
        inserts = fasta(show(INTEGRATION_COMMIT, f"constructs_l1_{rung}.fasta"))
        for r in table:
            key = r.get("construct_id") or r["candidate_id"]
            nt = inserts.get(key) or inserts.get(r["candidate_id"]) or inserts.get(r["design_id"])
            if nt is None:  # held / not built in this rung table
                continue
            prot = translate_insert(nt)
            stated = r.get("protein_sha256", "")
            n_cas = prot.count(CASSETTE)
            # residues the rung folds without: the cleaved Igk leader (multipass, type III Igk)
            mature = ""
            if r.get("mature_aa") and r.get("mature_aa") != r.get("protein_aa"):
                mature = prot[len(prot) - int(r["mature_aa"]) :]
            fold_id = folded.get(sha(prot)) or (folded.get(sha(mature)) if mature else None)
            e = enc.get(r["candidate_id"])
            row = {
                "rung": rung, "construct": key, "design_id": r["design_id"], "candidate_id": r["candidate_id"],
                "gene_symbol": r["gene_symbol"], "protein_sha256": sha(prot),
                "matches_stated_sha256": "yes" if stated and sha(prot) == stated else ("no_stated" if not stated else "NO"),
                "cassette_aa_constant": "yes" if n_cas == 1 else f"NO ({n_cas} exact copies)",
                "member_encoding": ("yes" if e in nt else "NO") if e else "none listed (shared p003 encoding)",
                "folded_form": ("precursor" if folded.get(sha(prot)) else "mature chain") if fold_id else "",
                "fold_id": fold_id or "",
                "status": "ok" if fold_id and n_cas == 1 and (not stated or sha(prot) == stated) and (not e or e in nt) else "CHECK",
            }  # fmt: skip
            out.append(row)
            if not fold_id and n_cas == 1:
                # native where the construct shows it: a snorkel is native + module + cassette
                # at the C-terminus; any other single insertion is the protein minus the cassette
                native = ""
                if "snorkel" in r.get("placement", ""):
                    native = prot[: int(r["native_aa"])]
                elif not mature and r.get("arm", "native_site") == "native_site":
                    i = prot.index(CASSETTE)
                    native = prot[:i] + prot[i + len(CASSETTE) :]
                unfolded.append({"rung": rung, "construct": key, "design_id": r["design_id"],
                                 "gene_symbol": r["gene_symbol"],
                                 "native_sha256": r.get("native_sha256") or (sha(native) if native else ""),
                                 "fold_form": "mature" if mature else "precursor",
                                 "native_seq": native, "tagged_seq": mature or prot})  # fmt: skip
    cols = list(out[0])
    with CHECK_OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    with UNFOLDED_OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rung", "construct", "design_id", "gene_symbol", "native_sha256", "fold_form", "native_seq", "tagged_seq"],
                           delimiter="\t", lineterminator="\n")  # fmt: skip
        w.writeheader()
        w.writerows(unfolded)
    print(f"{len(out)} built constructs checked")
    for k, n in sorted(Counter((o["rung"], o["status"]) for o in out).items()):
        print(f"  {k[0]:<10} {k[1]:<6} {n}")
    for col in ("matches_stated_sha256", "cassette_aa_constant", "member_encoding"):
        print(
            f"  {col}: {dict(Counter(o[col] if not o[col].startswith('none') else 'none listed' for o in out))}"
        )
    print(f"  not folded: {len(unfolded)} -> {UNFOLDED_OUT.name}")
    check_recodes()
    check_order_of_record(folds)


def check_order_of_record(folds: Path) -> None:
    """Every insert in the order of record: protein translated, cassette amino acids exactly
    once, and a LIVE (non-superseded) fold of that protein -- the precursor, or the mature
    chain of an Igk design -- with a verdict."""
    from deliverome_analysis.core.s3 import processed

    path = Path(processed.download(ORDER_KEY))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    queue = read_queue(EXP / "tag_refold_queue.tsv")
    live = {q["tagged_sha256"]: q["fold_id"] for q in queue
            if q["kind"] != "superseded" and (folds / f"{q['fold_id']}.pdb").exists()}  # fmt: skip
    verdict = {
        r["tagged_sha256"]: r["verdict"]
        for r in csv.DictReader((EXP / "tag_refold_metrics.tsv").open(), delimiter="\t")
        if r.get("tagged_sha256")
    }
    out = []
    text = path.read_text()
    heads = [line[1:] for line in text.splitlines() if line.startswith(">")]
    seqs = fasta(text)
    for head in heads:
        name = head.split("|")[0].split()[0]
        parts = head.split()
        prot = translate_insert(seqs[name])
        forms = [("precursor", prot)] + (
            [("mature chain", prot[len(IGK_LEADER) :])] if prot.startswith(IGK_LEADER) else []
        )
        form, fid = next(((f, live[sha(sq)]) for f, sq in forms if sha(sq) in live), ("", ""))
        hit = next((sq for f, sq in forms if f == form), prot)
        out.append({"insert": name, "design_id": parts[-2], "rung": parts[-1], "protein_sha256": sha(prot),
                    "cassette_aa_constant": "yes" if prot.count(CASSETTE) == 1 else "NO",
                    "folded_form": form, "fold_id": fid, "verdict": verdict.get(sha(hit), "") if fid else "",
                    "status": "ok" if fid and prot.count(CASSETTE) == 1 else "CHECK"})  # fmt: skip
    with ORDER_OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"order of record {ORDER_KEY.split('/')[-2]} (sha256 {digest[:12]}): {len(out)} inserts")
    for k, n in sorted(Counter(o["status"] for o in out).items()):
        print(f"  {k}: {n}")
    print(f"  folded as: {dict(Counter(o['folded_form'] or 'NOT FOLDED' for o in out))}")
    print(f"  verdicts: {dict(Counter(o['verdict'] or '-' for o in out))}")


def check_recodes() -> None:
    """Apply #303's Twist recodes to the #299 inserts and prove no residue changes.

    Each recode names a construct, a residue (1-based in the protein), the amino acid and
    the from/to codons. Checks: both codons encode that amino acid; the insert really has
    ``from_codon`` at that residue; and the re-translated protein's SHA-256 is unchanged.
    """
    recodes = rows(show(RECODES_COMMIT, "synthesis_recodes.tsv"))
    inserts: dict[tuple[str, str], str] = {}
    for rung in RUNGS:
        for key, nt in fasta(show(INTEGRATION_COMMIT, f"constructs_l1_{rung}.fasta")).items():
            inserts[(rung, key)] = nt
    by_construct: dict[tuple[str, str], list[dict[str, str]]] = {}
    for r in recodes:
        by_construct.setdefault((r["rung"], r["construct"].split("|")[0].split()[0]), []).append(r)
    out = []
    for (rung, key), rs in sorted(by_construct.items()):
        nt = inserts.get((rung, key))
        row = {"rung": rung, "construct": key, "design_id": rs[0]["design_id"], "gene_symbol": rs[0]["gene_symbol"],
               "n_recodes": len(rs)}  # fmt: skip
        if nt is None:
            out.append(
                {
                    **row,
                    "status": "construct not on #299",
                    "protein_unchanged": "",
                    "synonymous": "",
                    "from_codon_found": "",
                }
            )
            continue
        start = nt.find("GCCGCCACC") + 9
        before = translate_insert(nt)
        seq = list(nt)
        syn = found = True
        for r in rs:
            i = start + 3 * (int(r["residue"]) - 1)
            aa = r["amino_acid"]
            syn &= CODON.get(r["from_codon"]) == aa == CODON.get(r["to_codon"])
            found &= nt[i : i + 3] == r["from_codon"]
            seq[i : i + 3] = list(r["to_codon"])
        after = translate_insert("".join(seq))
        same = sha(after) == sha(before)
        out.append({**row, "synonymous": "yes" if syn else "NO", "from_codon_found": "yes" if found else "NO",
                    "protein_unchanged": "yes" if same else "NO",
                    "status": "ok" if syn and found and same else "CHECK"})  # fmt: skip
    cols = [
        "rung",
        "construct",
        "design_id",
        "gene_symbol",
        "n_recodes",
        "synonymous",
        "from_codon_found",
        "protein_unchanged",
        "status",
    ]
    with RECODES_OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(out)
    print(f"{len(recodes)} Twist recodes over {len(out)} constructs (#303 {RECODES_COMMIT[:8]}):")
    for k, n in sorted(Counter(o["status"] for o in out).items()):
        print(f"  {k}: {n}")


if __name__ == "__main__":
    main()
