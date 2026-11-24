"""Audit: L1 library references whose designed protein starts downstream of the gene's signal peptide.

Report only; nothing in the library changes. Requested by Becca on 2026-11-11, after the type III
second review (shared notes ``type3-risk-review.md``; PR 282) found three references (CD200R1L,
LRRC37B, CACHD1) whose ordered form starts after its gene's UniProt signal peptide while the
signal-peptide form sits in L2. Implemented by the Claude Code session l1-reviewer-package, PR 283.

For every active L1 ``library_reference`` member:

* the gene's UniProt canonical is ``current_accession`` in the pinned ``reference_selection.tsv``;
  its signal peptide (``ft_signal``) and sequence come from the UniProt REST API, queried once
  and saved beside the output as ``reference_signal_peptides_uniprot.tsv`` with the release
  UniProt reports, so the audit can be rerun against that saved copy (``--offline``);
* the member's protein is its designed CDS translated (``library_design.designed_proteins``);
* it is flagged when the canonical has a signal peptide ending at residue *e* and the member's
  protein does not contain the canonical's first *e* residues, i.e. it lacks the signal peptide.
  How it lacks it is reported: ``starts_downstream`` when the member's N-terminal 15 residues
  after Met1 occur in the canonical at or after residue *e* - 1 (an N-terminal truncation), and
  ``divergent_n_terminus`` otherwise (a different first exon, or a sequence that does not align).
  The member's DeepTMHMM v1 signal-peptide call (``l1_tag_paths.tsv``) is reported beside it.

    uv run --locked python audit_reference_signal_peptides.py
    uv run --locked python audit_reference_signal_peptides.py --offline   # from the saved UniProt copy
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import urllib.request
from pathlib import Path

from deliverome_analysis import library_design as ld
from deliverome_analysis.s3 import processed

HERE = Path(__file__).parent
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
TAG_PATHS = HERE / "l1_tag_paths.tsv"
REFSEL_KEY = "sprint4/s4e9-orf-design-fg-panel/library-design/reference_selection.tsv"
UNIPROT_COPY = HERE / "reference_signal_peptides_uniprot.tsv"
OUT = HERE / "l1_reference_signal_peptide_audit.tsv"
URL = (
    "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,ft_signal,sequence"
    "&query=(reviewed:true)%20AND%20(organism_id:9606)%20AND%20(ft_signal:*)"
)
PROBE = 15


def fetch_uniprot() -> list[dict]:
    with urllib.request.urlopen(URL, timeout=300) as resp:
        release = resp.headers.get("X-UniProt-Release", "")
        text = resp.read().decode()
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    out = []
    for r in rows:
        m = re.search(r"SIGNAL (\d+)\.\.(\d+)", r["Signal peptide"])
        if m:
            out.append(
                {
                    "accession": r["Entry"],
                    "sp_end": m.group(2),
                    "sequence": r["Sequence"],
                    "uniprot_release": release,
                }
            )
    with UNIPROT_COPY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--refresh",
        action="store_true",
        help="re-query the live UniProt REST API and overwrite the saved copy",
    )
    ap.add_argument(
        "--offline",
        action="store_true",
        help="use the saved UniProt copy (the default; kept for compatibility)",
    )
    args = ap.parse_args(argv)

    uni = (
        fetch_uniprot()
        if args.refresh and not args.offline
        else list(csv.DictReader(UNIPROT_COPY.open(), delimiter="\t"))
    )
    sp = {r["accession"]: r for r in uni}
    refsel = {
        r["symbol"]: r for r in csv.DictReader(open(processed.download(REFSEL_KEY)), delimiter="\t")
    }
    proteins = ld.designed_proteins(MANIFEST, library=ld.LIBRARY_L1)
    design = ld.design_ids(MANIFEST)
    v1 = {r["candidate_id"]: r for r in csv.DictReader(TAG_PATHS.open(), delimiter="\t")}
    rows = []
    for cid, p in sorted(proteins.items()):
        if p["candidate_set"] != "library_reference":
            continue
        acc = refsel.get(p["gene_symbol"], {}).get("current_accession", "")
        u = sp.get(acc)
        if not u:
            continue
        canon, end, seq = u["sequence"], int(u["sp_end"]), p["sequence"].rstrip("*")
        if canon[:end] in seq:
            continue
        probe = seq[1 : 1 + PROBE]
        pos = canon.find(probe)
        kind = "starts_downstream" if pos >= end - 1 else "divergent_n_terminus"
        rows.append(
            {
                "candidate_id": cid,
                "design_id": design.get(cid, cid),
                "gene_symbol": p["gene_symbol"],
                "uniprot_canonical": acc,
                "uniprot_sp_end": end,
                "member_length_aa": len(seq),
                "canonical_length_aa": len(canon),
                # canonical residue aligned to the member's Met1 (1-based), blank when it does not align
                "member_start_in_canonical": pos if pos >= 0 else "",
                "kind": kind,
                "deeptmhmm_v1_signal_peptide_aa": v1.get(cid, {}).get("signal_peptide_aa", ""),
                "deeptmhmm_v1_placement": v1.get(cid, {}).get("placement", ""),
                "uniprot_release": u["uniprot_release"],
                "decided_by": "no decision taken (report only)",
                "implemented_by": "Claude Code session l1-reviewer-package, PR #283",
                "reason": "Becca 2026-11-11: audit every L1 reference that starts downstream of its gene's UniProt signal peptide",
            }
        )
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            list(rows[0]) if rows else ["candidate_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)
    by = {}
    for r in rows:
        by[r["kind"]] = by.get(r["kind"], 0) + 1
    print(f"{len(rows)} L1 references lack their gene's UniProt signal peptide: {by}")
    for r in rows:
        print(
            f"  {r['candidate_id']} {r['gene_symbol']:<10} {r['kind']:<22} SP 1-{r['uniprot_sp_end']}, "
            f"member {r['member_length_aa']} aa vs {r['canonical_length_aa']}; "
            f"v1 SP {r['deeptmhmm_v1_signal_peptide_aa']}, {r['deeptmhmm_v1_placement']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
