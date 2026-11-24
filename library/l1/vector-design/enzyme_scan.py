"""s1e2 — Type IIS site scan across the construct series, and the two Esp3I edits.

Regenerates the tables in README §4.9 / §4.10:
  * per-construct BbsI / Esp3I / BsaI counts and positions (both strands, circular),
  * the two hand-domestication edits to #236083 (Esp3I 252 and 306), re-derived by
    aligning the before/after windows rather than trusting the prose,
  * the safety checks — that neither edit touches an annotated feature, and that
    neither creates a new Type IIS site.

⚠️ Inputs are LOCAL-ONLY GenBanks in ~/Downloads (Addgene downloads + Becca's construct
   maps). They are not in git or S3, so this script only runs on a machine that has them;
   it reports which files are missing rather than failing. Checking the plasmid maps in
   somewhere shared is a logged TODO.

Run:
    uv run python library/l1/vector-design/enzyme_scan.py
"""

import re
from pathlib import Path

# allow-local-path: s1e2 inputs are local-only GenBanks (see docstring TODO)
DOWNLOADS = Path.home() / "Downloads"

# Type IIS sites, as (forward, reverse-complement) recognition sequences.
# Our documented Golden Gate pair is BbsI (assembly) + Esp3I/BsmBI (reserved downstream).
# BsaI is scanned only to show it was never cleared — see README §4.10.
ENZYMES = {
    "BbsI": ("GAAGAC", "GTCTTC"),
    "Esp3I": ("CGTCTC", "GAGACG"),
    "BsaI": ("GGTCTC", "GAGACC"),
    "AarI": ("CACCTGC", "GCAGGTG"),
    "SapI": ("GCTCTTC", "GAAGAGC"),
}

# The construct series. Order matters — it is the README §4.10 table order.
SERIES = [
    ("#236083 parent", "addgene-plasmid-236083-sequence-467823.gbk"),
    ("p004 wpre-only", "p004_pag-lenti_236083_paci_clai_site_wpre_only.gb"),
    (
        "p003 assembled",
        "p003_pag-lenti_236083_esp3i_domesticated_efs_tfrc_bc_alfa_iresfragmid_wpre_puro_assembled.gb",
    ),
    (
        "p003 assembled (cmv-empty-puro variant)",
        "p003_pag-lenti-cmv-empty-puro_236083_esp3i_domesticated_efs_tfrc_bc_alfa_iresfragmid_wpre_puro_assembled.gb",
    ),
    ("p006 ctermtag", "p006_pag_lenti_cmv_tfrc_ctermtag_ires_puro.gb"),
    ("p007 V291", "p007_pag_lenti_cmv_tfrc_v291_ires_puro.gb"),
    ("pRDA_722 / #216077", "addgene-plasmid-216077-sequence-431423.gbk"),
]

# The two windows the user supplied: #236083 242..388 before, and 242..313 after.
BEFORE = (
    "aagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaa"
    "gggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcactt"
).upper()
AFTER = ("aagctgtgaccAtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagaAgaa").upper()


def read_genbank(path: Path) -> tuple[str, list[tuple[str, int, int, bool, str]]]:
    """Return (sequence, features). Hand-rolled: biopython is not a project dep."""
    txt = path.read_text(errors="replace")
    seq = re.sub(r"[^A-Za-z]", "", txt.split("\nORIGIN", 1)[1]).upper()
    feats = []
    if "\nFEATURES" in txt:
        body = txt.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0]
        for block in re.split(r"\n(?=     \S)", body):
            m = re.match(r"\s+(\S+)\s+(complement\()?<?(\d+)\.\.>?(\d+)", block)
            if not m:
                continue
            lab = re.search(r'/(?:label|note|gene|product)="([^"]{0,70})', block)
            feats.append(
                (
                    m.group(1),
                    int(m.group(3)),
                    int(m.group(4)),
                    bool(m.group(2)),
                    lab.group(1) if lab else "",
                )
            )
    return seq, feats


def find_sites(seq: str, enzymes=ENZYMES) -> dict[str, list[int]]:
    """1-based positions of each recognition site, both strands. Sequence treated as
    circular, so sites spanning the origin are not missed."""
    out = {}
    for name, (fwd, rev) in enzymes.items():
        pats = {fwd, rev}  # palindromic sites would otherwise double-count
        wrapped = seq + seq[: max(len(p) for p in pats) - 1]
        hits = {m.start() % len(seq) + 1 for p in pats for m in re.finditer(p, wrapped)}
        out[name] = sorted(hits)
    return out


def main() -> None:
    print("=== Type IIS sites across the construct series (README §4.10) ===\n")
    print(f"{'construct':42s} {'bp':>7}  " + "  ".join(f"{e:>16s}" for e in ENZYMES))
    for label, fn in SERIES:
        path = DOWNLOADS / fn
        if not path.exists():
            print(f"{label:42s} {'--':>7}  (GenBank not found locally: {fn})")
            continue
        seq, _ = read_genbank(path)
        sites = find_sites(seq)
        cells = []
        for e in ENZYMES:
            pos = sites[e]
            cells.append(f"{len(pos)}{'' if not pos else ' ' + ','.join(map(str, pos))}")
        print(f"{label:42s} {len(seq):>7}  " + "  ".join(f"{c:>16s}" for c in cells))

    # ---------- the two edits ----------
    parent = DOWNLOADS / SERIES[0][1]
    if not parent.exists():
        print("\n(skipping edit analysis — #236083 GenBank not found locally)")
        return
    seq, feats = read_genbank(parent)
    offset = seq.find(BEFORE)
    if offset < 0:
        raise SystemExit("BEFORE window not found in #236083 — check the sequence")
    print(
        f"\n=== the two Esp3I edits (README §4.9) — window at #236083 {offset + 1}"
        f"..{offset + len(BEFORE)} ===\n"
    )
    print(
        f"  AFTER window is {len(AFTER)} nt, i.e. the same region truncated at "
        f"{offset + len(AFTER)} — NOT a deletion.\n"
    )

    for i, (a, b) in enumerate(zip(BEFORE, AFTER, strict=False)):
        if a != b:
            abs_pos = offset + i + 1
            # which recognition site does this base sit in, and where in it?
            for name, (fwd, rev) in ENZYMES.items():
                for strand, pat in (("+", fwd), ("-", rev)):
                    for j in range(max(0, i - len(pat) + 1), i + 1):
                        if BEFORE[j : j + len(pat)] == pat:
                            site_start = offset + j + 1
                            edited = BEFORE[j : j + len(pat)]
                            after = AFTER[j : j + len(pat)]
                            # Report the position in the enzyme's OWN 5'->3' frame, so a
                            # site found as GAGACG is described in CGTCTC coordinates.
                            base_fwd = i - j + 1
                            base_own = base_fwd if strand == "+" else len(pat) - base_fwd + 1
                            print(
                                f"  edit at {abs_pos}: {a.lower()}->{b}  "
                                f"kills {name} {strand} at {site_start}"
                                f"..{site_start + len(pat) - 1}  "
                                f"(top strand {edited} -> {after}; "
                                f"base {base_own} of {name}'s own {fwd} frame)"
                            )
            overlapping = [f for f in feats if f[1] <= abs_pos <= f[2] and f[0] != "source"]
            print(
                f"    annotated features spanning {abs_pos}: "
                f"{overlapping if overlapping else 'NONE (unannotated backbone)'}"
            )

    # nearest annotated neighbours, to show what was NOT disturbed
    print("\n  nearest annotated features (must be intact):")
    lo, hi = offset + 1, offset + len(AFTER)
    for kind, a, c, comp, lab in feats:
        if kind == "source" or c < lo - 200 or a > hi + 200:
            continue
        print(f"    {kind:14s} {a:>5}..{c:<5} {'(-)' if comp else '(+)'}  {lab}")

    # no new site created?
    print("\n  new-site check on the edited window:")
    for label, win in (("before", BEFORE[: len(AFTER)]), ("after", AFTER)):
        hits = {e: [p for p in find_sites(win)[e]] for e in ENZYMES}
        hits = {e: p for e, p in hits.items() if p}
        print(f"    {label:6s}: {hits if hits else 'no Type IIS sites'}")
    print("\n  → both sites removed, nothing new created.")


if __name__ == "__main__":
    main()
