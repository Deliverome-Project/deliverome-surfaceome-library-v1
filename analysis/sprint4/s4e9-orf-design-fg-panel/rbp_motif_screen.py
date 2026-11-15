"""RBP-motif screen of the L1 library: shared elements once, ORF-intrinsic hits per member.

Every active L1 member is assembled in silico into the p0039 backbone with the order-time
insert grammar (Esp3I | Kozak | designed CDS + stop | pad | 5' handle | barcode | gap + RT |
10x CS1 | Esp3I), and its EFS transcript (EFS end to the 3' LTR R end) is scanned with
``deliverome_analysis.rbp_motifs`` -- the same code path the construct gate runs on a real
map. Hits are split three ways, because natural CDSs carry some motifs and the question for
the order is different for each:

- **shared**: outside the ORF, identical in every member (the 5' leader, the insert tail,
  the p0039 backbone transcript: IRES2, PuroR, WPRE, 3' LTR to R). Reported once.
- **junction**: overlapping the ORF boundary, so created by this ORF meeting the shared
  sequence. Designed sequence; reported per member.
- **ORF-intrinsic**: wholly inside the designed CDS. Reported per member as counts.

Barcodes are not yet assigned to L1 members, so the barcode slot is masked (``N``) in the
member scan and the s1e2 pool (``barcodes_20nt_30k.tsv``) is scanned on its own, in its
handle/RT context so motifs formed across the junctions are caught. Issued linker-only tag
cassettes (``l1_cassette_encodings.tsv``) are scanned on their own as designed sequence.

    uv run --locked python rbp_motif_screen.py

Inputs (S3): ``library_design.CONSTRUCTS_KEY`` (designed CDS),
``p0039_EFS_IRES_puro_backbone.gb``, ``sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv``.
Committed: ``candidate_manifest_tranches.tsv`` (active L1 = ``library == l1`` and not retired),
``l1_cassette_encodings.tsv``.

Outputs: ``l1_rbp_motif_shared.tsv``, ``l1_rbp_motif_members.tsv``,
``l1_rbp_motif_cassettes.tsv``, ``l1_rbp_motif_summary.json``.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from deliverome_analysis import library_design, plasmid_qc, rbp_motifs, s3

HERE = Path(__file__).parent
PREFIX = "sprint4/s4e9-orf-design-fg-panel/"
#: The current constructs record (append-only, versioned; the 2026-11-10 one adds the adopted
#: Compara 116 ortholog forms), so a member is never scanned from a superseded record.
CONSTRUCTS_KEY = library_design.CONSTRUCTS_KEY
BACKBONE_KEY = PREFIX + "p0039_EFS_IRES_puro_backbone.gb"
BARCODE_KEY = "sprint1/s1e2-fragmid-vector-design/barcodes_20nt_30k.tsv"
MANIFEST = HERE / "candidate_manifest_tranches.tsv"
CASSETTES = HERE / "l1_cassette_encodings.tsv"

OUT_SHARED = HERE / "l1_rbp_motif_shared.tsv"
OUT_MEMBERS = HERE / "l1_rbp_motif_members.tsv"
OUT_CASSETTES = HERE / "l1_rbp_motif_cassettes.tsv"
OUT_SUMMARY = HERE / "l1_rbp_motif_summary.json"

#: The p0039 slot (MCS + stuffer) the insert replaces; as in assemble_candidates.py.
SLOT_START, SLOT_LEN = 4808, 195
ESP3I_5, ESP3I_3 = "CGTCTCA", "AGAGAGACG"
MASK = "N" * plasmid_qc.BARCODE_LEN
#: Context either side of a barcode when scanning the pool: long enough to hold any motif
#: in the panel that crosses the junction, including the acceptor's 15-nt tract.
BARCODE_CONTEXT = 18

#: Choices this script makes on top of ``rbp_motifs.DECISIONS`` (convention: Becca,
#: 2026-11-10, note who decided, who implemented, and why).
DECISIONS = {
    "split_shared_junction_orf": (
        "Becca, 2026-11-08 (request relayed to the rbp-motif-screen session on 2026-11-10)",
        rbp_motifs.IMPLEMENTED_BY,
        "natural CDSs carry some motifs, so shared-element hits are reported apart from "
        "ORF-intrinsic ones; the junction class (hits spanning the ORF boundary) was added by "
        "the agent because those are designed sequence",
    ),
    "barcode_masked_pool_scanned_separately": (
        "Claude Code session rbp-motif-screen, 2026-11-10; no separate sign-off, shipped in "
        "PR #277 (merged by Becca 2026-11-10)",
        rbp_motifs.IMPLEMENTED_BY,
        "L1 members have no barcodes yet, so the slot is masked and the s1e2 pool is scanned "
        f"in a {BARCODE_CONTEXT}-nt handle/RT context on its own",
    ),
}

INSERT = (  # (feature type, label, sequence or None for the ORF / barcode)
    ("misc_feature", "Esp3I (shuttle out)", ESP3I_5),
    ("regulatory", "Kozak (AgeI-XbaI)", plasmid_qc.KOZAK),
    ("CDS", "ORF", None),
    ("misc_feature", "pad", plasmid_qc.SPACER),
    ("primer_bind", "5' handle", plasmid_qc.FWD_PRIMER),
    ("misc_feature", "barcode", MASK),
    ("primer_bind", "gap + RT site", plasmid_qc.REV_PRIMER),
    ("misc_feature", "10x Capture Sequence 1", plasmid_qc.CS1),
    ("misc_feature", "Esp3I (shuttle out)", ESP3I_3),
)


def active_l1() -> list[dict]:
    with MANIFEST.open() as fh:
        return [
            r
            for r in csv.DictReader(fh, delimiter="\t")
            if r["library"] == "l1" and not r["retired_reason"]
        ]


def member_record(backbone: SeqRecord, orf: str) -> SeqRecord:
    """p0039 with the slot replaced by this member's insert, features carried over."""
    bb = str(backbone.seq).upper()
    parts, feats, pos = [], [], SLOT_START
    for ftype, label, s in INSERT:
        s = orf if s is None else s
        feats.append(
            SeqFeature(
                FeatureLocation(pos, pos + len(s), 1), type=ftype, qualifiers={"label": [label]}
            )
        )
        parts.append(s)
        pos += len(s)
    insert = "".join(parts)
    shift = len(insert) - SLOT_LEN
    for f in backbone.features:
        s, e = int(f.location.start), int(f.location.end)
        if e <= SLOT_START:
            feats.append(f)
        elif s >= SLOT_START + SLOT_LEN:
            feats.append(
                SeqFeature(
                    FeatureLocation(s + shift, e + shift, f.location.strand),
                    type=f.type,
                    qualifiers=f.qualifiers,
                )
            )
    rec = SeqRecord(
        Seq(bb[:SLOT_START] + insert + bb[SLOT_START + SLOT_LEN :]),
        id="member",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    rec.features = feats
    return rec


def classify_hits(unit, orf_len: int):
    """Split one member's hits into (shared keyed by position, junction, intrinsic)."""
    orf = next(r for r in unit.unit.regions if r.label == "ORF")
    o_start, o_end = orf.start, orf.start + orf_len
    shared, junction, intrinsic = [], [], []
    for h in unit.hits:
        if h.start >= o_start and h.end <= o_end:
            intrinsic.append(h)
        elif h.end <= o_start:
            shared.append((("5'", h.start), h))
        elif h.start >= o_end:
            shared.append((("3'", h.start - o_end), h))
        else:
            junction.append(h)
    return shared, junction, intrinsic


def scan_barcodes(path: Path) -> dict:
    """Pool-level counts for every barcode in its handle/RT context."""
    left = plasmid_qc.FWD_PRIMER[-BARCODE_CONTEXT:]
    right = plasmid_qc.REV_PRIMER[:BARCODE_CONTEXT]
    with path.open() as fh:
        bcs = [r["barcode"].upper() for r in csv.DictReader(fh, delimiter="\t")]
    per_motif = Counter()
    flagged = Counter()
    lo, hi = len(left), len(left) + plasmid_qc.BARCODE_LEN
    for bc in bcs:
        seq = left + bc + right
        hits = [h for h in rbp_motifs.scan(seq) if h.end > lo and h.start < hi]
        for m in {h.motif for h in hits}:
            per_motif[m] += 1
        for m in {h.motif for h in hits if h.flagged}:
            flagged[m] += 1
    return {
        "pool": f"processed://{BARCODE_KEY}",
        "barcodes": len(bcs),
        "context": f"5' handle[-{BARCODE_CONTEXT}:] + barcode + gap+RT[:{BARCODE_CONTEXT}]",
        "barcodes_with_motif": {m: per_motif.get(m, 0) for m in rbp_motifs.MOTIF_NAMES},
        "barcodes_with_flagged_motif": {m: flagged.get(m, 0) for m in rbp_motifs.MOTIF_NAMES},
    }


def write(path: Path, rows: list[dict]) -> None:
    cols: list[str] = []
    for r in rows:
        cols += [c for c in r if c not in cols]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval="", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    constructs = json.loads(s3.processed.download(CONSTRUCTS_KEY).read_text())
    backbone = SeqIO.read(str(s3.processed.download(BACKBONE_KEY)), "genbank")
    members = active_l1()
    print(f"{len(members)} active L1 members")

    shared_seen: dict[tuple, list] = defaultdict(list)
    member_rows = []
    for r in members:
        orf = constructs[r["key"]]["construct"].upper()
        (uh,) = rbp_motifs.screen_record(member_record(backbone, orf))
        shared, junction, intrinsic = classify_hits(uh, len(orf))
        for key, h in shared:
            shared_seen[(key, h.motif)].append(h)
        row = {
            "candidate_id": r["candidate_id"],
            "design_id": r.get("design_id", ""),
            "key": r["key"],
            "orf_nt": len(orf),
        }
        counts = rbp_motifs.summarize(intrinsic)
        for m in rbp_motifs.MOTIF_NAMES:
            row[f"orf_{m}"] = counts[m]
        row["orf_flagged"] = sum(h.flagged for h in intrinsic)
        row["orf_PRE_per_kb"] = round(1000 * counts["PRE"] / len(orf), 3)
        row["orf_max_copy_number"] = max((h.copy_number for h in intrinsic), default=0)
        row["orf_tiled_motifs"] = ";".join(
            sorted({f"{h.motif}x{h.copy_number}" for h in intrinsic if h.copy_number >= 2})
        )
        row["junction_hits"] = ";".join(
            f"{h.motif}:{h.match}@{h.region}" + ("*" if h.flagged else "") for h in junction
        )
        row["junction_flagged"] = sum(h.flagged for h in junction)
        member_rows.append(row)

    n = len(members)
    shared_rows = []
    for ((side, off), motif), hs in sorted(
        shared_seen.items(), key=lambda x: (x[0][0][0] != "5'", x[0][0][1])
    ):
        h = hs[0]
        shared_rows.append(
            {
                "side_of_orf": side,
                "offset": off,  # 5': nt from transcript start; 3': nt after the ORF
                "motif": motif,
                "match": h.match,
                "region": h.region,
                "region_kind": h.region_kind,
                "designed": h.designed,
                "expected_pas": h.expected,
                "members_with_hit": len(hs),
                "members_flagged": sum(x.flagged for x in hs),
                "reasons": "; ".join(sorted({x.reason for x in hs if x.reason})),
            }
        )

    cassette_rows = []
    with CASSETTES.open() as fh:
        for c in csv.DictReader(fh, delimiter="\t"):
            hits = rbp_motifs.scan(c["cassette_nt"])
            counts = rbp_motifs.summarize(hits)
            cassette_rows.append(
                {
                    "candidate_id": c["candidate_id"],
                    "scheme": c["scheme"],
                    **{m: counts[m] for m in rbp_motifs.MOTIF_NAMES},
                    "flagged": ";".join(f"{h.motif}:{h.match}" for h in hits if h.flagged),
                }
            )

    barcodes = scan_barcodes(s3.processed.download(BARCODE_KEY))

    write(OUT_SHARED, shared_rows)
    write(OUT_MEMBERS, member_rows)
    write(OUT_CASSETTES, cassette_rows)

    def tot(col):
        return sum(r[col] for r in member_rows)

    full = [r for r in shared_rows if r["members_with_hit"] == n]
    summary = {
        "members": n,
        "inputs": {
            "designed_cds": f"processed://{CONSTRUCTS_KEY}",
            "backbone": f"processed://{BACKBONE_KEY}",
            "manifest": MANIFEST.name,
        },
        "transcript": "EFS end -> 3' LTR R end (AATAAA + 25 nt), insert as INSERT in this script",
        "decisions": {
            name: {"decided_by": d, "implemented_by": i, "reason": r}
            for name, (d, i, r) in {**rbp_motifs.DECISIONS, **DECISIONS}.items()
        },
        "shared": {
            "positions_in_every_member": len(full),
            "by_motif": dict(Counter(r["motif"] for r in full)),
            "flagged_in_designed_sequence": [
                f"{r['motif']} {r['match']} in {r['region']}"
                for r in full
                if r["designed"] and r["members_flagged"]
            ],
            "positions_in_some_members_only": len(shared_rows) - len(full),
        },
        "junction": {
            "members_with_junction_hit": sum(bool(r["junction_hits"]) for r in member_rows),
            "members_with_flagged_junction_hit": sum(
                r["junction_flagged"] > 0 for r in member_rows
            ),
        },
        "orf_intrinsic": {
            "hits_by_motif": {m: tot(f"orf_{m}") for m in rbp_motifs.MOTIF_NAMES},
            "members_with_motif": {
                m: sum(r[f"orf_{m}"] > 0 for r in member_rows) for m in rbp_motifs.MOTIF_NAMES
            },
            "members_with_PRE_ge2": sum(r["orf_PRE"] >= 2 for r in member_rows),
            "members_with_tiled_flagged_motif": sum(
                bool(r["orf_tiled_motifs"]) for r in member_rows
            ),
        },
        "cassettes": {
            "encodings": len(cassette_rows),
            "with_flagged_hit": sum(bool(r["flagged"]) for r in cassette_rows),
            "hits_by_motif": {m: sum(r[m] for r in cassette_rows) for m in rbp_motifs.MOTIF_NAMES},
        },
        "barcode_pool": barcodes,
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
