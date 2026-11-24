# S1E3 signal-peptide prioritization: frozen 2026-09-11 review

The exact main-library identities behind the **519 precedent-informed** and
**1,667 broad topology-based** candidates are reproducible from
`signal_peptide_review.py`. The script contains the classification rules,
gene-specific qualifications and literature URLs. No model calls or live catalog
updates are used during reproduction.

The scientific writeup is saved on [S3E1: GPCR surface expression / SP / BRIL](https://www.notion.so/3d4ab86c768381c1b4c9fc660927be16).
Library context: [S1E3](https://www.notion.so/3b8ab86c768381bc8a73dc79ec72d832).
This is AI-authored prioritization awaiting scientific review, not a validated
rescue list or an instruction to modify every candidate.

## Frozen results

| SP-confidence tier | Meaning | Main proteins | Tedman below 1,000 |
| --- | --- | ---: | ---: |
| A | Direct improvement precedent, with host/species/co-modification caveats | 4 | 1 |
| B | Family, construct-feasibility or translocation precedent | 515 | 223 |
| C | Mechanistic/topology inference | 1,148 | 28 |
| D | Hold for topology-specific design or annotation review | 677 | 1 |
| Total | Main proteins only | 2,344 | 253 |

A+B = 519; A+B+C = 1,667. These SP tiers are independent of the existing library
tiers (T1/T2/T3/LL/SV/OE). The broad pool contains 1,121 annotated native-SP
replacement hypotheses and 546 addition hypotheses. Counts express the output
of the stated rules, not a measured probability of benefit.

GLP1R is B, replacement, Tedman 560.3241 (SD 52.3678). FZD1 remains D because the
cached N-terminal topology conflicts with the GPCR family assignment; it is not
a biological rejection. MRAP/MRAP2 are held because of dual topology. Most B
assignments inherit family evidence. Native-SP compatibility does not establish
that native targeting is defective. All designs need selected-ORF/cleavage-site
review before synthesis.

## Durable data and complete membership lists

All files live under this fixed S3 prefix:

`s3://deliverome-processed/sprint1/s1e3-surface-library-filters/signal-peptide-review-2026-09-11-v1/`

| File under prefix | Contents |
| --- | --- |
| `inputs/library_membership.tsv` | Exact 2,344-gene input membership |
| `inputs/cohort.jsonl.gz` | Frozen full catalog annotations used by the review |
| `inputs/tedman-v1.xlsx` | Original Tedman v1 Canonical worksheet source |
| `outputs/signal_peptide_all.tsv` | All 2,344 identities, tiers, design hypotheses, Tedman values/SD, flags, source rows, rationale and evidence URLs |
| `outputs/signal_peptide_precedent.tsv` | Exact 519-member A+B set with all annotations |
| `outputs/signal_peptide_broad.tsv` | Exact 1,667-member A+B+C set with all annotations |
| `outputs/signal_peptide_low_tedman.tsv` | All 253 below-1,000 genes with all annotations |
| `outputs/s1e3-signal-peptide-candidates.json` | Full machine-readable results and class summary |
| `outputs/signal_peptide_summary.json` | Mutually exclusive class counts |
| `archive_manifest.json` | SHA-256, size and S3 key for each archived input/output/document |
| `s1e3-signal-peptide-review-2026-09-11-v1.zip` | Portable inputs, outputs, script, metadata and scientific report; also attached to Notion |

Data files are deliberately not committed: S3 is the repository's data home.
The committed `signal_peptide_snapshot.json` pins all input keys and SHA-256
checksums. New reviews must use a new version prefix; do not update this snapshot
in place. The ZIP and individual membership tables are attached to Notion, so
the user can retrieve them without resolving a local computer path.

## Reproduce

From the repository root, with dependencies installed and S3 access configured:

```bash
uv run python library/l1/surface-filters/signal_peptide_review.py \
  --fetch --input-dir /tmp/s1e3-sp-inputs --output-dir /tmp/s1e3-sp-results
```

For an offline run, extract the Notion ZIP, point `--input-dir` at its `inputs`
directory, and omit `--fetch`. Outputs go to the explicitly supplied directory.
The script refuses input checksum mismatches and checks the frozen counts.
Blank TSV measurements mean **not measured**, not zero/negative.

The raw Tedman gene-label join gives 524 measured/251 low. Three blank gene IDs
are restored from explicit construct names (GPR18, XCR1, GPR35), and the
[verified GPR182 to ACKR5 rename](https://www.ncbi.nlm.nih.gov/gene/11318) gives
528 measured/253 low. All other 1,816 main genes are unmeasured. The strict
1,000 threshold is a requested prioritization cutoff, separate from S3E1's older
approximately 1,025 cutoff for its Lv comparison. Intensities are relative
surface immunostaining, not copies per cell. This snapshot measures canonical
gene-level priors, not all additional isoforms and orthologs.

## Validation

- Reproduced all 2,344 prior gene records exactly, including classification,
  design, evidence notes, intensity, SD and original worksheet row.
- Verified exact subset sizes 519 / 1,667 / 253 and the 528 measured total.
- Checksums pin the original membership, compressed cohort and source workbook.
- Archive publication verifies each object by downloading it and comparing
  SHA-256. Notion attachment placement and complete list counts are checked
  after saving.

No library membership or trafficking classifications were changed by archiving
this review. The confidence tiers remain provisional scientific judgments.
