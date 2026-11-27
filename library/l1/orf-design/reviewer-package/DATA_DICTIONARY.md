# Data dictionary: `l1_order_catalog.tsv.gz`

`l1_order_catalog.tsv.gz` is a single flat, plate-ordered table of all **3,352** ordered
constructs (across **2,730** built CDSs and **35** 96-well plates) for the Deliverome v1
surface-protein ORF library (order of record `2026-11-25-pre-freeze-e`). It is generated
deterministically by [`build_l1_order_catalog.py`](build_l1_order_catalog.py) and every DNA
sequence is verified against `synthesis_compliance.tsv` (`recoded_insert_sha256`) and
`l1_order_inserts_recoded.fasta`.

Unpack with `gunzip -k l1_order_catalog.tsv.gz` (or open directly in Python with
`pandas.read_csv("l1_order_catalog.tsv.gz", sep="\t")`).

| Column | Type | Description |
|---|---|---|
| `insert_id` | string | Plate-ordered insert identifier (`l1-i0001`..`l1-i3352`), sorted shortest to longest. |
| `plate` | integer | 96-well synthesis stock plate (`1`..`35`). |
| `well` | string | Column-major 96-well coordinate (`A01`..`H12`). |
| `arm_id` | string | Stable construct arm identifier (`<design_id>-a1`, `-a2`, `-a3`, or `-sn1` for snorkels). |
| `design_id` | string | Stable member/CDS identifier (`l1-pNNNN`). Multiple arms of an internally tagged CDS share one `design_id`. |
| `gene_symbol` | string | HGNC human gene symbol (shared by human isoforms and mouse/cyno orthologs of that gene). |
| `candidate_set` | enum | `library_reference` (human reference proteoform), `isoform` (human alternate isoform), or `ortholog` (mouse/cyno). |
| `species` | enum | `human`, `mouse`, or `cyno` (*Macaca fascicularis*). |
| `source_key` | string | Candidate key linking back to upstream selection tables (`lib:<GENE>`, `iso:<GENE>:<ID>`, `orth:<SPECIES>:<GENE>:<ID>`). |
| `transcript_or_accession` | string | Source Ensembl transcript ID (`ENST`/`ENSMUST`/`ENSMFAT`) or UniProt accession used to build the CDS. |
| `also_encodes` | string | Additional identical-paralog human gene symbol(s) whose reference protein is 100% identical to this CDS (4 deduplicated paralog pairs). |
| `rung` | enum | Tag-placement topology rung: `nterm` (N after SP), `multipass` (N-out >1 TM), `type3` (single-pass N-out no SP), `cterm` (extracellular C-terminus), or `internal` (extracellular loop / snorkel). |
| `placement` | string | Specific placement rule applied (e.g., `native_site`, `igk_swap`, `multipass_igk`, `type3_tag_alone`, `cterm`, `internal_loop`, `snorkel`). |
| `insert_after_residue` | string | 1-based native protein residue after which the 51-aa `GGGGS–SpyTag003–GGGGSGGGGS–P–ALFA–P–GGGGS` cassette is inserted (`1` for residue-1 N-terminal tags; `C-term_snorkel` for appended 106-aa snorkel modules). |
| `n_tm` | integer | Transmembrane helix count predicted by DeepTMHMM 1.0.24 on the exact built native protein. |
| `signal_peptide_aa` | integer | Signal peptide length in residues predicted by DeepTMHMM 1.0.24 (`0` if none). |
| `native_aa` | integer | Native untagged protein length (amino acids, excluding stop). |
| `tagged_aa` | integer | Ordered tagged protein length (amino acids, translated from `insert_dna`, including Igκ leader or snorkel where used). |
| `insert_nt` | integer | Total synthesized DNA insert length in nucleotides (Esp3I shuttle sites, Kozak, tagged ORF, stop+CC pad, T7, JT832 handle, 20-N barcode placeholder, padlock, RT T9A, 10x CS1). |
| `vector_bp` | integer | Total assembled lentiviral plasmid length in base pairs inside destination vector `p0039` (`EFS–insert–IRES2–PuroR–WPRE`). |
| `synthesis_category` | enum | Vendor-calibrated synthesis complexity tier after synonymous recoding (`Low`, `Medium`, `High`; `0` constructs are `Not Accepted`). |
| `codons_recoded` | integer | Number of synonymous codons edited in `synthesis_compliance.py` to clear GC, repeat, or Tn5-ME synthesis thresholds (protein, cassette, and readout block unchanged). |
| `refold_verdict` | string | ColabFold (`alphafold2_ptm`) tag-refold check verdict (`PASS` or flagged/dispositioned; all 3,352 constructs refolded). |
| `refold_mean_plddt` | float | Mean per-residue pLDDT of the ColabFold-predicted tagged protein structure. |
| `shed_risk` | string | Extracellular proteolytic shedding annotation from UniProt/literature curation (`none` or risk label). |
| `decided_by` | string | Provenance label naming the rule, session, or human decision maker for this construct's placement. |
| `reason` | string | Human-readable rationale for the construct's placement, leader swap, or disposition. |
| `insert_sha256` | hex | SHA-256 digest of `insert_dna`, matching `synthesis_compliance.tsv` (`recoded_insert_sha256`). |
| `tagged_protein_aa` | string | Complete amino acid sequence encoded by `insert_dna`. |
| `insert_dna` | string | Complete synthesis-ready DNA insert sequence (uppercase A/C/G/T, with a 20-`N` placeholder at the barcode slot assigned at order synthesis). |
