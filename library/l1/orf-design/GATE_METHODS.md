# How four of the L1 gene-list gates are actually decided

The L1 gene list is produced by the S1E3 cohort filter,
[`final_library.py`](../surface-filters/final_library.py), run over
accessible-surfaceome release 1.3.0. The brief states each gate as a rule. This page states
how the call behind the rule is made: which attribute or predictor decides it, what version,
over which sequence, at what threshold, and where the answer is durably kept. It covers the
four gates an external reviewer cannot otherwise check — anchoring, the co-receptor
requirement, the membrane-anchoring test inside the low-literature tier, and the removal of
EPPIN.

One distinction governs all four, so it is worth stating once. The filter reads **annotations
from the catalog record**; it does not run predictors. Every gate below is a query against
fields of a per-gene record in release 1.3.0 (Zenodo
[doi:10.5281/zenodo.23445521](https://doi.org/10.5281/zenodo.23445521), cut 2026-10-18),
pinned for this experiment as `external://accessible-surfaceome/release-1.3.0/deep_dives_1.3.0.tar.gz`,
SHA-256 `4ba846e8f10c98625cfab9557ffd8b9c2c63438c3d1d6c972251d71b4ac0e967`, with every input
hash recorded in [`library_pipeline_inputs.json`](library_pipeline_inputs.json). The
predictors that could independently check those annotations — DeepTMHMM for topology, NetGPI
for GPI anchors — run upstream in the catalog or downstream in the tag-site design, never
inside the gene gate. Where the two disagree, the gate follows the annotation. That is the
mechanism behind the EPPIN removal in section 4, and it is the single most useful thing a
reviewer can know about these gates.

A second distinction, relevant to three of the four: the gates run on **the gene's UniProt
canonical sequence**, because that is what the catalog record describes. The designed CDS a
library member actually carries is chosen later, by the proteoform ladder in
[`library_pipeline.py`](library_pipeline.py), and for a substantial minority of genes it is
not the canonical sequence. No gate below was evaluated on the designed CDS.

## 1. Anchoring: "at least one transmembrane helix, or GPI-anchored"

The operative rule is the `anchored()` function, and it is a two-term disjunction over two
different kinds of evidence:

```python
def anchored(r):
    f = r["filters"]
    return (
        f.get("subcategory") == "GPI_anchored"
        or ((r.get("canonical_topology") or {}).get("tm_helix_count") or 0) >= 1
    )
```

| Term | What decides it | Sequence | Threshold |
|---|---|---|---|
| transmembrane helix | `canonical_topology.tm_helix_count`, the helix count from DeepTMHMM 1.0.24 (the v1 release; `topology.TOOL_V1`, `deeptmhmm_v1`) as recorded in the catalog's topology table | the gene's UniProt canonical sequence | `>= 1` helix |
| GPI anchor | `filters.subcategory == "GPI_anchored"` — the catalog's **architecture annotation**, one of seven mutually exclusive values (`single_pass_T1`, `single_pass_T2`, `multi_pass`, `GPCR`, `GPI_anchored`, `tetraspanin`, `other`) | none: a categorical claim about the gene, not a prediction over a sequence | none |

So **no predictor makes the GPI call in this gate, and there is no score or threshold behind
it.** The architecture annotation originates as `executive_summary.subcategory`, a field the
catalog's per-gene deep-dive agent emits, and is then copied verbatim into `filters.subcategory`
so the catalog can be indexed on it. It is a language-model judgment with cited evidence in the
record, reviewed in aggregate rather than per gene, and it is not the output of NetGPI, of a
UniProt keyword query, or of any sequence-based rule.

NetGPI 1.1 does enter L1, but downstream and for a different purpose. A GPI anchor matters for
tag placement because the ω-site deletes the C-terminus, so an ecto-C tag on a GPI-anchored
member is removed with the propeptide and reads as a non-expressor. The GPI set used for
placement is built by [`build_gpi_set.py`](build_gpi_set.py) into [`gpi_set.tsv`](gpi_set.tsv),
and it is a union of two sources, neither of which is the gene gate's architecture annotation:

| | Source | Detail |
|---|---|---|
| Predictor | NetGPI 1.1, run through the DTU Health Tech web service in `short` output format | 918 predictions over 148 genes: one reference form per gene plus 770 alternative forms (UniProt isoforms and Ensembl transcripts) across 126 of them. Raw output committed at [`netgpi_gpi_human_predictions.txt`](../../../library/l1/proteoform-selection/netgpi_gpi_human_predictions.txt) |
| Curation | the `gpi_anchored` triage reason, exported from the catalog's `triage_run_public` table | 130 genes, committed at [`gpi_anchored_genes.tsv`](gpi_anchored_genes.tsv) on 2026-10-15 |

The rule is `anchored = NetGPI calls any form anchored, OR curation calls the gene anchored`.
NetGPI's own binary call is taken as it comes, with no score cutoff applied to a positive. A
score cutoff applies only to **negatives**, and only to the reference form's negative:
`CONFIDENT = 0.55`, below which a NetGPI negative is treated as uninformative and curation
stands. The value is read off the control run rather than tuned — every genuine false negative
on 148 curated-GPI plus 100 single-pass reference sequences fell below it (CD14 0.388,
FCGR3B 0.445, ART4 0.534, NYX 0.549), and CD14 and FCGR3B are textbook GPI-anchored proteins.
On that control set the predictor reads sensitivity 89.9%, specificity 99.0%, MCC 0.875,
against the 0.895 the NetGPI paper reports. Applying the cutoff to the strongest negative
across all forms instead of to the reference form's reads an alternative form's confident
negative as the reference's, which is wrong because alternative forms are different proteins
that are legitimately unanchored.

The resulting set is 157 GPI-anchored genes of 159 evaluated: 138 called by NetGPI, 15 held by
curation under an uninformative negative, and 4 (CFC1, CRIPTO, GML, RAET1G) where curation and
a confident NetGPI negative conflict. Those 4 are flagged for review rather than silently
resolved; all 4 sit on the N-after-signal-peptide rung, so no placement turns on them.

**Provenance: durable.** The architecture annotation that decides the gate is in the pinned
1.3.0 deposit; the NetGPI predictions and the curated list are committed in this repository,
and `gpi_set.tsv` rebuilds from them byte for byte under `library_pipeline.py`'s verifier. The
NetGPI 1.1 installation itself is not preserved — the trained weights are not distributed and
the DTU service is the only way to run the published tool, so a rerun depends on that service
remaining available. Whether the predictor packages may be kept in a private S3 prefix is an
open licence question, recorded in the cover note.

## 2. The co-receptor gate: obligate partners, and the TPM that excuses them

The rule is `coreceptor_ok()`, and it has four layers, evaluated in this order:

| Order | Layer | Effect |
|---:|---|---|
| 1 | six hand-reviewed accessions in [`coreceptor_rules.py`](../surface-filters/coreceptor_rules.py) | decides outright, pass or fail, and is the only layer that distinguishes alternative partners from simultaneously required ones |
| 2 | the `INCLUDE_ALONE` allowlist in `coreceptor_calls.tsv` | passes a subunit whose saved deep-dive rationale shows it reaches the surface without a heterologous partner |
| 3 | `filters.co_receptor_dependency != "required"` | passes: the catalog does not consider the gene partner-dependent at all |
| 4 | the partner-expression test | passes if partners are measurably common |

**Where the obligate-partner relationships come from.** They are not curated in this
repository. Both the dependency flag (`filters.co_receptor_dependency`, which reads `required`
for 333 of the 5,130 cohort genes) and the partner list itself
(`accessibility_risks.co_receptor_requirements.partners`) are fields of the catalog's per-gene
record, produced by the same deep-dive agent that emits the architecture annotation, with
cited evidence in the record. They are language-model-curated catalog attributes, and their
free text shows it: some partner entries are phrases rather than symbols
(`"ATP4B-like beta subunit (ATP1B1/FXYD family)"`), which no identifier join resolves. Hand
curation enters only as the two override layers above — six accessions reviewed against the
literature (ATP1B3, CDC50A, P2RX6, HFE, ATP8B1, CDC50B; four of them also carrying a reviewed
exception to the complex-partner localization veto), plus the `INCLUDE_ALONE` allowlist.

**Which HPA product supplies the TPM.** The per-gene median across HPA cancer cell lines, from
the HPA RNA cell-line consensus file `rna_celline.tsv`
(`proteinatlas.org/download/tsv/rna_celline.tsv.zip`; 20,151 genes × 1,206 cell lines, values
in TPM). The per-gene median was computed in this repository on 2026-08-04 by
[`build_expression_reference.py`](../surface-filters/pipeline/build_expression_reference.py),
which streams the full per-cell-line matrix — the published per-gene summary carries only a
maximum, and a maximum answers "detected in at least one line," not "a standard cell line
supplies it." The gate reads the result as `median_tpm_celline` from
`external://expression-tpm/rederived-2026-08-04/hpa_celline_per_gene.csv`, SHA-256
`ba9735fa00dc31fb793e86af494cc63ed1054315ce9de2297c8a6740782aeeff`. The threshold is median
`>= 5` TPM (`PARTNER_MIN_TPM`), chosen as the point at which a partner is present in roughly
half the cell lines rather than merely detectable somewhere.

Two qualifications the brief's wording does not carry, both of which loosen the rule:

- **Partners whose median cannot be measured are dropped, not failed.** The test is
  `bool(meds) and all(m >= 5 for m in meds)` over partners with a measurable median, so a gene
  passes when at least one partner resolves and every resolving partner clears 5. A partner
  given as free text, or absent from the HPA matrix, is silently excluded from the test. This
  is the practical consequence of the free-text partner entries above.
- **Self partners are removed first.** A partner naming the gene itself, or containing "homo"
  or "self", is a homo-oligomer requirement the gene satisfies by being expressed, so a gene
  whose only partners are self passes (BDKRB1, OTOP1/2/3, PCDHGB7).

The general partner lookup resolves a partner through the gene resolver to an Ensembl base ID
and falls back to a display-symbol join when no Ensembl ID matches. The six reviewed routes use
a stricter lookup that requires exactly one matching HPA matrix gene and never falls back to a
symbol.

**Provenance: durable for the data, with one unrecorded upstream version.** The partner lists
and dependency flags are in the pinned 1.3.0 deposit; the reviewed rules are committed;
`coreceptor_calls.tsv` and the saved rationales are pinned with SHA-256 on
`processed://sprint4/s4e9-orf-design-fg-panel/library-design/pipeline-inputs/2026-11-02/s1e3/`,
copied from `processed://sprint1/s1e3-surface-library-filters/partner-audit-2026-09-10/inputs/`.
The gap is the HPA release, and it is not recoverable. The download URL is unversioned and the
re-derivation records only the fetch date, so the exact HPA release behind these medians is
identified by `2026-08-04` and the derived file's SHA-256, not by a version number. The
`download_traceability.json` that `build_expression_reference.py` writes alongside the table —
which holds the per-source URL, SHA-256, byte count and timestamp that would name the release —
is not committed, is not on the `external://expression-tpm/` prefix (which holds only the three
derived tables), and is not on local disk. Re-running the download would fetch whatever HPA
serves today, not the release used. So the medians are pinned and reproducible as values, and
the release they were computed from can now be established only by matching them against HPA's
own release archive. The practical consequence is bounded: this is one threshold on a partner's
breadth of expression, not a per-gene measurement, and the three other gates that use expression
read the separately pinned genome-wide table.

## 3. "Membrane-anchored" in the low-literature UniProt tier

The low-literature UniProt surface tier admits genes whose catalog surface call rests on
literature the catalog has little of, letting UniProt's own surface annotation stand in. Becca's
question is what "membrane-anchored" means there, and the answer is the shortest in this
document: **it is the same `anchored()` function as section 1, called unchanged.** One
transmembrane helix or more on the canonical sequence by DeepTMHMM 1.0.24, or the catalog
architecture annotation `GPI_anchored`. There is no separate definition, no additional
predictor, and no threshold beyond the `>= 1` helix count.

The full tier test, as the code implements it, is five conditions:

| Condition | As implemented |
|---|---|
| UniProt annotates the protein at the cell surface | bit 0 of the catalog row's `db` mask; equivalently the release's `low_lit_uniprot` badge, which is this bit and the paper count below taken together |
| little literature | `filters.n_papers_found < 100`, and the count must be present |
| not already canonical-tier | `not passes_canonical(filters)` |
| membrane-anchored | `anchored(r)` — identical to section 1 |
| expressed, and synthesizable | maximum TPM `>= 5` (`EXPR_MIN_TPM`) across the four expression references, and ORF `<= 5,000` nt (`MAX_ORF_NT`) |

Two things worth stating because the tier is easy to read as a weaker version of the core
filter rather than a different filter. First, a gene already admitted to the core, or admitted
by the screen-validated GPCR tier, is skipped before these tests run, so the tier cannot
double-count. Second, the tier applies **only** those five conditions: the ectodomain-size,
co-receptor, localization-reason, accessibility and specificity gates of the core filter are
not evaluated. A tier member can therefore carry an ectodomain below the 30-residue minimum
that the core filter imposes.

**Provenance: durable.** Every field is in the pinned 1.3.0 deposit, and the deposit publishes
the `low_lit_uniprot` badge directly, so a reviewer can reproduce the tier membership from the
deposit and the committed filter alone. A stale comment above this block in `final_library.py`
still reads "max TPM >= 10"; the code uses `EXPR_MIN_TPM`, which is 5, and 5 is what the brief
states.

## 4. Why EPPIN is removed

EPPIN passed the anchoring gate on the GPI term, and on nothing else. Read directly out of the
pinned deposit, its record carries `filters.subcategory = GPI_anchored`,
`executive_summary.subcategory = GPI_anchored` (the field the former is copied from) and
`canonical_topology.tm_helix_count = 0`. Its canonical sequence is 133 residues with a signal
peptide and **zero** transmembrane helices, so the helix term of `anchored()` is false and the
architecture annotation is the only reason it entered the library. It clears every other gate on
its own — `co_receptor_dependency` is `none` and `surface_specificity` is `surface_dominant` —
so no other gate would have caught it. That is what the exclusion record describes,
and it is correct: the architecture annotation, as section 1 establishes, is a
language-model-emitted field of the catalog record, carried into `filters.subcategory` for
indexing. It is the gate's sole authority on GPI anchoring, and here it is the sole source that
calls EPPIN anchored at all.

Four independent sources disagree with it:

| Source | On EPPIN | Detail |
|---|---|---|
| UniProt | not GPI-anchored | no KW-0336 keyword on O95925 |
| NetGPI 1.1 | not GPI-anchored | all five human forms, likelihood 0.994–0.996; reference form O95925 at 0.996 — far above the 0.55 cutoff, so this is a confident negative, not an uninformative one |
| DeepTMHMM | no anchor | signal peptide, 0 transmembrane helices, 112 extracellular residues of 133 |
| the catalog's own GPI curation | not GPI-anchored | EPPIN is absent from the 130-gene `gpi_anchored` triage export, so the curated attribute and the architecture annotation disagree with each other |

That last row is the sharpest point for a reviewer: the catalog carries two GPI-relevant
attributes, and only the architecture annotation calls EPPIN anchored. [`gpi_set.tsv`](gpi_set.tsv)
records EPPIN as `gpi_anchored = 0`, basis `netgpi_not_anchored`, `in_curated_set = 0` — one of
only two genes in that table the union rule rejects.

EPPIN does reach the sperm surface. It is a component of the eppin protein complex with
lactotransferrin, clusterin and semenogelin, isolated from human semen and ejaculate spermatozoa
([Wang et al., 2007](https://pubmed.ncbi.nlm.nih.gov/17567961/)); that work proposes the
partners' receptors, rather than eppin itself, as the complex's plasma-membrane attachment, and
does not demonstrate an anchor on eppin. Taken with four sources that find no anchor on its own
sequence, this is consistent with EPPIN carrying none. A surface-display library needs the
protein to carry its own anchor; a member with none would read as a non-expressor for a reason
that has nothing to do with the tag. It is therefore removed outright rather than flagged, recorded with
its date and evidence in [`gene_list_exclusions.tsv`](gene_list_exclusions.tsv).

**Provenance: durable, and verified end to end.** The NetGPI output for all five EPPIN forms is
committed and checkable line by line; the DeepTMHMM topology is both in the pinned deposit and
in this repository's committed canonical-topology table; and the architecture annotation EPPIN
was admitted on was read out of `deep_dives_1.3.0.tar.gz` whose SHA-256 matches the value
recorded in `library_pipeline_inputs.json`, so every step of this removal can be reproduced from
pinned artifacts. One item is unestablished, and it does not affect the removal: on the catalog
side, the append-only `netgpi_prediction` and `gpi_curation_record` tables that would hold these
predictions with their `tool_version`, source URI and source checksum exist in the schema but
are not part of the schema pinned to the live public database, so whether they are populated is
not established here; a row count against both databases would settle it. The NetGPI 1.1 install
record is covered by the same open licence question as section 1.

## What a reviewer should take from this

| Gate | Decided by | Independent check available |
|---|---|---|
| transmembrane helix | DeepTMHMM 1.0.24 on the canonical sequence, via the catalog | yes — topology is published per gene in the deposit |
| GPI anchor (gene gate) | the catalog's architecture annotation; no predictor, no threshold | only indirectly, by comparing against the committed NetGPI set and the curated triage export |
| obligate partner | the catalog's partner list and dependency flag, with six hand-reviewed overrides | partially — partner lists are published; the free-text entries resolve to nothing |
| partner ubiquity | median HPA cell-line TPM `>= 5`, re-derived 2026-08-04 | yes for the values, by SHA-256; the HPA release they came from is not recoverable |
| membrane-anchored (low-literature tier) | the same two-term rule as the anchoring gate | yes |
