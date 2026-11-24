# Human references whose chosen form changes the tag lane

**Question.** For which L1 human library references is the problem the human proteoform we
chose, rather than how a member was drawn? The lane check on the "N at residue 1, multipass"
rung raised the question through SLC22A17, where the cyno ortholog matches the classic
receptor but the human reference does not.

**Method.** [`screen_reference_forms.py`](screen_reference_forms.py) compares each of the
2,360 active human `library_reference` members with its gene's reviewed UniProt canonical
(release 2026_04). It asks two questions:

- **Lane.** Each protein gets the placement lane `l1_tag_paths` would give it on its own
  DeepTMHMM v1 string (by sequence hash). Does the member's lane differ from the canonical's?
- **Start codon.** Does the member's CDS start on a non-ATG codon?

Per-member results, with UniProt entry and sequence versions, are in
[`reference_form_screen.tsv`](reference_form_screen.tsv).

## Result

2,147 references are identical to their canonical, and 152 more share its first 20 residues.
18 are flagged: 15 whose lane differs from the canonical's, and 3 with a non-ATG start. They
fall into five groups.

| Group | Members | What is wrong with the human form | Status elsewhere |
|---|---|---|---|
| A. Canonical signal peptide dropped | ADGRF3 `l1-p0067`, CACHD1 `l1-p0215`, CD200R1L `l1-p0272`, LRRC37B `l1-p1136`, LMLN `l1-p1106` | The form starts downstream of, or diverges from, the canonical's signal peptide. Its ectodomain then has to cross the membrane with no leader, the same defect as the PTCH1 and SLC9D1 isoforms. | Lane check (PR #296) recommends replacing ADGRF3 with Q8IZF5-1. The type III review (PR #282) recommends swapping CACHD1, CD200R1L and LRRC37B, and leaves LMLN undecided. |
| B. Sequence added upstream of the signal peptide | EPHA6 `l1-p0578` (+94 aa, Ensembl canonical EPHA6-201), RGMB `l1-p1731` (+41 aa, RGMB-201, not the Ensembl canonical) | The canonical's signal peptide becomes internal. v1 then reads it as a membrane helix with a cytoplasmic N-terminus, so the member's lane and topology are not the protein's. | Not flagged before. EPHA6 is in the TAG_PATHS reconciliation as "ordered proteoform not canonical"; RGMB only as a GPI reroute. |
| C. Non-ATG start | SLC22A17 `l1-p1882` (CTG), TRPV6 `l1-p2281` (ACG), OR8G3 `l1-p1525` (GTG) | The CDS opens on a codon that initiates poorly in a construct. Translation is then expected to begin mostly at a downstream AUG, which is a different protein. | `library_cds_boundary_decisions.tsv`: SLC22A17 and TRPV6 "replace start with ATG"; OR8G3 "prepend ATG". None is approved for canonical activation yet. |
| D. Different TM region, already reconciled | SERINC5 `l1-p1799`, SLC17A2 `l1-p1861`, TSPAN10 `l1-p2284`, LSR `l1-p1160`, AQP12B `l1-p0128` | Each is lane-assigned on its own topology on purpose (TAG_PATHS, "Reconciling against the 515"). | AQP12B is on hold in the lane check. |
| E. Member gains or loses a signal-peptide call | ADGRD2 `l1-p0061`, TMEM132B `l1-p2168`, OR10R2 `l1-p1374` | ADGRD2 and TMEM132B gain a signal peptide the canonical lacks. OR10R2's canonical gets an unusual v1 signal-peptide call for an olfactory receptor. In each case the member's lane is at least as plausible as the canonical's. | Low priority. |

## SLC22A17: the start codon is the problem, not the form

The human reference is Ensembl's canonical transcript SLC22A17-202 (ENST00000354772, 631 aa),
which UniProt describes as isoform 2 (NgalR-2) of Q8WUG5. It differs from the displayed
649-aa isoform 1 only at residues 499–517.

- **UniProt treats the long N-terminus as real.** It records five cDNAs as "erroneous
  initiation: truncated N-terminus", so the 538-aa receptor that starts at Met112 is the
  truncation, not the reference.
- **The defect is the CTG start.** The record opens on CTG, so in a lentiviral construct most
  ribosomes would be expected to skip it and initiate at the downstream AUG. That would make
  the classic receptor, the form the cyno ortholog already matches.

The fix is already proposed: replace the CTG with ATG, without adding a residue
(`library_cds_boundary_decisions.tsv`, status "annotation supported, pending review"). The
cyno mismatch then needs no action, because the Compara 116 adoption (PR #283) retires that
ortholog.

## Decisions

Becca approved all four recommendations on 2026-11-12. Each flagged row in the table carries
the decision in `decision` and `decided_by`, and they are applied on PR #283.

| Decision | Members | Recommendation |
|---|---|---|
| Approve the ATG replacements | SLC22A17, TRPV6 | Approve: UniProt supports both full-length forms, and TRPV6's ACG initiation is experimentally supported |
| OR8G3 start | OR8G3 | Curation is still required before choosing between prepending ATG and recoding GTG to ATG |
| Replace references that drop the signal peptide | ADGRF3, CACHD1, CD200R1L, LRRC37B; LMLN undecided | Replace each with the form that starts at the canonical signal peptide, via named rows in `reference_overrides.tsv` |
| Replace references with an upstream extension | EPHA6, RGMB | Replace each with the form that starts at the signal peptide (UniProt Q9UF33 and Q6NW40); recheck the lanes afterward |

Every change would carry `decided_by` / `implemented_by` / `reason`. Reference overrides
belong to the reviewer-package session (PR #283), which owns the reference ladder.

## Reproducing

    uv run --locked python screen_reference_forms.py            # offline, from the snapshot
    uv run --locked python screen_reference_forms.py --refresh  # re-fetch UniProt

UniProt entries for the library genes are snapshotted in `uniprot_reference_screen.tsv`, and
the release in `uniprot_reference_screen_snapshot.json`. Topology is DeepTMHMM v1 from the
snapshot pinned in `topology.SNAPSHOT_KEY`.

## Decision log

| Date | Decision | Reason | Decided by / implemented by |
|---|---|---|---|
| 2026-11-12 | Screen every human reference against its UniProt canonical by DeepTMHMM v1 lane and start codon | Find where the chosen human form, not the draw, is the problem | Becca ("look at where the human form is the problem") / Claude Code session l1-nterm-res1-multipass, PR #280 |
| 2026-11-12 | Fixes approved: ATG starts for SLC22A17 and TRPV6; OR8G3 stays in curation; ADGRF3, CACHD1, CD200R1L, LRRC37B, EPHA6 and RGMB replaced by the form that starts at the canonical signal peptide; LMLN undecided | The recommendations in "Decisions for Becca" | Becca ("Send all fixes") / applied by Claude Code session l1-reviewer-package, PR #283; recorded here by l1-nterm-res1-multipass, PR #280 |
