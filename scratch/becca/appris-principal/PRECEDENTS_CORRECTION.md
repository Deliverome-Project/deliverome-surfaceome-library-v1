# Correction owed to the s1e3 precedents doc

`analysis/sprint1/s1e3-surface-library-filters/oversized_truncation_precedents.md`,
plexin section, currently reads:

> **Modern, multi-family:** PLXNA1ΔECD + PLXNB2ΔECD both surface-express + signal
> (Pinskey/Giger 2022, eLife, PMID 36169302) — the trick transfers across A- and B-class.

**The PLXNB2 half is wrong.** Verified here against the PMC full text (PMC9553217,
fetched from NCBI, superscripts preserved): every ΔECD construct in that paper is
PLXNA1 — `PLXNA1ΔECD` appears 8 times plus six tagged variants. PLXNB2 appears only
as full-length and once as `PLXNB2R1`, the GAP-dead arginine point mutant. There is
no PLXNB2ΔECD in it.

Suggested replacement:

> **Modern:** PLXNA1ΔECD surface-expresses + signals (Pinskey/Giger 2022, eLife,
> PMID 36169302) — **this is the only ΔECD plexin in any gene whose surface
> localisation has actually been measured.** The B-class ΔECD that does exist is
> Huang/Friedel 2021 (PMID 33514835), a human lentiviral ∆Ecto deleting aa 30–1189,
> Addgene #86238 — functional, but its surface expression was never measured either.

**Why this matters beyond the citation:** it means PLXNA1 is the single gene in the
whole plexin arm with a measured surface localisation for a ΔECD construct. Every
other gene's evidence is functional. A surface-stain QC arm is not optional.

**Not applied.** Editing that file pulls all of s1e3 into the pre-push S3 check, and
40 of its outputs are not on S3 — so landing this one line would mean publishing 40
files from an experiment this work does not own. Left for whoever owns s1e3 to apply,
or run `scripts/upload-experiment-outputs.sh analysis/sprint1/s1e3-surface-library-filters/`
first and then apply it.
