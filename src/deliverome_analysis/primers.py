"""Primer checks: NEB-equivalent Tm, secondary structure, binding sites, BLAST.

Why this exists
---------------
Checking a primer by hand means opening NEB's Tm calculator in a browser,
eyeballing IDT OligoAnalyzer for dimers, and running a BLAST — three tools, no
record of what was checked, and nothing a reviewer can rerun. This module does
all three from one call and returns a dataclass that serialises, so the check
lands in the experiment folder alongside the sequences it cleared.

Tm
--
``neb_tm`` reproduces NEB's Tm calculator offline. NEB computes Tm client-side
with the SantaLucia (1998) unified nearest-neighbour set (initiation
dH 0.2 kcal/mol, dS -5.7 cal/mol/K, symmetry dS -1.4, terminal-AT penalty)
followed by the Owczarzy (2008) salt correction, using a per-buffer monovalent
salt equivalent and a per-product primer concentration. Both tables are read
off NEB's own ``tmcalculatordata.json`` and reproduced in
``NEB_BUFFER_SALT_MM`` / ``NEB_PRODUCT_CT_NM``.

The Q5 buffer's 150 mM is the single parameter that matters most: running the
same nearest-neighbour model at a textbook 50 mM puts every Q5 primer about
11 degC low, which is large enough to change an annealing temperature decision.

``tests/test_primers.py`` pins this implementation against values read from
NEB's web calculator, so a drift in either direction fails CI rather than
quietly moving a Tm.

**Q5 only.** NEB selects a Tm method per product group and quarters the primer
concentration for Phusion, so method 4 is not a general answer -- applied to
Phusion it reads about 8 degC high. Products outside the validated Q5 family
raise rather than return an unchecked number; ``UNSUPPORTED_PRODUCTS`` says what
each would need.

Dimers and hairpins
-------------------
All secondary-structure energies come from primer3's ``ntthal`` engine via
``primer3-py`` -- a full thermodynamic alignment with gaps and loops, which is
the same calculation IDT's OligoAnalyzer reports. ``PRIMER3_CONDITIONS`` pins
the reaction to 150 mM monovalent, 0.2 mM dNTP and 500 nM primer so the numbers
correspond to a Q5 reaction; primer3's own defaults are 50 mM / 50 nM and read
several kcal/mol different.

``self_dimer`` and ``cross_dimer`` give the most stable duplex overall;
``end_stability`` gives the 3'-end configuration, which is the one that extends
into primer-dimer product rather than merely sequestering primer. Each returns
the structure's own melting temperature alongside dG, which is what decides
whether it is populated at the annealing temperature -- a -13 kcal/mol dimer
whose structure melts at 38 degC is not present during a 68 degC anneal.

Vector context
--------------
``vector_homology`` reports the longest stretch of a primer that recurs anywhere
in its own template. A second copy is a mispriming risk, and in a lentiviral
construct a long direct repeat also templates the recombination that swaps
barcodes between ORFs, so it matters for the vector and not only the PCR. It
runs automatically whenever ``check_primer`` is given a template.

``transcript_seed_sites`` reports conserved miRNA seed target sites, from the
TargetScanHuman table under s1e2. Seeds act on the transcript, so the scan is
only meaningful on the plus strand -- ``check_pair`` therefore scans the
amplicon span rather than the two primers, and a reverse primer's own sequence
is the wrong strand. Most sites in a real construct are pre-existing, so compare
a proposed edit against the same span before it.

Specificity
-----------
``blast_primers`` submits to NCBI and ``summarise_blast`` reduces the report to
3'-anchored matches, since an alignment that stops short of the 3' terminus
cannot be extended. ``genomic_amplicon_risk`` then does what Primer-BLAST does
and this module's per-primer view cannot: pairs up forward and reverse hits on
the same subject in the correct orientation and within a plausible distance, so
an off-target *product* is caught rather than just an off-target site.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import primer3

from deliverome_analysis.seqtools import gc_content, reverse_complement

# --------------------------------------------------------------------------- #
# NEB parameters, read off https://tmcalculator.neb.com/data/tmcalculatordata.json
# --------------------------------------------------------------------------- #

#: Monovalent salt equivalent (mM) per NEB buffer.
NEB_BUFFER_SALT_MM: dict[str, int] = {
    "standard_taq": 55,
    "thermopol": 40,
    "hemo_klentaq": 70,
    "crimson_taq": 55,
    "longamp": 100,
    "longamphs": 100,
    "multiplex": 90,
    "phusion_hf": 222,
    "phusion_gc": 222,
    "phusionflex_hf": 222,
    "phusionflex_gc": 222,
    "onetaq_gc": 80,
    "onetaq_std": 54,
    "q5": 150,
    "q5mm": 150,
    "q5u": 170,
    "q5bd": 150,
}

#: Product group -> (buffer key, recommended primer concentration in nM).
NEB_PRODUCTS: dict[str, tuple[str, int]] = {
    "Q5": ("q5", 500),
    "Q5 Hot Start": ("q5", 500),
    "Q5 Master Mix": ("q5mm", 500),
    "Q5 Blood Direct": ("q5bd", 500),
}

#: Product groups NEB supports that this module does not, and why.
#:
#: NEB does not use one formula for everything. Its calculator selects a method
#: per product group, and Phusion additionally quarters the primer
#: concentration::
#:
#:     case "Phusion": case "Phusion Hot Start Flex": s = 5; v /= 4;
#:     default: s = 4;
#:
#: This module implements method 4 -- SantaLucia 1998 nearest-neighbour with the
#: Owczarzy 2008 salt correction -- which is what the Q5 family uses. Applying it
#: to Phusion reads about 8 degC high: measured against NEB, not assumed. Rather
#: than ship numbers for products nobody has checked, these raise.
UNSUPPORTED_PRODUCTS: dict[str, str] = {
    "Phusion": "NEB uses method 5 and Ct/4; method 4 reads ~8 degC high",
    "Phusion GC": "as Phusion: method 5 and Ct/4",
    "Phusion Hot Start Flex": "as Phusion: method 5 and Ct/4",
    "Q5U Hot Start": "buffer q5u is 170 mM, not validated against NEB",
    "OneTaq": "method 4 but unvalidated; buffer onetaq_std 54 mM, Ct 200 nM",
    "OneTaq Hot Start": "unvalidated",
    "Taq": "unvalidated; buffer standard_taq 55 mM, Ct 200 nM",
    "Hot Start Taq": "unvalidated",
    "LongAmp Taq": "unvalidated; buffer longamp 100 mM, Ct 400 nM",
    "Master Mix": "method depends on whether the mix is Phusion-based",
    "Vent": "unvalidated",
    "Deep Vent": "unvalidated",
    "Hemo KlenTaq": "unvalidated",
    "EpiMark Hot Start": "unvalidated",
}

#: NEB refuses to recommend an annealing temperature above this for Q5/Phusion.
NEB_ANNEAL_CEILING_C = 72.0

# SantaLucia (1998) unified nearest-neighbour parameters.
# dH kcal/mol, dS cal/mol/K.
NN: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

INIT_DH, INIT_DS = 0.2, -5.7  # duplex initiation
TERM_AT_DH, TERM_AT_DS = 2.2, 6.9  # penalty per terminal A/T
SYMMETRY_DS = -1.4  # self-complementary duplexes only
R_CAL = 1.987  # gas constant, cal/mol/K (NEB uses 1.987)

_VALID = re.compile(r"^[ACGT]+$")


def _clean(seq: str) -> str:
    s = re.sub(r"\s", "", seq).upper().replace("U", "T")
    if not s:
        raise ValueError("empty sequence")
    if not _VALID.match(s):
        bad = sorted(set(s) - set("ACGT"))
        raise ValueError(f"sequence has non-ACGT characters: {''.join(bad)}")
    return s


def nn_thermodynamics(seq: str) -> tuple[float, float]:
    """Unified NN dH (kcal/mol) and dS (cal/mol/K) for a perfect duplex."""
    s = _clean(seq)
    dh, ds = INIT_DH, INIT_DS
    for a, b in itertools.pairwise(s):
        step_dh, step_ds = NN[a + b]
        dh += step_dh
        ds += step_ds
    for end in (s[0], s[-1]):
        if end in "AT":
            dh += TERM_AT_DH
            ds += TERM_AT_DS
    if s == reverse_complement(s):
        ds += SYMMETRY_DS
    return dh, ds


def _owczarzy_monovalent(fgc: float, salt_m: float) -> float:
    """Owczarzy (2008) reciprocal salt correction, monovalent-only branch."""
    ln_salt = math.log(salt_m)
    return 1e-5 * (4.29 * fgc - 3.95) * ln_salt + 9.4e-6 * ln_salt * ln_salt


def neb_tm(seq: str, product: str = "Q5", ct_nm: float | None = None) -> float:
    """Melting temperature in degC, matching NEB's Tm calculator.

    ``product`` selects the buffer salt and the recommended primer
    concentration; ``ct_nm`` overrides the latter.
    """
    s = _clean(seq)
    if product not in NEB_PRODUCTS:
        if product in UNSUPPORTED_PRODUCTS:
            raise NotImplementedError(
                f"{product!r} is not validated against NEB here: "
                f"{UNSUPPORTED_PRODUCTS[product]}. Validated products: "
                f"{', '.join(sorted(NEB_PRODUCTS))}."
            )
        raise KeyError(f"unknown product {product!r}; validated: {', '.join(sorted(NEB_PRODUCTS))}")
    buffer_key, default_ct = NEB_PRODUCTS[product]
    salt_m = NEB_BUFFER_SALT_MM[buffer_key] / 1000.0
    ct_m = (default_ct if ct_nm is None else ct_nm) * 1e-9

    dh, ds = nn_thermodynamics(s)
    # NEB uses the primer concentration directly (primerc = 1e-6 * ct), with no
    # Ct/4 divisor -- the textbook divisor assumes two strands at equal
    # concentration, whereas a PCR primer is in large excess over its template.
    # Using Ct/4 here puts every Tm about 2 degC low.
    tm_1m = (dh * 1000.0) / (ds + R_CAL * math.log(ct_m))
    inv = 1.0 / tm_1m + _owczarzy_monovalent(gc_content(s), salt_m)
    return 1.0 / inv - 273.15


def anneal_temp(fwd: str, rev: str, product: str = "Q5") -> tuple[float, bool]:
    """NEB's recommended annealing temperature, and whether it hit the ceiling.

    NEB uses Ta = Tm(lower) + 1 for Q5 and Phusion, clamped to 72 degC. A
    clamped value is a ceiling, not a match -- it leaves no room to raise
    stringency if the first reaction is dirty, which is worth knowing before
    ordering rather than after.
    """
    low = min(neb_tm(fwd, product), neb_tm(rev, product))
    ta = low + 1.0
    if ta > NEB_ANNEAL_CEILING_C:
        return NEB_ANNEAL_CEILING_C, True
    return ta, False


# --------------------------------------------------------------------------- #
# Secondary structure -- primer3 (ntthal) does all of it
# --------------------------------------------------------------------------- #

#: Reaction conditions handed to primer3. Matched to Q5 as closely as its model
#: allows: 150 mM monovalent (NEB's Q5 buffer equivalent), 0.2 mM dNTPs,
#: 500 nM primer. primer3's own defaults are 50 mM / 50 nM and will read several
#: kcal/mol different, so always pass these rather than relying on defaults.
PRIMER3_CONDITIONS: dict[str, float] = {
    "mv_conc": 150.0,
    "dv_conc": 0.0,
    "dntp_conc": 0.2,
    "dna_conc": 500.0,
}


@dataclass
class Structure:
    """A hairpin or dimer as primer3 reports it."""

    dg: float  # kcal/mol at 37 degC; 0.0 when no structure is found
    tm: float  # melting temperature of the structure itself, degC
    found: bool
    ascii: str = ""

    @property
    def concerning(self) -> bool:
        """IDT's convention: flag below about -6 kcal/mol."""
        return self.found and self.dg <= -6.0


def _structure(result: Any, want_ascii: bool = False) -> Structure:
    if not result.structure_found:
        return Structure(0.0, float("nan"), False, "")
    return Structure(
        dg=result.dg / 1000.0,
        tm=result.tm,
        found=True,
        ascii=(result.ascii_structure or "") if want_ascii else "",
    )


def hairpin(seq: str, **kwargs: Any) -> Structure:
    """Most stable hairpin, via primer3's thermodynamic alignment."""
    opts: dict[str, Any] = PRIMER3_CONDITIONS | kwargs
    return _structure(primer3.calc_hairpin(_clean(seq), output_structure=True, **opts), True)


def self_dimer(seq: str, **kwargs: Any) -> Structure:
    """Most stable homodimer."""
    opts: dict[str, Any] = PRIMER3_CONDITIONS | kwargs
    return _structure(primer3.calc_homodimer(_clean(seq), output_structure=True, **opts), True)


def cross_dimer(a: str, b: str, **kwargs: Any) -> Structure:
    """Most stable heterodimer between two primers."""
    opts: dict[str, Any] = PRIMER3_CONDITIONS | kwargs
    return _structure(
        primer3.calc_heterodimer(_clean(a), _clean(b), output_structure=True, **opts), True
    )


def end_stability(a: str, b: str, **kwargs: Any) -> Structure:
    """Stability of ``a``'s 3' end against ``b``.

    This is the configuration that matters: a duplex involving the 3' terminus
    can be extended into primer-dimer product, where an equally stable internal
    duplex only sequesters primer.
    """
    opts: dict[str, Any] = PRIMER3_CONDITIONS | kwargs
    return _structure(primer3.calc_end_stability(_clean(a), _clean(b), **opts))


def primer3_tm(seq: str, **kwargs: Any) -> float:
    """primer3's own Tm, as an independent cross-check on :func:`neb_tm`.

    It will not equal ``neb_tm``: primer3 applies a different salt correction
    and does not implement NEB's product-specific conventions. Use ``neb_tm``
    for anything that has to agree with the calculator the lab quotes, and this
    to notice if the two ever diverge badly.
    """
    opts: dict[str, Any] = {k: v for k, v in (PRIMER3_CONDITIONS | kwargs).items()}
    tm: float = primer3.calc_tm(_clean(seq), **opts)  # primer3 is untyped for mypy
    return tm


def longest_palindrome(seq: str) -> tuple[int, str, int]:
    """Longest reverse-complement palindrome, as (length, sequence, 1-based pos).

    A palindrome inside a primer makes it anneal to a copy of itself. Primers
    placed on polylinkers inherit one from the restriction site they sit on,
    which is the usual cause of a severe self-dimer.
    """
    s = _clean(seq)
    best = (0, "", 0)
    for i in range(len(s)):
        for j in range(i + 2, len(s) + 1, 2):
            sub = s[i:j]
            if sub == reverse_complement(sub) and len(sub) > best[0]:
                best = (len(sub), sub, i + 1)
    return best


# --------------------------------------------------------------------------- #
# Template binding
# --------------------------------------------------------------------------- #


@dataclass
class BindingSite:
    strand: str  # "+" if the primer matches the plus strand as written
    start: int  # 1-based inclusive, on the plus strand
    end: int
    mismatches: int
    three_prime_intact: bool


def find_binding_sites(
    primer: str,
    template: str,
    circular: bool = True,
    max_mismatches: int = 3,
    anchor: int = 5,
) -> list[BindingSite]:
    """Every site where ``primer`` could anneal, on both strands.

    ``three_prime_intact`` records whether the extendable 3' ``anchor`` bases
    match perfectly; a site with a mismatched 3' terminus will not prime, so a
    hit without it is not a specificity problem.
    """
    p = _clean(primer)
    ref = _clean(template)
    n = len(ref)
    scan = ref + ref[: len(p) - 1] if circular else ref
    sites: list[BindingSite] = []
    for strand in ("+", "-"):
        query = p if strand == "+" else reverse_complement(p)
        for i in range(n if circular else n - len(query) + 1):
            window = scan[i : i + len(query)]
            if len(window) < len(query):
                continue
            mm = sum(1 for x, y in zip(window, query, strict=True) if x != y)
            if mm > max_mismatches:
                continue
            if strand == "+":
                intact = query[-anchor:] == window[-anchor:]
            else:
                intact = query[:anchor] == window[:anchor]
            sites.append(BindingSite(strand, i + 1, (i + len(query) - 1) % n + 1, mm, intact))
    return sites


def amplicons(
    fwd: str, rev: str, template: str, circular: bool = True, max_length: int = 20000
) -> list[tuple[int, int, int]]:
    """``(start, end, length)`` for every product the pair makes on a template."""
    ref = _clean(template)
    n = len(ref)
    f_sites = [s for s in find_binding_sites(fwd, ref, circular) if s.mismatches == 0]
    r_sites = [s for s in find_binding_sites(rev, ref, circular) if s.mismatches == 0]
    out = []
    for f, r in itertools.product(f_sites, r_sites):
        if f.strand != "+" or r.strand != "-":
            continue
        length = (r.end - f.start) % n + 1
        if 0 < length <= max_length:
            out.append((f.start, r.end, length))
    return sorted(out, key=lambda t: t[2])


# --------------------------------------------------------------------------- #
# Vector homology and miRNA seed sites
# --------------------------------------------------------------------------- #

_REPO_ROOT = Path(__file__).resolve().parents[2]

#: TargetScanHuman 8.0 seed target sites, tier-annotated. It lives under the
#: s1e2 experiment for historical reasons -- it was generated there and has
#: three callers now. Relocating it to package data is a sensible follow-up;
#: until then this mirrors the repo-root convention used by env and plotting.
MIRNA_SEED_TSV = _REPO_ROOT / "library/l1/vector-design/mirna_seed_sites.tsv"

#: Tiers worth excluding by default. "poorly_conserved" and "unassigned" are
#: numerous and weakly supported; screening on them rejects almost everything.
DEFAULT_SEED_TIERS = frozenset({"broadly_conserved", "conserved"})


@dataclass
class SeedHit:
    site: str  # the 7mer target site found
    family: str  # the miRNA family's seed, as RNA
    tier: str
    position: int  # 1-based, within the scanned sequence


def load_mirna_seed_sites(
    tiers: frozenset[str] | set[str] = DEFAULT_SEED_TIERS,
    path: Path | None = None,
) -> dict[str, set[tuple[str, str]]]:
    """``{7mer target site: {(miRNA family, tier)}}`` for the requested tiers."""
    table = path or MIRNA_SEED_TSV
    if not table.exists():
        raise FileNotFoundError(f"{table} not found; regenerate it with s1e2 barcode_design.py")
    sites: dict[str, set[tuple[str, str]]] = {}
    for line in table.read_text().splitlines():
        if line.startswith(("#", "target_site")) or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 3 and parts[2] in tiers:
            sites.setdefault(parts[0], set()).add((parts[1], parts[2]))
    return sites


def transcript_seed_sites(
    sequence: str,
    tiers: frozenset[str] | set[str] = DEFAULT_SEED_TIERS,
    sites: dict[str, set[tuple[str, str]]] | None = None,
) -> list[SeedHit]:
    """miRNA seed target sites in a transcribed sequence.

    **Pass the plus strand of the construct, not a primer.** A seed site acts on
    the transcript, so it only means anything read in the transcribed
    orientation. Scanning a reverse primer's own sequence reads the minus strand
    and invents hits that do not exist -- which is exactly the mistake this
    function exists to stop someone repeating.

    For a primer pair, scan the amplicon's plus-strand span once rather than the
    two primers, and compare against the same span before the change so that
    pre-existing sites are not attributed to it.
    """
    table = sites if sites is not None else load_mirna_seed_sites(tiers)
    s = _clean(sequence)
    return [
        SeedHit(s[i : i + 7], family, tier, i + 1)
        for i in range(len(s) - 6)
        for family, tier in sorted(table.get(s[i : i + 7], ()))
    ]


@dataclass
class VectorHomology:
    """How much of a primer is echoed elsewhere in its own template."""

    longest_shared: int  # longest k-mer of the primer found in the template
    shared_sequence: str
    copies: int  # occurrences of that k-mer, both strands
    extra_copies: int  # copies beyond the primer's own intended site

    @property
    def concerning(self) -> bool:
        """A second copy is a mispriming and recombination risk.

        Long direct repeats in a lentiviral construct also template the
        recombination that swaps barcodes between ORFs, so this is worth
        knowing for the vector and not only for the PCR.
        """
        return self.extra_copies > 0 and self.longest_shared >= 12


def vector_homology(
    primer: str, template: str, min_k: int = 8, expected_sites: int = 1
) -> VectorHomology:
    """Longest stretch of ``primer`` that occurs anywhere in ``template``.

    ``expected_sites`` is how many times the primer is supposed to bind; copies
    beyond that are the risk. Pass 0 for a handle that is not yet installed.
    """
    p, ref = _clean(primer), _clean(template)

    def occurrences(sub: str) -> int:
        total = 0
        for pattern in {sub, reverse_complement(sub)}:
            start = 0
            while (idx := ref.find(pattern, start)) != -1:
                total += 1
                start = idx + 1
        return total

    for k in range(len(p), min_k - 1, -1):
        for i in range(len(p) - k + 1):
            sub = p[i : i + k]
            n = occurrences(sub)
            if n:
                return VectorHomology(k, sub, n, max(0, n - expected_sites))
    return VectorHomology(0, "", 0, 0)


# --------------------------------------------------------------------------- #
# Reports
# --------------------------------------------------------------------------- #


@dataclass
class PrimerReport:
    name: str
    sequence: str
    length: int
    gc_percent: float
    tm_c: float  # NEB-equivalent, the number the lab quotes
    tm_primer3_c: float  # primer3's own Tm, as a divergence check
    product: str
    hairpin: Structure
    self_dimer: Structure
    self_end_stability: Structure
    palindrome_length: int
    palindrome_seq: str
    palindrome_pos: int
    binding_sites: list[BindingSite] = field(default_factory=list)
    homology: VectorHomology | None = None

    def concerns(self) -> list[str]:
        """Flags worth a human look. Empty means nothing tripped a threshold."""
        out: list[str] = []
        if self.self_dimer.concerning:
            out.append(
                f"self-dimer {self.self_dimer.dg:.1f} kcal/mol "
                f"(structure Tm {self.self_dimer.tm:.0f} degC)"
            )
        if self.self_end_stability.found and self.self_end_stability.dg <= -5:
            out.append(
                f"3'-end self-stability {self.self_end_stability.dg:.1f} kcal/mol "
                f"(Tm {self.self_end_stability.tm:.0f} degC)"
            )
        if self.hairpin.found and self.hairpin.dg <= -3:
            out.append(f"hairpin {self.hairpin.dg:.1f} kcal/mol")
        if self.palindrome_length >= 6:
            out.append(f"{self.palindrome_length} bp palindrome {self.palindrome_seq}")
        if not 40 <= self.gc_percent <= 65:
            out.append(f"GC {self.gc_percent:.0f}%")
        if self.sequence[-1] not in "GC":
            out.append("no 3' G/C clamp")
        if abs(self.tm_c - self.tm_primer3_c) > 8:
            out.append(
                f"NEB and primer3 Tm disagree by {abs(self.tm_c - self.tm_primer3_c):.0f} degC"
            )
        perfect = [s for s in self.binding_sites if s.mismatches == 0]
        if self.binding_sites and len(perfect) != 1:
            out.append(f"{len(perfect)} perfect sites on the template")
        if self.homology is not None and self.homology.concerning:
            h = self.homology
            out.append(
                f"{h.longest_shared} nt ({h.shared_sequence}) recurs "
                f"{h.extra_copies}x elsewhere in the template"
            )
        return out


def check_primer(
    sequence: str,
    name: str = "primer",
    product: str = "Q5",
    template: str | None = None,
    circular: bool = True,
) -> PrimerReport:
    """Full check of one primer. Thermodynamics from primer3, Tm from NEB."""
    s = _clean(sequence)
    pal_len, pal_seq, pal_pos = longest_palindrome(s)
    return PrimerReport(
        name=name,
        sequence=s,
        length=len(s),
        gc_percent=gc_content(s) * 100,
        tm_c=neb_tm(s, product),
        tm_primer3_c=primer3_tm(s),
        product=product,
        hairpin=hairpin(s),
        self_dimer=self_dimer(s),
        self_end_stability=end_stability(s, s),
        palindrome_length=pal_len,
        palindrome_seq=pal_seq,
        palindrome_pos=pal_pos,
        binding_sites=find_binding_sites(s, template, circular) if template else [],
        homology=vector_homology(s, template) if template else None,
    )


@dataclass
class PairReport:
    forward: PrimerReport
    reverse: PrimerReport
    tm_difference: float
    anneal_c: float
    anneal_at_ceiling: bool
    cross_dimer: Structure
    cross_end_stability: Structure
    products: list[tuple[int, int, int]] = field(default_factory=list)
    amplicon_seed_sites: list[SeedHit] = field(default_factory=list)

    def concerns(self) -> list[str]:
        out = [f"forward: {c}" for c in self.forward.concerns()]
        out += [f"reverse: {c}" for c in self.reverse.concerns()]
        if self.tm_difference > 5:
            out.append(f"Tm difference {self.tm_difference:.0f} degC exceeds NEB's 5 degC limit")
        if self.anneal_at_ceiling:
            out.append(
                "annealing temperature is NEB's 72 degC ceiling, not a match - "
                "no headroom to raise stringency"
            )
        if self.cross_dimer.concerning:
            out.append(
                f"cross-dimer {self.cross_dimer.dg:.1f} kcal/mol "
                f"(structure Tm {self.cross_dimer.tm:.0f} degC)"
            )
        if self.cross_end_stability.found and self.cross_end_stability.dg <= -5:
            out.append(f"3'-end cross-stability {self.cross_end_stability.dg:.1f} kcal/mol")
        if self.products and len(self.products) != 1:
            out.append(f"{len(self.products)} products on the template")
        if self.amplicon_seed_sites:
            families = sorted({h.family for h in self.amplicon_seed_sites})
            out.append(
                f"{len(self.amplicon_seed_sites)} conserved miRNA seed site(s) in the "
                f"amplicon ({', '.join(families)}) - compare against the same span "
                f"before any edit, most are pre-existing"
            )
        return out

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _amplicon_seeds(
    forward: str, reverse: str, template: str | None, circular: bool
) -> list[SeedHit]:
    """Seed sites in the amplicon, read on the plus strand.

    Scanning here rather than on the two primers is deliberate: a seed site acts
    on the transcript, so it is only meaningful in the transcribed orientation,
    and a reverse primer's own sequence is the wrong strand.
    """
    if not template:
        return []
    products = amplicons(forward, reverse, template, circular)
    if len(products) != 1:
        return []
    start, end, _ = products[0]
    ref = _clean(template)
    span = ref[start - 1 : end] if start <= end else ref[start - 1 :] + ref[:end]
    try:
        return transcript_seed_sites(span)
    except FileNotFoundError:
        return []


def check_pair(
    forward: str,
    reverse: str,
    template: str | None = None,
    product: str = "Q5",
    circular: bool = True,
    fwd_name: str = "fwd",
    rev_name: str = "rev",
) -> PairReport:
    """Full check of a primer pair, including products on a template."""
    f = check_primer(forward, fwd_name, product, template, circular)
    r = check_primer(reverse, rev_name, product, template, circular)
    ta, capped = anneal_temp(forward, reverse, product)
    return PairReport(
        forward=f,
        reverse=r,
        tm_difference=abs(f.tm_c - r.tm_c),
        anneal_c=ta,
        anneal_at_ceiling=capped,
        cross_dimer=cross_dimer(forward, reverse),
        cross_end_stability=end_stability(forward, reverse),
        products=amplicons(forward, reverse, template, circular) if template else [],
        amplicon_seed_sites=_amplicon_seeds(forward, reverse, template, circular),
    )


# --------------------------------------------------------------------------- #
# BLAST
# --------------------------------------------------------------------------- #

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


#: PDB-style accessions: ``9SYR_S2``, ``4V6X_BA``. These are structure entries,
#: and NCBI classifies a viral RNA inside a human ribosome complex under
#: taxid 9606 -- so an organism-restricted search still returns them. They are
#: not genomic sequence and cannot be a gDNA off-target.
_STRUCTURE_ACCESSION = re.compile(r"^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]{1,3}$")


def is_structure_entry(subject: str) -> bool:
    """Whether a BLAST subject is a PDB structure rather than a sequence record."""
    if not subject:
        return False
    accession = subject.split()[0]
    if _STRUCTURE_ACCESSION.match(accession):
        return True
    rest = subject[len(accession) :].strip()
    return rest.startswith("Chain ")


@dataclass
class BlastSummary:
    name: str
    length: int
    alignments: int
    full_length_hits: int  # sequence records only; structures excluded
    structure_hits: int  # PDB entries, reported separately
    best_three_prime_run: int  # sequence records only
    best_subject: str

    def primes_off_target(self, min_run: int = 17) -> bool:
        """Whether any hit anchors enough 3' sequence to plausibly prime."""
        return self.best_three_prime_run >= min_run


def blast_primers(
    primers: dict[str, str],
    entrez_query: str = "txid9606[ORGN]",
    database: str = "core_nt",
    hitlist_size: int = 100,
    timeout_s: int = 1800,
    poll_s: int = 20,
) -> tuple[str, str]:
    """Submit primers to NCBI BLAST. Returns ``(rid, text_report)``.

    Short primers need ``task=blastn`` with a small word size and a permissive
    E-value; megablast simply misses them. Defaults restrict to human, which is
    the relevant background when amplifying from a human cell line's gDNA.
    """
    query = "\n".join(f">{n}\n{_clean(s)}" for n, s in primers.items())
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "MEGABLAST": "off",
        "DATABASE": database,
        "QUERY": query,
        "WORD_SIZE": "7",
        "EXPECT": "1000",
        "HITLIST_SIZE": str(hitlist_size),
        "FILTER": "F",
    }
    if entrez_query:
        params["ENTREZ_QUERY"] = entrez_query
    req = urllib.request.Request(BLAST_URL, data=urllib.parse.urlencode(params).encode())
    with urllib.request.urlopen(req, timeout=120) as fh:
        html = fh.read().decode("utf-8", "replace")
    match = re.search(r"RID = (\S+)", html)
    if not match:
        raise RuntimeError("NCBI did not return a RID; the submission was rejected")
    rid = match.group(1)

    deadline = time.time() + timeout_s
    while time.time() < deadline:
        time.sleep(poll_s)
        status_url = f"{BLAST_URL}?{urllib.parse.urlencode({'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid})}"
        with urllib.request.urlopen(status_url, timeout=120) as fh:
            status = fh.read().decode("utf-8", "replace")
        if "Status=WAITING" in status:
            continue
        if "Status=READY" in status:
            break
        raise RuntimeError(f"BLAST {rid} failed: {status[:200]}")
    else:
        raise TimeoutError(f"BLAST {rid} did not finish within {timeout_s}s")

    fetch = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "Text",
        "ALIGNMENTS": str(hitlist_size),
        "DESCRIPTIONS": str(hitlist_size),
    }
    with urllib.request.urlopen(f"{BLAST_URL}?{urllib.parse.urlencode(fetch)}", timeout=600) as fh:
        return rid, fh.read().decode("utf-8", "replace")


def parse_blast_text(report: str) -> dict[str, list[dict[str, Any]]]:
    """Parse a BLAST text report into per-query alignment records.

    Kept separate from ``blast_primers`` so it can be tested on a fixture
    without a network round trip.
    """
    out: dict[str, list[dict[str, Any]]] = {}
    current: str | None = None
    subject = ""
    ident = ""
    lines = report.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("Query= "):
            current = line[len("Query= ") :].strip()
            out[current] = []
            subject = ""
        elif line.startswith(">"):
            subject = line[1:].strip()
        elif line.strip().startswith("Identities ="):
            ident = line.strip()
        elif line.startswith("Query ") and current is not None:
            parts = line.split()
            if len(parts) == 4 and parts[1].isdigit():
                qstart, qseq, qend = int(parts[1]), parts[2], int(parts[3])
                col = line.index(qseq)
                bars = lines[i + 1][col : col + len(qseq)] if i + 1 < len(lines) else ""
                # Subject coordinates are needed to pair hits into an amplicon.
                # On a Plus/Minus alignment send < sstart, which is how the
                # orientation is recovered.
                sline = lines[i + 2] if i + 2 < len(lines) else ""
                sparts = sline.split()
                has_coords = len(sparts) == 4 and sparts[1].isdigit() and sparts[3].isdigit()
                out[current].append(
                    {
                        "subject": subject,
                        "identities": ident,
                        "qstart": qstart,
                        "qend": qend,
                        "qseq": qseq,
                        "bars": bars,
                        "sstart": int(sparts[1]) if has_coords else None,
                        "send": int(sparts[3]) if has_coords else None,
                    }
                )
                i += 2
        i += 1
    return out


def three_prime_run(record: dict[str, Any], query_length: int) -> int:
    """Perfectly matched bases counting back from the primer's 3' end.

    Zero when the alignment stops short of the 3' terminus -- such a hit cannot
    be extended and so cannot prime, however high its identity.
    """
    if record["qend"] < query_length:
        return 0
    run = 0
    for q, b in zip(reversed(record["qseq"]), reversed(record["bars"]), strict=False):
        if q == "-":
            continue
        if b == "|":
            run += 1
        else:
            break
    return run


def summarise_blast(report: str, primers: dict[str, str]) -> list[BlastSummary]:
    """Reduce a BLAST report to what matters for priming, per primer."""
    parsed = parse_blast_text(report)
    out = []
    for name, seq in primers.items():
        qlen = len(_clean(seq))
        records = parsed.get(name, [])
        runs = [(three_prime_run(r, qlen), r) for r in records]
        runs = [(n, r) for n, r in runs if n > 0]
        structural = [(n, r) for n, r in runs if is_structure_entry(r["subject"])]
        sequence = [(n, r) for n, r in runs if not is_structure_entry(r["subject"])]
        sequence.sort(key=lambda t: -t[0])
        best_run, best_subject = (
            (sequence[0][0], sequence[0][1]["subject"]) if sequence else (0, "")
        )
        out.append(
            BlastSummary(
                name=name,
                length=qlen,
                alignments=len(records),
                full_length_hits=sum(1 for n, _ in sequence if n >= qlen),
                structure_hits=sum(1 for n, _ in structural if n >= qlen),
                best_three_prime_run=best_run,
                best_subject=best_subject[:100],
            )
        )
    return out


@dataclass
class GenomicAmplicon:
    """An off-target product a primer pair could make on a BLAST subject."""

    subject: str
    forward_primer: str
    reverse_primer: str
    start: int
    end: int
    length: int
    forward_3p_match: int
    reverse_3p_match: int

    @property
    def plausible(self) -> bool:
        """Both ends anchored enough 3' sequence to prime at high stringency."""
        return self.forward_3p_match >= 15 and self.reverse_3p_match >= 15


def genomic_amplicon_risk(
    report: str,
    primers: dict[str, str],
    max_length: int = 4000,
    min_three_prime: int = 12,
) -> list[GenomicAmplicon]:
    """Find off-target products, the way Primer-BLAST does.

    A single off-target primer site makes no product. A product needs a forward
    and a reverse site on the same subject, pointing at each other, close enough
    to amplify -- which is why a per-primer specificity table can look clean
    while the pair is not. This pairs up BLAST hits to catch that case.

    Only hits whose 3' end is anchored by at least ``min_three_prime`` matched
    bases are considered, since the rest cannot be extended.
    """
    parsed = parse_blast_text(report)
    # subject -> list of (primer_name, orientation, position, 3' run)
    by_subject: dict[str, list[tuple[str, str, int, int, int]]] = {}
    for name, seq in primers.items():
        qlen = len(_clean(seq))
        for record in parsed.get(name, []):
            run = three_prime_run(record, qlen)
            if run < min_three_prime:
                continue
            if is_structure_entry(record["subject"]):
                continue  # a PDB chain is not a genomic locus
            sstart, send = record.get("sstart"), record.get("send")
            if sstart is None or send is None:
                continue
            orientation = "+" if send >= sstart else "-"
            by_subject.setdefault(record["subject"], []).append(
                (name, orientation, sstart, send, run)
            )

    out: list[GenomicAmplicon] = []
    for subject, hits in by_subject.items():
        for (n1, o1, s1, e1, r1), (n2, o2, s2, e2, r2) in itertools.combinations(hits, 2):
            if n1 == n2 or o1 == o2:
                continue
            # The plus-strand hit must sit 5' of the minus-strand hit.
            plus = (n1, s1, e1, r1) if o1 == "+" else (n2, s2, e2, r2)
            minus = (n2, s2, e2, r2) if o1 == "+" else (n1, s1, e1, r1)
            start = min(plus[1], plus[2])
            end = max(minus[1], minus[2])
            length = end - start + 1
            if not 0 < length <= max_length:
                continue
            out.append(
                GenomicAmplicon(
                    subject=subject[:100],
                    forward_primer=plus[0],
                    reverse_primer=minus[0],
                    start=start,
                    end=end,
                    length=length,
                    forward_3p_match=plus[3],
                    reverse_3p_match=minus[3],
                )
            )
    return sorted(out, key=lambda a: -min(a.forward_3p_match, a.reverse_3p_match))


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _format_pair(rep: PairReport) -> str:
    lines = []
    lines.append(
        f"{'':<10} {'sequence':<28} {'nt':>3} {'GC%':>6} {'Tm':>6} {'p3 Tm':>6} "
        f"{'hairpin':>8} {'selfdim':>9} {'dimer Tm':>9} {'3p end':>8}"
    )
    lines.append("-" * 105)
    for label, p in (("forward", rep.forward), ("reverse", rep.reverse)):
        dimer_tm = f"{p.self_dimer.tm:.0f}" if p.self_dimer.found else "-"
        lines.append(
            f"{label:<10} {p.sequence:<28} {p.length:>3} {p.gc_percent:>5.1f}% "
            f"{p.tm_c:>6.1f} {p.tm_primer3_c:>6.1f} {p.hairpin.dg:>8.2f} "
            f"{p.self_dimer.dg:>9.2f} {dimer_tm:>9} {p.self_end_stability.dg:>8.2f}"
        )
    lines.append("")
    ceiling = "  (NEB ceiling, not a match)" if rep.anneal_at_ceiling else ""
    lines.append(
        f"Tm difference {rep.tm_difference:.1f} degC    anneal {rep.anneal_c:.0f} degC{ceiling}"
    )
    lines.append(
        f"cross-dimer {rep.cross_dimer.dg:.2f} kcal/mol "
        f"(Tm {rep.cross_dimer.tm:.0f} degC), 3'-end {rep.cross_end_stability.dg:.2f}"
    )
    if rep.products:
        for start, end, length in rep.products:
            lines.append(f"product {start}..{end} = {length} bp")
    lines.append("")
    concerns = rep.concerns()
    if concerns:
        lines.append("Concerns:")
        lines.extend(f"  - {c}" for c in concerns)
    else:
        lines.append("No threshold tripped.")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m deliverome_analysis.primers",
        description="Check a primer or primer pair: Tm, dimers, template sites, BLAST.",
    )
    parser.add_argument("--fwd", required=True, help="forward primer sequence")
    parser.add_argument("--rev", help="reverse primer sequence")
    parser.add_argument("--template", type=Path, help="FASTA or GenBank file to map sites on")
    parser.add_argument("--linear", action="store_true", help="template is linear, not circular")
    parser.add_argument("--product", default="Q5", choices=sorted(NEB_PRODUCTS))
    parser.add_argument("--blast", action="store_true", help="also run NCBI BLAST (slow)")
    parser.add_argument("--entrez", default="txid9606[ORGN]", help="BLAST organism filter")
    parser.add_argument("--json", type=Path, help="write the full report as JSON")
    args = parser.parse_args(argv)

    template = None
    if args.template:
        from Bio import SeqIO

        fmt = "genbank" if args.template.suffix.lower() in {".gb", ".gbk"} else "fasta"
        template = str(SeqIO.read(args.template, fmt).seq).upper()

    if args.rev:
        rep = check_pair(args.fwd, args.rev, template, args.product, circular=not args.linear)
        print(_format_pair(rep))
        payload = rep.to_dict()
    else:
        p = check_primer(args.fwd, "primer", args.product, template, not args.linear)
        print(f"{p.sequence}  {p.length} nt  GC {p.gc_percent:.1f}%  Tm {p.tm_c:.1f} degC")
        for c in p.concerns() or ["No threshold tripped."]:
            print(f"  - {c}")
        payload = asdict(p)

    if args.blast:
        primers = {"fwd": args.fwd} | ({"rev": args.rev} if args.rev else {})
        rid, report = blast_primers(primers, entrez_query=args.entrez)
        print(f"\nBLAST {rid}")
        for s in summarise_blast(report, primers):
            verdict = "could prime" if s.primes_off_target() else "cannot prime"
            print(
                f"  {s.name}: {s.full_length_hits} full-length hits, "
                f"best 3' run {s.best_three_prime_run} nt -> {verdict}"
            )
            if s.best_subject:
                print(f"      {s.best_subject}")
        payload["blast"] = {
            "rid": rid,
            "summaries": [asdict(s) for s in summarise_blast(report, primers)],
        }

    if args.json:
        args.json.write_text(json.dumps(payload, indent=2, default=str))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
