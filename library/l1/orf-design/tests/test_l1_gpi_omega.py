"""No GPI-anchored member carries its cassette after the omega site.

The omega site deletes everything after it with the propeptide, so a cassette placed there
leaves the cell with the propeptide. ``gpi_set.tsv`` records ``propeptide_aa``, the residues
after omega in the anchored form NetGPI scored best. The GPI signal is C-terminal, so the
distance from omega to the C-terminus carries over to a member that shares that end, even
where ``omega_site`` itself is counted on a fragment (CPM 15, GPC5 27). A cassette is kept
only if more native residues follow it than the propeptide holds.

Reads committed tables only. Members join on design_id.
"""

from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]

#: Arm tables, and the column giving the native residue the cassette follows.
ARM_TABLES = {
    "constructs_l1_cterm.tsv": "cassette_after_residue",
    "constructs_l1_internal.tsv": "cassette_after_residue",
    "constructs_l1_nterm.tsv": "native_after_residue",
    "constructs_l1_multipass.tsv": None,
    "constructs_l1_rescue.tsv": None,
    "constructs_l1_type3.tsv": None,
}


def _tsv(name: str) -> list[dict]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _gpi() -> dict[str, dict]:
    return {r["gene_symbol"]: r for r in _tsv("gpi_set.tsv") if r["gpi_anchored"] == "1"}


def test_gpi_set_records_full_length_omega():
    for gene, r in _gpi().items():
        if not r["omega_site"]:
            continue
        assert r["omega_form"] and r["omega_form_length"], f"{gene}: omega has no form"
        assert int(r["propeptide_aa"]) == int(r["omega_form_length"]) - int(r["omega_site"])
        if r["anchored_in_reference_form"] == "1":
            assert r["omega_site_reference"], f"{gene}: reference anchored, no reference omega"


def test_no_gpi_member_is_tagged_after_its_omega_site():
    gpi = _gpi()
    checked, bad = 0, []
    for table, col in ARM_TABLES.items():
        for r in _tsv(table):
            g = gpi.get(r["gene_symbol"])
            if not g or not g["propeptide_aa"]:
                continue
            assert col, f"{table}: GPI member {r['design_id']} in a table with no site column"
            after = int(r["native_aa"]) - int(r[col])
            checked += 1
            if after <= int(g["propeptide_aa"]):
                bad.append((table, r["design_id"], r["gene_symbol"], after, g["propeptide_aa"]))
    assert checked, "no GPI arms found"
    assert not bad, f"cassette after omega: {bad}"
