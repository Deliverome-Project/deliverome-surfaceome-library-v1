"""Reviewed S1E3 partner requirements; see CORECEPTOR_AUDIT.md for evidence."""

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import isfinite

PARTNER_MIN_TPM = 5


@dataclass(frozen=True)
class PartnerRule:
    symbol: str
    observed_partners: frozenset[str]
    required_groups: tuple[tuple[str, ...], ...]
    integral_surface_complex: bool = False


REVIEWED_RULES = {
    "P54709": PartnerRule(
        "ATP1B3",
        frozenset({"ATP1A1", "ATP1A2", "ATP1A3", "ATP1A4"}),
        (("ATP1A1",),),
        True,
    ),
    "Q9NV96": PartnerRule(
        "CDC50A",
        frozenset(
            {
                "ATP8A1",
                "ATP8A2",
                "ATP8B1",
                "ATP8B2",
                "ATP10A",
                "ATP10B",
                "ATP10D",
                "ATP11A",
                "ATP11B",
                "ATP11C",
            }
        ),
        (("ATP8B1", "ATP8B2"),),
        True,
    ),
    "O15547": PartnerRule(
        "P2RX6",
        frozenset({"P2RX2", "P2RX4"}),
        (("P2RX4",),),
        True,
    ),
    "Q30201": PartnerRule("HFE", frozenset({"B2M"}), (("B2M",),), True),
    "O43520": PartnerRule(
        "ATP8B1",
        frozenset({"CDC50A", "CDC50B"}),
        (("CDC50A", "CDC50B"),),
    ),
    "Q3MIR4": PartnerRule(
        "CDC50B",
        frozenset({"ATP8B1", "ATP8B4"}),
        (("ATP8B1",),),
    ),
}


def groups_satisfied(
    groups: Sequence[Sequence[str]],
    partner_median: Callable[[str], float | None],
) -> bool:
    """Require every partner role, with at least one measured alternative per role."""
    return bool(groups) and all(
        any(
            (median := partner_median(partner)) is not None
            and isfinite(median)
            and median >= PARTNER_MIN_TPM
            for partner in group
        )
        for group in groups
    )


def reviewed_partner_ok(accession, partners, partner_median):
    """Return None for unreviewed targets and False for changed/missing partner lists."""
    rule = REVIEWED_RULES.get(accession)
    if rule is None:
        return None
    if frozenset(partners) != rule.observed_partners:
        return False
    return groups_satisfied(rule.required_groups, partner_median)


def integral_complex_ok(record, partners, partner_median):
    """A reviewed, anchored surface complex may satisfy the localization reason gate."""
    rule = REVIEWED_RULES.get(record.get("uniprot"))
    topology = record.get("canonical_topology") or {}
    return bool(
        rule
        and rule.integral_surface_complex
        and (topology.get("tm_helix_count") or 0) >= 1
        and reviewed_partner_ok(record.get("uniprot"), partners, partner_median)
    )
