"""
condition_filter.py
-------------------
Filters parsed variants by user-defined condition keywords at runtime.
Supports partial / case-insensitive matching so users don't need exact names.
"""

from src.vcf_parser import Variant

# ── Built-in condition aliases (maps friendly names → search tokens) ──────────
CONDITION_ALIASES: dict[str, list[str]] = {
    "diabetes":            ["diabetes", "t2d", "maturity onset"],
    "cad":                 ["coronary artery disease", "cad", "familial hypercholesterolemia"],
    "breast cancer":       ["breast cancer", "brca"],
    "ovarian cancer":      ["ovarian cancer"],
    "alzheimer":           ["alzheimer"],
    "sickle cell":         ["sickle cell"],
    "cystic fibrosis":     ["cystic fibrosis", "cftr"],
    "stroke":              ["stroke"],
    "hypertension":        ["hypertension"],
    "obesity":             ["obesity"],
}


def resolve_aliases(condition: str) -> list[str]:
    """
    Given a user-supplied condition string, return a list of search tokens.
    Falls back to the raw string if no alias match is found.
    """
    lower = condition.lower().strip()
    for key, tokens in CONDITION_ALIASES.items():
        if lower == key or lower in tokens:
            return tokens
    return [lower]


def filter_by_condition(
    variants: list[Variant],
    conditions: list[str],
    include_benign: bool = False,
    only_carriers: bool = True,
) -> list[Variant]:
    """
    Filter variants to those matching any of the requested conditions.

    Args:
        variants:        Full list of Variant objects from the parser.
        conditions:      List of condition strings supplied by the user.
        include_benign:  If False, drops Benign / VUS variants from results.
        only_carriers:   If True, drops homozygous-ref (non-carrier) records.

    Returns:
        Filtered list of Variant objects.
    """
    # Build a flat set of search tokens from all requested conditions
    search_tokens: list[str] = []
    for cond in conditions:
        search_tokens.extend(resolve_aliases(cond))

    results = []
    for v in variants:
        # Carrier filter
        if only_carriers and not v.is_variant:
            continue

        # Benign filter
        if not include_benign and v.risk_tier == "Benign / VUS":
            continue

        # Condition match — check against all disease strings for this variant
        disease_text = v.disease_name.lower()
        cond_text    = " ".join(v.conditions).lower()
        combined     = disease_text + " " + cond_text

        if any(tok in combined for tok in search_tokens):
            results.append(v)

    return results


def group_by_condition(variants: list[Variant], conditions: list[str]) -> dict[str, list[Variant]]:
    """
    Group filtered variants by the matched condition for report sectioning.
    A single variant may appear under multiple conditions if it is pleiotropic.
    """
    groups: dict[str, list[Variant]] = {}
    for cond in conditions:
        tokens = resolve_aliases(cond)
        label  = cond.title()
        groups[label] = []
        for v in variants:
            combined = (v.disease_name + " " + " ".join(v.conditions)).lower()
            if any(tok in combined for tok in tokens):
                groups[label].append(v)
    return groups


def list_available_conditions() -> list[str]:
    """Return all built-in condition alias keys."""
    return sorted(CONDITION_ALIASES.keys())
