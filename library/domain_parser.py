"""Utilities for parsing domain descriptions.

This module provides functions to normalize textual descriptions of
protein domains and to map them to canonical domain type and
localization tags.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Union
import logging
import re

from .domain_dictionaries import DOMAIN_TYPE_MAP, DOMAIN_LOC_MAP, normalize_text

LOGGER = logging.getLogger(__name__)


@dataclass
class DomainParseResult:
    """Parsed information for a domain description."""

    domain_type: Optional[str] = None
    domain_loc: List[str] = field(default_factory=list)
    domain_index: Optional[Union[int, List[int]]] = None
    domain_source_variant: Optional[str] = None
    domain_reason: List[str] = field(default_factory=list)


_TYPE_PRIORITY = {
    "BROMO_BD1": 1,
    "BROMO_BD2": 2,
    "BROMO_TANDEM": 3,
    "BROMO": 4,
    "KD": 5,
    "LBD": 6,
    "DBD": 7,
    "MOTOR": 8,
    "GYRB": 9,
    "PKD": 10,
    "KI": 11,
    "TMD": 12,
    "ECD": 13,
    "ICD": 14,
    "N_TERM": 15,
    "C_TERM": 16,
}

_LOC_PRIORITY = {"ICD": 1, "ECD": 2, "C_TERM": 3, "N_TERM": 4}

_DROPPED_WORDS = {"receptor", "mutant", "protein", "subtype"}


def _match_from_dict(
    text: str, lookup: dict[str, str], priority: dict[str, int]
) -> Optional[str]:
    """Return the best canonical tag found in *lookup* for *text*."""
    matches: List[str] = []
    for alias, canonical in lookup.items():
        tokens = [re.escape(tok) for tok in alias.split()]
        pattern = r"\b" + r"\s+".join(tokens) + r"\b"
        if re.search(pattern, text):
            matches.append(canonical)
    if not matches:
        return None
    return min(matches, key=lambda c: priority.get(c, 100))


def parse_domain(text: str) -> DomainParseResult:
    """Parse domain information from *text*.

    Parameters
    ----------
    text: str
        Raw string possibly describing a protein domain.

    Returns
    -------
    DomainParseResult
        Structured information about the detected domain.
    """
    result = DomainParseResult()
    if not text:
        return result

    normalized, norm_reasons = normalize_text(text)
    variants = [v.strip() for v in normalized.split("|") if v.strip()]
    for variant in variants:
        cleaned = re.sub(
            r"\b(" + "|".join(_DROPPED_WORDS) + r")\b", "", variant
        ).strip()
        if not cleaned:
            continue

        domain_type = _match_from_dict(cleaned, DOMAIN_TYPE_MAP, _TYPE_PRIORITY)
        domain_loc = _match_from_dict(cleaned, DOMAIN_LOC_MAP, _LOC_PRIORITY)
        index: Optional[Union[int, List[int]]] = None


        has1 = bool(
            re.search(r"\bbd1\b|bromodomain 1|bromo domain 1|domain 1|\b1\b", cleaned)
        )
        has2 = bool(
            re.search(r"\bbd2\b|bromodomain 2|bromo domain 2|domain 2|\b2\b", cleaned)
        )
        if domain_type in {"BROMO", "BROMO_BD1", "BROMO_BD2"}:
            if has1 and has2:
                domain_type = "BROMO_TANDEM"
                index = [1, 2]
                result.domain_reason.append("DERIVED_FROM_BROMODOMAIN12")
            elif domain_type == "BROMO" and has1:
                domain_type = "BROMO_BD1"
                index = 1
                result.domain_reason.append("DERIVED_FROM_BROMODOMAIN1")
            elif domain_type == "BROMO" and has2:
                domain_type = "BROMO_BD2"
                index = 2
                result.domain_reason.append("DERIVED_FROM_BROMODOMAIN2")
            elif domain_type == "BROMO_BD1" and index is None:
                index = 1
            elif domain_type == "BROMO_BD2" and index is None:
                index = 2

        if (
            domain_type == "BROMO_TANDEM"
            and index is None
            and re.search(r"\bbrd4\b", cleaned)
        ):

            index = [1, 2]

        if index is None:
            m = re.search(r"\bdomain\s*([0-9]+)\b", cleaned)
            if m:
                index = int(m.group(1))

        if domain_type or domain_loc or index is not None:
            result.domain_source_variant = variant
            if norm_reasons:
                result.domain_reason.extend(norm_reasons)
            if domain_type:
                if not result.domain_type or _TYPE_PRIORITY[
                    domain_type
                ] < _TYPE_PRIORITY.get(result.domain_type, 100):
                    result.domain_type = domain_type
                result.domain_reason.append(f"MATCH_TYPE_{domain_type}")
            if domain_loc:
                if domain_loc not in result.domain_loc:
                    result.domain_loc.append(domain_loc)
                result.domain_reason.append(f"MATCH_LOC_{domain_loc}")
            if index is not None:
                result.domain_index = index

    if (
        not result.domain_type
        and not result.domain_loc
        and re.search(r"\bdomain\b", normalized)
    ):
        result.domain_reason.append("FALLBACK_DOMAIN_WORD")

    return result


__all__ = ["parse_domain", "DomainParseResult"]
