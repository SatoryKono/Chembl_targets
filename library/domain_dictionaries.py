"""Dictionaries for domain type and localization normalization.

The dictionaries map various textual representations of domains and
localizations to canonical tags. Keys are normalized using
:func:`normalize_text` to ensure consistent matching.
"""

from __future__ import annotations

from typing import Dict, List, Tuple
import re


def normalize_text(text: str) -> Tuple[str, List[str]]:
    """Normalize a domain description fragment.

    Parameters
    ----------
    text: str
        Raw text describing a domain.

    Returns
    -------
    tuple[str, list[str]]
        A tuple of the normalized text and a list of normalization
        reason codes applied during processing.
    """
    reasons: List[str] = []
    s = re.sub(r"[-_]+", " ", text.lower())
    s = re.sub(r"(?:^| )(c|n)terminal\b", r" \1 terminal", s)
    s = s.replace("cterminal", "c terminal").replace("ndomain", "n domain")
    s = s.replace("cmet", "c met").replace("cdomain", "c domain")

    def _replace_bromo(match: re.Match[str]) -> str:
        digit = match.group(2)
        if digit == "1":
            reasons.append("NORM_SLUG_BROMO1")
        elif digit == "2":
            reasons.append("NORM_SLUG_BROMO2")
        return f"{match.group(1)} {digit}"

    s = re.sub(r"(bromodomain)(\d)", _replace_bromo, s)
    s = re.sub(r"(domain)(\d)", r"\1 \2", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s, reasons


_RAW_DOMAIN_TYPES: Dict[str, str] = {
    "KD": "kinase domain|tyrosine kinase domain|ptk domain|catalytic domain",
    "LBD": "ligand binding domain|lbd|hormone binding domain",
    "DBD": "dna binding domain|dbd|zf domain|zinc finger",
    "BROMO": "bromodomain|bromo domain|bd",
    "BROMO_BD1": "bromodomain1|bromo domain 1|bromo domain-1|bd1|domain1|domain-1",
    "BROMO_BD2": "bromodomain2|bromo domain 2|bromo domain-2|bd2|domain2|domain-2",
    "MOTOR": "motor domain",
    "ECD": "ectodomain|ecd|extracellular domain",
    "BROMO_TANDEM": "tandem domain|tandem bromodomain",
    "GYRB": "gyrb domain",
    "KI": "kinase insert|ki domain",
    "PKD": "pseudokinase domain",
    "TMD": "transmembrane domain|tm domain|tmd",
}

_RAW_DOMAIN_LOCS: Dict[str, str] = {
    "ICD": "cytoplasmic domain|intracellular domain|icd",
    "ECD": "extracellular domain|ecto domain|ectodomain|ecd",
    "C_TERM": "c-terminal domain|c terminal domain|c-domain|cdomain|cterminal",
    "N_TERM": "n-terminal domain|n terminal domain|n-domain|ndomain|nterminal",
}


def _build_dict(raw_map: Dict[str, str]) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for canonical, aliases in raw_map.items():
        for alias in aliases.split("|"):
            normalized, _ = normalize_text(alias)
            result[normalized] = canonical
    return result


DOMAIN_TYPE_MAP: Dict[str, str] = _build_dict(_RAW_DOMAIN_TYPES)
DOMAIN_LOC_MAP: Dict[str, str] = _build_dict(_RAW_DOMAIN_LOCS)

__all__ = ["DOMAIN_TYPE_MAP", "DOMAIN_LOC_MAP", "normalize_text"]
