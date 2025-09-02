"""Text normalization transformations for protein target names."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
import unicodedata
from typing import Dict, List, Sequence, Tuple

logger = logging.getLogger(__name__)

# Default dictionaries -----------------------------------------------------

GREEK_LETTERS: Dict[str, str] = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "ε": "epsilon",
    "ζ": "zeta",
    "η": "eta",
    "θ": "theta",
    "ι": "iota",
    "κ": "kappa",
    "λ": "lambda",
    "μ": "mu",
    "ν": "nu",
    "ξ": "xi",
    "ο": "omicron",
    "π": "pi",
    "ρ": "rho",
    "σ": "sigma",
    "τ": "tau",
    "υ": "upsilon",
    "φ": "phi",
    "χ": "chi",
    "ψ": "psi",
    "ω": "omega",
}

SUPERSCRIPTS: Dict[str, str] = {
    "¹": "1",
    "²": "2",
    "³": "3",
    "⁴": "4",
    "⁵": "5",
    "⁶": "6",
    "⁷": "7",
    "⁸": "8",
    "⁹": "9",
    "⁰": "0",
    "₁": "1",
    "₂": "2",
    "₃": "3",
    "₄": "4",
    "₅": "5",
    "₆": "6",
    "₇": "7",
    "₈": "8",
    "₉": "9",
    "₀": "0",
}

STOP_WORDS: Sequence[str] = (
    "protein",
    "receptor",
    "isoform",
    "fragment",
    "subunit",
    "chain",
    "precursor",
    "like",
    "putative",
    "probable",
    "predicted",
    "family",
)

RECEPTOR_RULES: Sequence[Tuple[re.Pattern[str], str, str]] = (
    # pattern, replacement, gene-like candidate
    (re.compile(r"beta2\s+adrenergic\s+receptor"), "beta2 adrenergic", "adrb2"),
    (re.compile(r"dopamine\s+d2\s+receptor"), "dopamine d2", "drd2"),
    (re.compile(r"serotonin\s+5-ht1a\s+receptor"), "5-ht1a serotonin", "htr1a"),
    (re.compile(r"histamine\s+h3\s+receptor"), "histamine h3", "hrh3"),
)

CANDIDATE_PREFIXES: Dict[str, str] = {
    "h3": "hrh3",
    "d2": "drd2",
    "beta2": "adrb2",
    "5-ht1a": "htr1a",
}

ROMAN_NUMERALS: Dict[str, str] = {"ii": "2", "iii": "3", "iv": "4"}

CONTROL_CHARS_RE = re.compile(r"[\u0000-\u001F\u007F]")
MULTI_SPACE_RE = re.compile(r"[\s\t]+")
TYPO_QUOTES_RE = re.compile(r"[“”«»„]|’")
LONG_DASH_RE = re.compile(r"[–—]")
PAREN_RE = re.compile(r"\(([^(]*)\)|\[([^\]]*)\]|\{([^}]*)\)")
TOKEN_SPLIT_RE = re.compile(r"[\s\-_/,:;\.]+")
LETTER_NUM_SPACE_RE = re.compile(r"(?<=\b)([a-z])\s+([0-9])(?=\b)", re.I)
NUM_LETTER_SPACE_RE = re.compile(r"(?<=\b)([0-9])\s+([a-z])(?=\b)", re.I)
HYPHEN_SPACE_RE = re.compile(r"\s*-\s*")


def sanitize_text(text: str) -> str:
    """Remove control characters and normalize whitespace."""
    text = CONTROL_CHARS_RE.sub("", text)
    text = text.replace("\ufeff", "").replace("\xa0", " ")
    text = MULTI_SPACE_RE.sub(" ", text)
    return text.strip()


def normalize_unicode(text: str) -> str:
    """Apply Unicode NFKC normalization and lowercase."""
    text = unicodedata.normalize("NFKC", text)
    text = text.lower()
    text = TYPO_QUOTES_RE.sub("'", text)
    text = LONG_DASH_RE.sub("-", text)
    return text


def replace_specials(
    text: str,
    greek: Dict[str, str] | None = None,
    superscripts: Dict[str, str] | None = None,
) -> str:
    """Replace greek letters and superscripts using mappings."""
    greek = greek or GREEK_LETTERS
    superscripts = superscripts or SUPERSCRIPTS
    for k, v in greek.items():
        text = text.replace(k, v)
    for k, v in superscripts.items():
        text = text.replace(k, v)
    return text


def replace_roman_numerals(text: str) -> str:
    """Replace standalone Roman numerals with digits."""

    def repl(match: re.Match[str]) -> str:
        token = match.group(0)
        return ROMAN_NUMERALS[token]

    pattern = re.compile(r"\b(ii|iii|iv)\b")
    return pattern.sub(repl, text)


def extract_parenthetical(text: str) -> Tuple[str, List[str]]:
    """Extract text within brackets into a hint field."""
    hints: List[str] = []

    def repl(match: re.Match[str]) -> str:
        # match groups correspond to (), [], {}
        for group in match.groups():
            if group:
                hints.append(group.strip())
        return ""

    text = PAREN_RE.sub(repl, text)
    return text, hints


def pretoken_cleanup(text: str) -> str:
    """Glue separated tokens before splitting."""
    text = LETTER_NUM_SPACE_RE.sub(r"\1\2", text)
    text = NUM_LETTER_SPACE_RE.sub(r"\1\2", text)
    text = HYPHEN_SPACE_RE.sub("-", text)
    return text


def tokenize(text: str) -> List[str]:
    """Split text into tokens using configured delimiters."""
    return [t for t in TOKEN_SPLIT_RE.split(text) if t]


def remove_weak_words(tokens: Sequence[str]) -> Tuple[List[str], List[str]]:
    """Remove weak/stop words from token list."""
    dropped: List[str] = []
    result: List[str] = []
    for tok in tokens:
        if tok in STOP_WORDS:
            dropped.append(tok)
        else:
            result.append(tok)
    return result, dropped


def apply_receptor_rules(text: str) -> Tuple[str, List[str], List[str]]:
    """Apply receptor family normalization rules.

    Returns
    -------
    Tuple of (new_text, gene_like_candidates, rules_applied)
    """
    candidates: List[str] = []
    applied: List[str] = []
    for pattern, replacement, candidate in RECEPTOR_RULES:
        if pattern.search(text):
            text = pattern.sub(replacement, text)
            candidates.append(candidate)
            applied.append(pattern.pattern)
    return text, candidates, applied


def generate_candidates(tokens: Sequence[str]) -> List[str]:
    """Generate gene-like candidates from tokens."""
    candidates: List[str] = []
    for tok in tokens:
        if tok in CANDIDATE_PREFIXES:
            candidates.append(CANDIDATE_PREFIXES[tok])
    return candidates


def final_cleanup(tokens: Sequence[str]) -> List[str]:
    """Remove duplicate tokens while preserving order."""
    seen = set()
    result: List[str] = []
    for tok in tokens:
        if tok not in seen:
            seen.add(tok)
            result.append(tok)
    return result


@dataclass
class NormalizationResult:
    raw: str
    clean_text: str
    query_tokens: List[str]
    gene_like_candidates: List[str]
    hint_taxon: int
    hints: Dict[str, List[str]]
    rules_applied: List[str]


def normalize_target_name(name: str) -> NormalizationResult:
    """Normalize a single target name.

    Parameters
    ----------
    name:
        Raw target name.

    Returns
    -------
    NormalizationResult
        Structured normalization information.
    """
    raw = name
    stage = sanitize_text(name)
    stage = normalize_unicode(stage)
    stage = replace_specials(stage)
    stage = replace_roman_numerals(stage)
    stage, parenthetical = extract_parenthetical(stage)
    stage = pretoken_cleanup(stage)
    stage, rule_candidates, rules_applied = apply_receptor_rules(stage)
    tokens = tokenize(stage)
    tokens, dropped = remove_weak_words(tokens)
    candidates = generate_candidates(tokens)
    candidates.extend(rule_candidates)
    tokens = final_cleanup(tokens)
    clean_text = " ".join(tokens)
    hints = {"parenthetical": parenthetical, "dropped": dropped}
    return NormalizationResult(
        raw=raw,
        clean_text=clean_text,
        query_tokens=tokens,
        gene_like_candidates=final_cleanup(candidates),
        hint_taxon=9606,
        hints=hints,
        rules_applied=rules_applied,
    )
