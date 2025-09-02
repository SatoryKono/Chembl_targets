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
    "channel",
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
CANDIDATE_REGEX_RULES: Sequence[Tuple[re.Pattern[str], str]] = (
    (re.compile(r"histamine\s+h(\d+)"), r"hrh\1"),
    (re.compile(r"dopamine\s+d(\d+)"), r"drd\1"),
    (re.compile(r"adrenergic\s+beta(\d+)"), r"adrb\1"),
    (re.compile(r"p2x(\d+)"), r"p2rx\1"),
    (re.compile(r"5[- ]?ht(\d+[a-z]?)"), r"htr\1"),
    (re.compile(r"gaba\s+a\s+alpha(\d+)"), r"gabra\1"),
)

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
WORD_NUM_SPACE_RE = re.compile(r"(?<=\b)([a-z]+)\s+([0-9]+)(?=\b)", re.I)
HYPHEN_TOKEN_RE = re.compile(r"\b[a-z0-9]+(?:-[a-z0-9]+)+\b")
SHORT_TOKEN_RE = re.compile(r"^[a-z0-9]{1,3}$")
INDEX_TOKEN_RE = re.compile(r"^(?:[a-z]\d|5-?ht\d+[a-z]?)$")


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


def extract_parenthetical(text: str) -> Tuple[str, List[str], List[str]]:
    """Extract bracketed text into hints and retain certain short tokens.

    Parameters
    ----------
    text:
        Input string with possible parenthetical segments.

    Returns
    -------
    Tuple[str, List[str], List[str]]
        The text with brackets removed, list of extracted strings for hints,
        and tokens that should remain in the main text.
    """

    hints: List[str] = []
    keep_tokens: List[str] = []

    def repl(match: re.Match[str]) -> str:
        for group in match.groups():
            if not group:
                continue
            token = group.strip()
            hints.append(token)
            compact = re.sub(r"[\s_\-]", "", token)
            if SHORT_TOKEN_RE.fullmatch(compact) or INDEX_TOKEN_RE.fullmatch(compact):
                keep_tokens.append(compact)
        return ""

    text = PAREN_RE.sub(repl, text)
    return text, hints, keep_tokens


def pretoken_cleanup(text: str) -> str:
    """Glue separated tokens before splitting."""
    text = LETTER_NUM_SPACE_RE.sub(r"\1\2", text)
    text = NUM_LETTER_SPACE_RE.sub(r"\1\2", text)
    text = HYPHEN_SPACE_RE.sub("-", text)
    return text


def detect_space_variants(text: str) -> List[str]:
    """Return variants for letter-number pairs separated by space."""

    variants: List[str] = []
    for letters, digits in WORD_NUM_SPACE_RE.findall(text):
        variants.append(f"{letters}-{digits}")
        variants.append(f"{letters}{digits}")
    return variants


def detect_hyphen_variants(text: str) -> List[str]:
    """Return tokens with hyphen and their concatenated counterparts."""

    variants: List[str] = []
    for token in HYPHEN_TOKEN_RE.findall(text):
        variants.append(token)
        variants.append(token.replace("-", ""))
    return variants


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


def generate_regex_candidates(text: str) -> List[str]:
    """Generate gene-like candidates from regex rules.

    Parameters
    ----------
    text:
        Concatenated normalized tokens.

    Returns
    -------
    List[str]
        Candidate gene symbols inferred from the text.
    """

    candidates: List[str] = []
    for pattern, repl in CANDIDATE_REGEX_RULES:
        for match in pattern.finditer(text):
            candidate = match.expand(repl)
            if candidate not in candidates:
                candidates.append(candidate)
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
    clean_text_alt: str
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
    stage, parenthetical, paren_tokens = extract_parenthetical(stage)
    if paren_tokens:
        stage = f"{stage} {' '.join(paren_tokens)}".strip()
    space_variants = detect_space_variants(stage)
    stage = pretoken_cleanup(stage)
    stage, rule_candidates, rules_applied = apply_receptor_rules(stage)
    hyphen_variants = detect_hyphen_variants(stage)
    tokens_raw = tokenize(stage)
    tokens_raw.extend(space_variants + hyphen_variants)
    tokens_no_stop, dropped = remove_weak_words(tokens_raw)
    tokens_no_stop = final_cleanup(tokens_no_stop)
    tokens_alt = final_cleanup(tokens_raw)
    text_for_candidates = " ".join(tokens_no_stop)
    regex_candidates = generate_regex_candidates(text_for_candidates)
    candidates = rule_candidates + regex_candidates
    clean_text = text_for_candidates
    clean_text_alt = " ".join(tokens_alt)
    hints = {"parenthetical": parenthetical, "dropped": dropped}
    return NormalizationResult(
        raw=raw,
        clean_text=clean_text,
        clean_text_alt=clean_text_alt,
        query_tokens=tokens_no_stop,
        gene_like_candidates=final_cleanup(candidates),
        hint_taxon=9606,
        hints=hints,
        rules_applied=rules_applied,
    )
