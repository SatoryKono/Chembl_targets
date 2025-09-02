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
HYPHEN_SPACE_RE = re.compile(r"\s*-\s*")
HYPHEN_TOKEN_RE = re.compile(r"\b[a-z0-9]+(?:-[a-z0-9]+)+\b")
SHORT_TOKEN_RE = re.compile(r"^[a-z0-9]{1,3}$")
INDEX_TOKEN_RE = re.compile(r"^(?:[a-z]\d(?:[a-z]\d+)?|5-?ht\d+[a-z]?)$")
LETTER_DIGIT_SPLIT_RE = re.compile(r"\b(?:[a-z]\s+\d+|\d+\s+[a-z])\b")

# Patterns for mutation extraction ------------------------------------------------
MUTATION_PATTERNS: Sequence[re.Pattern[str]] = (
    re.compile(r"p\.[A-Z][0-9]+[A-Z]", re.IGNORECASE),
    re.compile(r"p\.[A-Z][0-9]+(?:\*|Ter)", re.IGNORECASE),
    re.compile(r"p\.[A-Z][0-9]+(?:_[A-Z][0-9]+)?del", re.IGNORECASE),
    re.compile(r"p\.[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+", re.IGNORECASE),
    re.compile(r"p\.[A-Z][0-9]+(?:_[A-Z][0-9]+)?dup", re.IGNORECASE),
    re.compile(r"p\.[A-Z][0-9]+fs(?:\*[0-9]+)?", re.IGNORECASE),
    re.compile(r"p\.Met1\?", re.IGNORECASE),
    re.compile(r"p\.\*[0-9]+[A-Z]", re.IGNORECASE),
    re.compile(r"p\.[A-Z][0-9]+(?:_[A-Z][0-9]+)?delins[A-Z]+", re.IGNORECASE),
    re.compile(r"p\.[A-Z][a-z]{2}[0-9]+(?:[A-Z][a-z]{2}|\*|Ter)", re.IGNORECASE),
    re.compile(r"\b[A-Z][0-9]+[A-Z]\b", re.IGNORECASE),
    re.compile(
        r"(?<!p\.)[A-Z][a-z]{2}[0-9]+(?:[A-Z][a-z]{2}|\*|Ter)\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\b[pcgnmr]\.[0-9]+[+-]?[0-9]*(?:_[+-]?[0-9]+)?(?:[ACGT]>[ACGT]|delins|del|ins|dup|inv|fs\*?[0-9]*)\b",
        re.IGNORECASE,
    ),
    re.compile(r"\b(?:[A-Z][0-9]+[A-Z])(?:/[A-Z][0-9]+[A-Z])+\b", re.IGNORECASE),
    re.compile(r"(?:Δ|delta)\s?[A-Z][0-9]+", re.IGNORECASE),
    re.compile(r"\b(mutant|variant|mut\.)\b", re.IGNORECASE),
)

MUTATION_WHITELIST = {
    "m2",
    "h3",
    "d2",
    "p2x7",
    "p2x",
    "5-ht1a",
    "alpha1",
    "beta2",
}


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
                keep_tokens.append(token.strip())
        return ""

    text = PAREN_RE.sub(repl, text)
    return text, hints, keep_tokens


def pretoken_cleanup(text: str) -> str:
    """Normalize spaces around hyphens before splitting."""
    text = HYPHEN_SPACE_RE.sub("-", text)
    return text


def generate_letter_digit_variants(tokens: Sequence[str]) -> List[Tuple[str, str]]:
    """Create joined variants for adjacent letter-digit tokens.

    Parameters
    ----------
    tokens:
        Tokens obtained after initial splitting.

    Returns
    -------
    List[Tuple[str, str]]
        Pairs of (variant, base_pattern) where ``base_pattern`` is the
        space-separated form in the original text.
    """

    variants: List[Tuple[str, str]] = []
    for i in range(len(tokens) - 1):
        left, right = tokens[i], tokens[i + 1]
        if left.isalpha() and right.isdigit():
            base = f"{left} {right}"
            variants.append((f"{left}{right}", base))
            variants.append((f"{left}-{right}", base))
    return variants


def detect_hyphen_variants(text: str) -> List[Tuple[str, str]]:
    """Return hyphenated tokens and their space-separated base pattern.

    Each result is a pair of (variant, base_pattern). Variants include the
    original hyphenated token and its concatenated counterpart.
    """

    variants: List[Tuple[str, str]] = []
    for token in HYPHEN_TOKEN_RE.findall(text):
        base = token.replace("-", " ")
        variants.append((token, base))
        variants.append((token.replace("-", ""), base))
    return variants


def build_variant_strings(
    base: str,
    substitutions: Sequence[Tuple[str, str]],
    extra: Sequence[str] | None = None,
) -> List[str]:
    """Generate variant strings from a base text and substitution patterns.

    Parameters
    ----------
    base:
        Base string joined from primary tokens.
    substitutions:
        Sequence of ``(variant, base_pattern)`` tuples used to replace parts of
        the base string.
    extra:
        Additional standalone tokens to include as separate variants.

    Returns
    -------
    List[str]
        Unique variant strings with empty values removed.
    """

    variants: List[str] = []
    base = base.strip()
    if base and not LETTER_DIGIT_SPLIT_RE.search(base):
        variants.append(base)
    for var, pattern in substitutions:
        if base and pattern in base:
            variants.append(base.replace(pattern, var))
        variants.append(var)
    if extra:
        variants.extend(extra)
    seen: List[str] = []
    for v in variants:
        v = v.strip()
        if v and v not in seen:
            seen.append(v)
    return seen


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


def find_mutations(text: str) -> List[str]:
    """Extract mutation-like substrings from text.

    Parameters
    ----------
    text:
        Input string in its original form.

    Returns
    -------
    List[str]
        Unique mutation substrings in order of appearance.
    """

    found: List[str] = []
    for pattern in MUTATION_PATTERNS:
        for match in pattern.finditer(text):
            token = match.group(0)
            lower = token.lower()
            if lower in MUTATION_WHITELIST:
                continue
            if any(lower in f.lower() for f in found):
                continue
            found = [f for f in found if f.lower() not in lower]
            if token not in found:
                found.append(token)
    return found


def mutation_token_set(mutations: Sequence[str]) -> set[str]:
    """Normalize and tokenize mutation strings for removal.

    Parameters
    ----------
    mutations:
        Mutation substrings captured from the raw text.

    Returns
    -------
    set[str]
        Lowercase tokens representing mutations to exclude from results.
    """

    tokens: set[str] = set()
    for mut in mutations:
        norm = normalize_unicode(mut)
        norm = replace_specials(norm)
        norm = replace_roman_numerals(norm)
        norm_tokens = tokenize(norm)
        for tok in norm_tokens:
            tokens.add(tok)
    return tokens


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
    hints: Dict[str, List[str] | bool]
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
    mutations = find_mutations(stage)
    stage = normalize_unicode(stage)
    stage = replace_specials(stage)
    stage = replace_roman_numerals(stage)
    stage, parenthetical, paren_tokens = extract_parenthetical(stage)
    if paren_tokens:
        stage = f"{stage} {' '.join(paren_tokens)}".strip()
    stage = pretoken_cleanup(stage)
    stage, rule_candidates, rules_applied = apply_receptor_rules(stage)
    hyphen_subs = detect_hyphen_variants(stage)
    tokens_base = tokenize(stage)
    letter_digit_subs = generate_letter_digit_variants(tokens_base)
    substitutions = hyphen_subs + letter_digit_subs
    tokens_base_no_stop, dropped = remove_weak_words(tokens_base)
    tokens_base_alt = final_cleanup(tokens_base)
    tokens_no_stop = tokens_base_no_stop + [v for v, _ in substitutions]
    tokens_alt = tokens_base_alt + [v for v, _ in substitutions]

    mutation_tokens = mutation_token_set(mutations)
    tokens_no_stop_orig = list(tokens_no_stop)
    tokens_alt_orig = list(tokens_alt)
    tokens_base_no_stop_orig = list(tokens_base_no_stop)
    tokens_base_alt_orig = list(tokens_base_alt)

    tokens_no_stop = [
        t for t in tokens_no_stop if t not in mutation_tokens or t in MUTATION_WHITELIST
    ]
    tokens_alt = [
        t for t in tokens_alt if t not in mutation_tokens or t in MUTATION_WHITELIST
    ]
    tokens_base_no_stop = [
        t
        for t in tokens_base_no_stop
        if t not in mutation_tokens or t in MUTATION_WHITELIST
    ]
    tokens_base_alt = [
        t
        for t in tokens_base_alt
        if t not in mutation_tokens or t in MUTATION_WHITELIST
    ]

    clean_tokens_check = [
        t for t in tokens_no_stop if not re.fullmatch(r"[a-z]$|\d+$", t)
    ]
    hints_mutations_only = False
    if not clean_tokens_check:
        tokens_no_stop = tokens_no_stop_orig
        tokens_alt = tokens_alt_orig
        tokens_base_no_stop = tokens_base_no_stop_orig
        tokens_base_alt = tokens_base_alt_orig
        hints_mutations_only = True

    tokens_no_stop = final_cleanup(tokens_no_stop)
    tokens_alt = final_cleanup(tokens_alt)
    base_no_stop_str = " ".join(tokens_base_no_stop)
    base_alt_str = " ".join(tokens_base_alt)

    clean_variants = build_variant_strings(
        base_no_stop_str, substitutions, paren_tokens
    )
    clean_variants_alt = build_variant_strings(
        base_alt_str, substitutions, paren_tokens
    )

    clean_text = (
        "|".join(clean_variants)
        if len(clean_variants) > 1
        else (clean_variants[0] if clean_variants else "")
    )
    clean_text_alt = (
        "|".join(clean_variants_alt)
        if len(clean_variants_alt) > 1
        else (clean_variants_alt[0] if clean_variants_alt else "")
    )

    text_for_candidates = " ".join(tokens_no_stop)
    regex_candidates = generate_regex_candidates(text_for_candidates)
    candidates = rule_candidates + regex_candidates
    hints: Dict[str, List[str] | bool] = {
        "parenthetical": parenthetical,
        "dropped": dropped,
        "mutations": mutations,
    }
    if hints_mutations_only:
        hints["mutations_only"] = True
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
