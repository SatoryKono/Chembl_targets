"""Text normalization transformations for protein target names."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
import unicodedata
from typing import Callable, Dict, List, Sequence, Tuple

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

_DEFAULT_SPECIALS_TABLE = str.maketrans({**GREEK_LETTERS, **SUPERSCRIPTS})

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

GENE_CANDIDATE_RULES: Sequence[Tuple[re.Pattern[str], str | Sequence[str]]] = (
    # Histamine, dopamine, adrenergic, etc.
    (re.compile(r"histamine\s+h(\d+)"), r"hrh\1"),
    (re.compile(r"dopamine\s+d(\d+)"), r"drd\1"),
    (re.compile(r"adrenergic\s+beta(\d+)"), r"adrb\1"),
    (re.compile(r"p2x(\d+)"), r"p2rx\1"),
    (re.compile(r"5[- ]?ht(\d+[a-z]?)"), r"htr\1"),
    (re.compile(r"gaba\s+a\s+alpha(\d+)"), r"gabra\1"),
    # TRP channels (trpv/trpm/trpc/trpa/trpk)
    (re.compile(r"trp\s*([vmcak])\s*(\d+)"), r"trp\1\2"),
    # Ionotropic glutamate receptors
    (re.compile(r"glua(\d)"), r"gria\1"),
    (re.compile(r"gluk(\d)"), r"grik\1"),
    (re.compile(r"nr(1|2[a-d]|3[a-b])"), r"grin\1"),
    # Metabotropic glutamate receptors
    (re.compile(r"mglur(\d)"), r"grm\1"),
    # Chemokine receptors and full forms
    (re.compile(r"ccr\s*(\d+)"), r"ccr\1"),
    (re.compile(r"cxcr\s*(\d+)"), r"cxcr\1"),
    (re.compile(r"chemokine\s+cc\s*(\d+)"), r"ccr\1"),
    (re.compile(r"chemokine\s+cxc\s*(\d+)"), r"cxcr\1"),
)

ALIAS_CANDIDATE_RULES: Sequence[Tuple[re.Pattern[str], Sequence[str]]] = (
    (re.compile(r"sdf[- ]?1"), ["cxcr4"]),
    (re.compile(r"il[- ]?8"), ["cxcr1", "cxcr2"]),
    (re.compile(r"rantes"), ["ccr1", "ccr3", "ccr5"]),
    (re.compile(r"fractalkine"), ["cx3cr1"]),
)

# GPCR-specific candidate rules -------------------------------------------------

RULES_GPCR: Sequence[
    Tuple[re.Pattern[str], Sequence[str] | Callable[[re.Match[str]], Sequence[str]]]
] = [
    # --- Adenosine ---
    (
        re.compile(r"\badenosine\s+receptor\b"),
        ["adora1", "adora2a", "adora2b", "adora3"],
    ),
    (
        re.compile(r"\badenosine\s+a\s*1\b|\ba1\b(?=.*adenosine)"),
        ["a1", "adora1"],
    ),
    (
        re.compile(r"\badenosine\s+a\s*2\s*a\b|\ba2a\b(?=.*adenosine)"),
        ["a2a", "adora2a"],
    ),
    (
        re.compile(r"\badenosine\s+a\s*2\s*b\b|\ba2b\b(?=.*adenosine)"),
        ["a2b", "adora2b"],
    ),
    (
        re.compile(r"\badenosine\s+a\s*3\b|\ba3\b(?=.*adenosine)"),
        ["a3", "adora3"],
    ),
    # --- Nociceptin / ORL1 ---
    (
        re.compile(
            r"\bnociceptin\s+receptor\b|\borphanin\s*fq\s+receptor\b|\bnop\b|\borl1\b"
        ),
        ["nop", "orl1", "oprl1"],
    ),
    # --- Neuropeptide Y (NPY) ---
    (
        re.compile(r"\bneuropeptide\s*y\s+receptor\b|\bnpy\s+receptor\b"),
        ["npy1r", "npy2r", "npy4r", "npy5r"],
    ),
    (re.compile(r"\b(y\s*1|npy\s*1)\b"), ["y1", "npy1r"]),
    (re.compile(r"\b(y\s*2|npy\s*2)\b"), ["y2", "npy2r"]),
    (re.compile(r"\b(y\s*4|npy\s*4)\b"), ["y4", "npy4r"]),
    (re.compile(r"\b(y\s*5|npy\s*5)\b"), ["y5", "npy5r"]),
    # --- Melanocortin (MC) ---
    (
        re.compile(r"\bmelanocortin\s+receptor\b|\bmcr\b"),
        ["mc1r", "mc2r", "mc3r", "mc4r", "mc5r"],
    ),
    (
        re.compile(r"\bmc\s*([1-5])\s*r?\b|\bmelanocortin\s*-\s*([1-5])\s*receptor\b"),
        lambda m: [f"mc{(m.group(1) or m.group(2))}r"],
    ),
    # --- Prostaglandin (EP/DP/FP/IP/TP) ---
    (
        re.compile(r"\bprostaglandin\s+receptor\b"),
        [
            "ptger1",
            "ptger2",
            "ptger3",
            "ptger4",
            "ptgdr",
            "ptgdr2",
            "ptgfr",
            "ptgir",
            "tbxa2r",
        ],
    ),
    # EP1-EP4
    (
        re.compile(r"\bep\s*([1-4])\b"),
        lambda m: [f"ep{m.group(1)}", f"ptger{m.group(1)}"],
    ),
    # DP1, DP2/CRTH2
    (re.compile(r"\bdp\s*1\b"), ["dp1", "ptgdr"]),
    (re.compile(r"\bdp\s*2\b|\bcrth2\b|\bgpr44\b"), ["dp2", "crth2", "ptgdr2"]),
    # FP / PGF2A receptor
    (re.compile(r"\bfp\b|\bpgf\s*2\s*a\b"), ["fp", "ptgfr"]),
    # IP
    (
        re.compile(r"\bip\b(?!(v|3))|\bprostacyclin\s+receptor\b"),
        ["ip", "ptgir"],
    ),
    # TP / thromboxane receptor
    (re.compile(r"\btp\b|\bthromboxane\s+receptor\b"), ["tp", "tbxa2r"]),
]

# Additional GPCR-specific rules ------------------------------------------------

RULES_GPCR_EXTRA: Sequence[
    Tuple[re.Pattern[str], Sequence[str] | Callable[[re.Match[str]], Sequence[str]]]
] = [
    # 1) Calcitonin/CGRP/Amylin (CALCR, CALCRL + RAMP)
    (
        re.compile(r"\b(calcitonin|cgrp|amylin)\s+receptor\b(?!\s*\d)|\bcalcrl?\b"),
        ["calcr", "calcrl", "ramp1", "ramp2", "ramp3"],
    ),
    (re.compile(r"\bcgrp\b"), ["calcrl", "ramp1", "ramp2", "ramp3"]),
    (re.compile(r"\bamylin\b"), ["calcr", "ramp1", "ramp2", "ramp3"]),
    # 2) Parathyroid hormone (PTH1R/PTH2R)
    (
        re.compile(
            r"\bparathyroid\s+hormone\s+receptor\b(?!\s*\d)|\bpth\s*receptor\b(?!\s*\d)"
        ),
        ["pth1r", "pth2r"],
    ),
    (re.compile(r"\bpth\s*1\s*r?\b"), ["pth1r"]),
    (re.compile(r"\bpth\s*2\s*r?\b"), ["pth2r"]),
    # 3) Neuropeptide S (NPS)
    (
        re.compile(
            r"\bneuropeptide\s*s\s+receptor\b(?!\s*\d)|\bnps\s*receptor\b(?!\s*\d)|\bnpsr1\b"
        ),
        ["npsr1"],
    ),
    # 4) Neuropeptide FF (NPFF)
    (
        re.compile(r"\bneuropeptide\s*ff\s+receptor\b(?!\s*\d)|\bnpffr\b"),
        ["npffr1", "npffr2"],
    ),
    (re.compile(r"\bnpffr\s*([12])\b"), lambda m: [f"npffr{m.group(1)}"]),
    # 5) Neuropeptide B/W (NPB/W)
    (
        re.compile(r"\bneuropeptide\s*(b|w)\s+receptor\b(?!\s*\d)|\bnpbwr\b"),
        ["npbwr1", "npbwr2"],
    ),
    (re.compile(r"\bnpbwr\s*([12])\b"), lambda m: [f"npbwr{m.group(1)}"]),
    # 6) Neuromedin U (NMU)
    (
        re.compile(r"\bneuromedin\s*u\s+receptor\b(?!\s*\d)|\bnmur\b"),
        ["nmur1", "nmur2"],
    ),
    (re.compile(r"\bnmur\s*([12])\b"), lambda m: [f"nmur{m.group(1)}"]),
    # 7) Kisspeptin (KISS1R)
    (
        re.compile(r"\bkisspeptin\s+receptor\b|\bgpr54\b|\bkiss1r\b"),
        ["kiss1r", "gpr54"],
    ),
    # 8) Ghrelin (GHSR)
    (re.compile(r"\bghrelin\s+receptor\b|\bghsr\b"), ["ghsr"]),
    # 9) Motilin (MLNR)
    (re.compile(r"\bmotilin\s+receptor\b|\bmlnr\b|\bgpr38\b"), ["mlnr", "gpr38"]),
    # 10) Prolactin-releasing peptide (PRLHR)
    (
        re.compile(
            r"\bprolactin-?releasing\s+peptide\s+receptor\b|\bprlhr\b|\bgpr10\b"
        ),
        ["prlhr", "gpr10"],
    ),
    # 11) Melanin-concentrating hormone (MCHR1/2)
    (
        re.compile(
            r"\bmelanin-?concentrating\s+hormone\s+receptor\b(?!\s*\d)|\bmchr\b"
        ),
        ["mchr1", "mchr2"],
    ),
    (re.compile(r"\bmchr\s*([12])\b"), lambda m: [f"mchr{m.group(1)}"]),
    # 12) Fractalkine (CX3CR1) и XCR1
    (re.compile(r"\bfractalkine\s+receptor\b(?!\s*\d)|\bcx3cr1\b"), ["cx3cr1"]),
    (
        re.compile(r"\bxcr\s*1\b|\bxcr1\b|\b(xc)\s*chemokine\s+receptor\s*1\b"),
        ["xcr1"],
    ),
    # 13) Platelet-activating factor (PTAFR)
    (
        re.compile(r"\bplatelet-?activating\s+factor\s+receptor\b(?!\s*\d)|\bptafr\b"),
        ["ptafr"],
    ),
    # 14) Formyl peptide receptors (FPR1-3)
    (
        re.compile(r"\bformyl\s+peptide\s+receptor\b(?!\s*\d)|\bfpr\b"),
        ["fpr1", "fpr2", "fpr3"],
    ),
    (re.compile(r"\bfpr\s*([1-3])\b"), lambda m: [f"fpr{m.group(1)}"]),
    (re.compile(r"\balx\b"), ["fpr2"]),  # FPR2/ALX
    # 15) Free fatty acid receptors (FFAR/GPR40/41/43/120 + GPR84)
    (
        re.compile(r"\bfree\s+fatty\s+acid\s+receptor\b(?!\s*\d)|\bffar\b"),
        ["ffar1", "ffar2", "ffar3", "ffar4", "gpr84"],
    ),
    (re.compile(r"\bffar\s*([1-4])\b"), lambda m: [f"ffar{m.group(1)}"]),
    (re.compile(r"\bgpr\s*120\b"), ["ffar4", "gpr120"]),
    (re.compile(r"\bgpr\s*40\b"), ["ffar1", "gpr40"]),
    (re.compile(r"\bgpr\s*41\b"), ["ffar3", "gpr41"]),
    (re.compile(r"\bgpr\s*43\b"), ["ffar2", "gpr43"]),
    (re.compile(r"\bgpr\s*84\b"), ["gpr84"]),
    # 16) Hydroxycarboxylic acid (HCAR/GPR81/109A/109B)
    (
        re.compile(r"\bhydroxycarboxylic\s+acid\s+receptor\b(?!\s*\d)|\bhcar\b"),
        ["hcar1", "hcar2", "hcar3"],
    ),
    (re.compile(r"\bhcar\s*([1-3])\b"), lambda m: [f"hcar{m.group(1)}"]),
    (re.compile(r"\bgpr\s*81\b"), ["hcar1", "gpr81"]),
    (re.compile(r"\bgpr\s*109\s*a\b|\bhcar2\b"), ["hcar2", "gpr109a"]),
    (re.compile(r"\bgpr\s*109\s*b\b|\bhcar3\b"), ["hcar3", "gpr109b"]),
    # 17) Trace amine-associated receptors (TAAR)
    (
        re.compile(r"\btrace\s+amine-?associated\s+receptor\b(?!\s*\d)|\btaar\b"),
        [
            "taar1",
            "taar2",
            "taar3",
            "taar4",
            "taar5",
            "taar6",
            "taar7",
            "taar8",
            "taar9",
        ],
    ),
    (re.compile(r"\btaar\s*([1-9])\b"), lambda m: [f"taar{m.group(1)}"]),
    # 18) Bile acid receptor (TGR5/GPBAR1)
    (
        re.compile(r"\bbile\s+acid\s+receptor\b(?!\s*\d)|\btgr5\b|\bgpbar1\b"),
        ["gpbar1", "tgr5"],
    ),
    # 19) Urotensin II receptor (UTS2R)
    (re.compile(r"\burotensin\s+ii\s+receptor\b|\buts2r\b"), ["uts2r"]),
    # 20) Apelin receptor (APLNR)
    (re.compile(r"\bapelin\s+receptor\b|\baplnr\b|\bagtrl1\b"), ["aplnr"]),
]


def gen_candidates(clean_text: str) -> List[str]:
    """Generate GPCR gene-like candidates based on ``clean_text``.

    Parameters
    ----------
    clean_text:
        Normalized text (ideally with stop words) used for pattern matching.

    Returns
    -------
    list[str]
        Candidate gene symbols inferred from GPCR-specific rules. Duplicates
        are removed while preserving the order of first occurrence.
    """

    s = clean_text.lower()
    out: List[str] = []
    for pat, val in list(RULES_GPCR) + list(RULES_GPCR_EXTRA):
        m = pat.search(s)
        if m:
            adds = val(m) if callable(val) else val
            out.extend(adds)

    # dedupe preserving order
    seen: set[str] = set()
    res: List[str] = []
    for x in out:
        if x and x not in seen:
            seen.add(x)
            res.append(x)
    return res


ROMAN_NUMERALS: Dict[str, str] = {
    "ii": "2",
    "iii": "3",
    "iv": "4",
    "vi": "6",
    "vii": "7",
    "viii": "8",
    "ix": "9",
}

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

# Pattern capturing simple letter-number-letter mutations like "V600E"
LETTER_DIGIT_LETTER_RE = re.compile(r"\b([A-Z])(\d+)([A-Z])\b", re.IGNORECASE)

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
    LETTER_DIGIT_LETTER_RE,
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
    """Replace Greek letters and superscripts using translation tables.

    Parameters
    ----------
    text:
        Input string to process.
    greek:
        Optional mapping of Greek characters to ASCII strings. If ``None``,
        :data:`GREEK_LETTERS` is used.
    superscripts:
        Optional mapping of superscript/subscript digits to ASCII digits. If
        ``None``, :data:`SUPERSCRIPTS` is used.

    Notes
    -----
    Characters are replaced using :meth:`str.translate` for efficiency. Unicode
    normalization should be applied beforehand to collapse variants such as the
    micro sign ``µ`` into ``μ``.
    """
    if greek is None and superscripts is None:
        translation = _DEFAULT_SPECIALS_TABLE
    else:
        translation = str.maketrans({**(greek or GREEK_LETTERS), **(superscripts or SUPERSCRIPTS)})
    return text.translate(translation)


def replace_roman_numerals(text: str) -> str:
    """Replace standalone Roman numerals with digits.

    The mapping covers numerals from II to IX. Single-letter numerals such as
    ``v`` or ``x`` are intentionally excluded to avoid altering valid gene
    symbols. The input should already be lower-cased.
    """

    def repl(match: re.Match[str]) -> str:
        token = match.group(0)
        return ROMAN_NUMERALS[token]

    pattern = re.compile(
        r"\b(" + "|".join(sorted(ROMAN_NUMERALS.keys(), key=len, reverse=True)) + r")\b"
    )
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
            # Skip letter-digit-letter where the letters are identical (e.g., A123A)
            if pattern is LETTER_DIGIT_LETTER_RE:
                if match.group(1).upper() == match.group(3).upper():
                    continue
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
    """Generate gene-like candidates from regex and alias rules.

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
    tokens = text.split()

    def apply_rules(
        rules: Sequence[Tuple[re.Pattern[str], str | Sequence[str]]],
    ) -> None:
        for pattern, repl in rules:
            for match in pattern.finditer(text):
                if isinstance(repl, str):
                    cand = match.expand(repl).lower()
                    if cand not in candidates:
                        candidates.append(cand)
                else:
                    for cand in repl:
                        cand = cand.lower()
                        if cand not in candidates:
                            candidates.append(cand)

    apply_rules(GENE_CANDIDATE_RULES)

    # Family-level inference when subtype tokens are absent
    if "ampa" in tokens and not any(t.startswith("glua") for t in tokens):
        for i in range(1, 5):
            cand = f"gria{i}"
            if cand not in candidates:
                candidates.append(cand)
    if "nmda" in tokens and not any(t.startswith("nr") for t in tokens):
        for cand in [
            "grin1",
            "grin2a",
            "grin2b",
            "grin2c",
            "grin2d",
            "grin3a",
            "grin3b",
        ]:
            if cand not in candidates:
                candidates.append(cand)
    if "kainate" in tokens and not any(t.startswith("gluk") for t in tokens):
        for i in range(1, 6):
            cand = f"grik{i}"
            if cand not in candidates:
                candidates.append(cand)
    if (
        "metabotropic" in tokens
        and "glutamate" in tokens
        and not any(t.startswith("mglur") for t in tokens)
    ):
        for i in range(1, 9):
            cand = f"grm{i}"
            if cand not in candidates:
                candidates.append(cand)

    apply_rules(ALIAS_CANDIDATE_RULES)

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
    gpcr_candidates = gen_candidates(clean_text_alt or clean_text)
    regex_candidates = generate_regex_candidates(text_for_candidates)
    candidates = rule_candidates + gpcr_candidates + regex_candidates
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
