"""Text normalization transformations for protein target names."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import re
import unicodedata
from typing import Callable, Dict, List, Sequence, Tuple

try:
    from hgvs.exceptions import HGVSParseError  # type: ignore[import-not-found]
    from hgvs.parser import Parser as _HgvsParser  # type: ignore[import-not-found]

    _HGVS_PARSER: _HgvsParser | None = _HgvsParser()
except Exception:  # pragma: no cover - optional dependency
    _HGVS_PARSER = None
    HGVSParseError = Exception

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
    "µ": "mu",  # micro sign
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

# ---------------------------------------------------------------------------
# Additional candidate rules supplied via domain-specific mapping

CandidateRule = Tuple[
    re.Pattern[str],
    Sequence[str] | Callable[[re.Match[str]], Sequence[str]],
    bool,
]

# Flexible separator "word - word" etc.
CUSTOM_SP = r"\s*[- ]\s*"

RULES_CUSTOM: Sequence[CandidateRule] = [
    # ===== HISTAMINE =====
    # H3
    (
        re.compile(
            rf"\b(h3r|h3{CUSTOM_SP}receptor|histamine{CUSTOM_SP}3{CUSTOM_SP}receptor)\b"
        ),
        ["h3r", "hrh3"],
        False,
    ),
    # H4
    (
        re.compile(
            rf"\b(h4r|h4{CUSTOM_SP}receptor|histamine{CUSTOM_SP}4{CUSTOM_SP}receptor)\b"
        ),
        ["h4r", "hrh4"],
        False,
    ),
    # ===== PROSTAGLANDIN / CRTH2 =====
    (
        re.compile(
            rf"\b(chemoattractant{CUSTOM_SP}receptor{CUSTOM_SP}homologous{CUSTOM_SP}molecule|crth2|gpr44|dp2|prostaglandin{CUSTOM_SP}d2{CUSTOM_SP}receptor\b.*(type|2)?)\b"
        ),
        ["crth2", "gpr44", "ptgdr2", "dp2"],
        False,
    ),
    (
        re.compile(
            rf"\b(prostanoid{CUSTOM_SP}dp{CUSTOM_SP}receptor|dp{CUSTOM_SP}receptor|dp\b(?!\d))\b"
        ),
        ["dp1", "ptgdr"],
        False,
    ),
    (
        re.compile(rf"\bprostaglandin{CUSTOM_SP}d2{CUSTOM_SP}receptor\b"),
        ["dp1", "ptgdr", "dp2", "ptgdr2"],
        False,
    ),
    # ===== NPFF =====
    (re.compile(rf"\bnpff2{CUSTOM_SP}receptor\b"), ["npffr2"], False),
    (re.compile(rf"\bnpff1{CUSTOM_SP}receptor\b"), ["npffr1"], False),
    # ===== SIGMA =====
    (
        re.compile(
            rf"\b(sigma{CUSTOM_SP}1{CUSTOM_SP}opioid{CUSTOM_SP}receptor|sigma1{CUSTOM_SP}opioid{CUSTOM_SP}receptor|sigma1{CUSTOM_SP}receptor|sigma{CUSTOM_SP}1{CUSTOM_SP}receptor|sigma-?1{CUSTOM_SP}receptor|sigma1-?receptor)\b"
        ),
        ["sigmar1"],
        False,
    ),
    (
        re.compile(
            rf"\b(sigma{CUSTOM_SP}2{CUSTOM_SP}receptor|sigma2{CUSTOM_SP}receptor|sigma-?2{CUSTOM_SP}receptor|sigma{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}2)\b"
        ),
        ["tmem97", "sigma2r"],
        False,
    ),
    (re.compile(rf"\bsigma{CUSTOM_SP}receptor\b"), ["sigmar1", "tmem97"], True),
    # ===== IMIDAZOLINE I1 =====
    (
        re.compile(
            rf"\b(imidazoline{CUSTOM_SP}i\s*1{CUSTOM_SP}receptor|imidazoline{CUSTOM_SP}receptor{CUSTOM_SP}1|i1{CUSTOM_SP}imidazoline{CUSTOM_SP}receptor|i1{CUSTOM_SP}receptor{CUSTOM_SP}imidazoline)\b"
        ),
        ["i1r", "nisch"],
        True,
    ),
    # ===== UROTENSIN II =====
    (re.compile(rf"\burotensin{CUSTOM_SP}2{CUSTOM_SP}receptor\b"), ["uts2r"], False),
    # ===== LEUKOTRIENE / CysLT =====
    (
        re.compile(
            rf"\b(cyslt2{CUSTOM_SP}receptor|cysteinyl{CUSTOM_SP}leukotriene{CUSTOM_SP}receptor{CUSTOM_SP}2)\b"
        ),
        ["cysltr2"],
        False,
    ),
    (
        re.compile(
            rf"\b(cyslt1{CUSTOM_SP}receptor|cysteinyl{CUSTOM_SP}leukotriene{CUSTOM_SP}receptor{CUSTOM_SP}1|leukotriene{CUSTOM_SP}d4{CUSTOM_SP}receptor|cysteinyl{CUSTOM_SP}leukotriene{CUSTOM_SP}d4{CUSTOM_SP}receptor|ltd4{CUSTOM_SP}receptor)\b"
        ),
        ["cysltr1"],
        False,
    ),
    (
        re.compile(
            rf"\b(leukotriene{CUSTOM_SP}b4{CUSTOM_SP}receptor|ltb4{CUSTOM_SP}receptor)\b"
        ),
        ["ltb4r", "ltb4r2"],
        False,
    ),
    # ===== P2Y / P2Y-like =====
    (re.compile(rf"\bp2y1{CUSTOM_SP}receptor\b"), ["p2ry1"], False),
    (re.compile(rf"\bp2y2{CUSTOM_SP}receptor\b"), ["p2ry2"], False),
    (re.compile(rf"\bp2y4{CUSTOM_SP}receptor\b"), ["p2ry4"], False),
    (re.compile(rf"\bp2y6{CUSTOM_SP}receptor\b"), ["p2ry6"], False),
    (re.compile(rf"\bp2y12{CUSTOM_SP}receptor\b"), ["p2ry12"], False),
    (re.compile(rf"\bp2y14{CUSTOM_SP}receptor\b"), ["p2ry14"], False),
    # ===== S1P / LPA =====
    (
        re.compile(
            rf"\b(s1p5{CUSTOM_SP}receptor|sphingosine{CUSTOM_SP}1{CUSTOM_SP}phosphate{CUSTOM_SP}receptor{CUSTOM_SP}5)\b"
        ),
        ["s1pr5"],
        False,
    ),
    (
        re.compile(
            rf"\b(s1p3{CUSTOM_SP}receptor|sphingosine{CUSTOM_SP}1{CUSTOM_SP}phosphate{CUSTOM_SP}receptor{CUSTOM_SP}3)\b"
        ),
        ["s1pr3"],
        False,
    ),
    (
        re.compile(
            rf"\b(s1p1r{CUSTOM_SP}receptor|s1p1{CUSTOM_SP}receptor|sphingosine{CUSTOM_SP}1{CUSTOM_SP}phosphate{CUSTOM_SP}receptor{CUSTOM_SP}1)\b"
        ),
        ["s1pr1"],
        False,
    ),
    (re.compile(rf"\bs1p{CUSTOM_SP}receptor\b"), ["s1pr1"], True),
    (
        re.compile(
            rf"\b(s1p4{CUSTOM_SP}receptor|sphingosine{CUSTOM_SP}1{CUSTOM_SP}phosphate{CUSTOM_SP}receptor{CUSTOM_SP}4)\b"
        ),
        ["s1pr4"],
        False,
    ),
    (
        re.compile(
            rf"\b(s1p2{CUSTOM_SP}receptor|sphingosine{CUSTOM_SP}1{CUSTOM_SP}phosphate{CUSTOM_SP}receptor{CUSTOM_SP}2)\b"
        ),
        ["s1pr2"],
        False,
    ),
    (re.compile(rf"\b(lpa1{CUSTOM_SP}receptor)\b"), ["lpar1"], False),
    (re.compile(rf"\b(lpa5{CUSTOM_SP}receptor)\b"), ["lpar5"], False),
    # ===== LIGAND GPCRs =====
    (re.compile(rf"\bgpr35{CUSTOM_SP}receptor\b"), ["gpr35"], False),
    (
        re.compile(rf"\bgpr103{CUSTOM_SP}receptor\b"),
        ["gpr103", "qrfpr"],
        False,
    ),
    (
        re.compile(rf"\bgpr10a{CUSTOM_SP}receptor\b"),
        ["prlhr", "gpr10"],
        False,
    ),
    (
        re.compile(rf"\bgpr154{CUSTOM_SP}receptor\b"),
        ["npsr1", "gpr154"],
        False,
    ),
    (re.compile(rf"\bebi2{CUSTOM_SP}receptor\b"), ["gpr183", "ebi2"], False),
    (
        re.compile(rf"\b(ffar1{CUSTOM_SP}receptor|ffa1{CUSTOM_SP}receptor)\b"),
        ["ffar1", "gpr40"],
        False,
    ),
    (re.compile(rf"\bffa2{CUSTOM_SP}receptor\b"), ["ffar2", "gpr43"], False),
    # ===== OXE / HCA (niacin) =====
    (
        re.compile(
            rf"\b(oxe{CUSTOM_SP}receptor|5{CUSTOM_SP}oxo{CUSTOM_SP}ete{CUSTOM_SP}receptor)\b"
        ),
        ["oxer1"],
        False,
    ),
    (
        re.compile(
            rf"\b(grp109a{CUSTOM_SP}receptor|hca2{CUSTOM_SP}receptor|niacin{CUSTOM_SP}receptor)\b"
        ),
        ["hcar2", "gpr109a"],
        False,
    ),
    # ===== OPIOID FAMILY =====
    (
        re.compile(
            rf"\b(mu|mop|mu-?type|mu1|mu2){CUSTOM_SP}(opioid|receptor)\b|\b(opioid{CUSTOM_SP}receptor{CUSTOM_SP}mu(\s*1|\s*2)?)\b|\bmu{CUSTOM_SP}opioid{CUSTOM_SP}receptor\b"
        ),
        ["mor", "oprm1"],
        False,
    ),
    (
        re.compile(
            rf"\b(kappa|kop|k-?opioid|kappa1){CUSTOM_SP}(opioid|receptor)\b|\bopioid{CUSTOM_SP}receptor{CUSTOM_SP}kappa\b"
        ),
        ["kor", "oprk1"],
        False,
    ),
    (
        re.compile(
            rf"\b(delta){CUSTOM_SP}(opioid|receptor)\b|\bopioid{CUSTOM_SP}receptor{CUSTOM_SP}delta\b"
        ),
        ["dor", "oprd1"],
        False,
    ),
    (
        re.compile(
            rf"\b(nociceptin{CUSTOM_SP}opioid{CUSTOM_SP}receptor|nociceptin{CUSTOM_SP}receptor|orphanin{CUSTOM_SP}fq{CUSTOM_SP}receptor|horl1{CUSTOM_SP}receptor|orl1{CUSTOM_SP}receptor)\b"
        ),
        ["nop", "orl1", "oprl1"],
        False,
    ),
    (
        re.compile(rf"\bopioid{CUSTOM_SP}receptor(s)?\b"),
        ["oprm1", "oprk1", "oprd1", "oprl1"],
        True,
    ),
    # ===== MELANIN CONCENTRATING HORMONE (MCHR1) =====
    (
        re.compile(
            rf"\b(machr1{CUSTOM_SP}receptor|mch1\){CUSTOM_SP}receptor|mch1r{CUSTOM_SP}receptor|mch1{CUSTOM_SP}receptor|melanin{CUSTOM_SP}concentrating{CUSTOM_SP}hormone(\b|{CUSTOM_SP})1{CUSTOM_SP}receptor|melanin{CUSTOM_SP}-{CUSTOM_SP}concentrating{CUSTOM_SP}hormone(\b|{CUSTOM_SP})1{CUSTOM_SP}receptor)\b"
        ),
        ["mchr1"],
        False,
    ),
    # ===== GHSR =====
    (
        re.compile(
            rf"\b(ghsr1a{CUSTOM_SP}receptor|growth{CUSTOM_SP}hormone{CUSTOM_SP}secretagogue{CUSTOM_SP}receptor|ghs1a{CUSTOM_SP}receptor|ghs{CUSTOM_SP}receptor|grln{CUSTOM_SP}receptor)\b"
        ),
        ["ghsr", "ghsr1a"],
        False,
    ),
    # ===== ESTROGEN RECEPTORS =====
    (
        re.compile(
            rf"\b(erbeta{CUSTOM_SP}receptor|estrogen{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}beta|estrogen{CUSTOM_SP}receptor{CUSTOM_SP}2)\b"
        ),
        ["er beta", "esr2"],
        False,
    ),
    (
        re.compile(
            rf"\b(estrogen{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}alpha|eralpha{CUSTOM_SP}receptor)\b"
        ),
        ["er alpha", "esr1"],
        False,
    ),
    # ===== MUSCARINIC =====
    (
        re.compile(
            rf"\bmuscarinic{CUSTOM_SP}(acetylcholine{CUSTOM_SP})?receptor{CUSTOM_SP}m3\b|\bmuscarinic{CUSTOM_SP}m3{CUSTOM_SP}receptor\b|\bm3{CUSTOM_SP}muscarinic{CUSTOM_SP}receptor\b"
        ),
        ["m3", "chrm3"],
        False,
    ),
    (
        re.compile(
            rf"\bmuscarinic{CUSTOM_SP}(acetylcholine{CUSTOM_SP})?receptor{CUSTOM_SP}m2\b|\bmuscarinic{CUSTOM_SP}m2{CUSTOM_SP}receptor\b"
        ),
        ["m2", "chrm2"],
        False,
    ),
    (
        re.compile(
            rf"\bmuscarinic{CUSTOM_SP}(acetylcholine{CUSTOM_SP})?receptor{CUSTOM_SP}m1\b|\bmuscarinic{CUSTOM_SP}m1{CUSTOM_SP}receptor\b"
        ),
        ["m1", "chrm1"],
        False,
    ),
    (
        re.compile(
            rf"\bmuscarinic{CUSTOM_SP}(acetylcholine{CUSTOM_SP})?receptor{CUSTOM_SP}m4\b|\bmuscarinic{CUSTOM_SP}m4{CUSTOM_SP}receptor\b"
        ),
        ["m4", "chrm4"],
        False,
    ),
    (
        re.compile(
            rf"\bmuscarinic{CUSTOM_SP}(acetylcholine{CUSTOM_SP})?receptor{CUSTOM_SP}m5\b|\bmuscarinic{CUSTOM_SP}m5{CUSTOM_SP}receptor\b"
        ),
        ["m5", "chrm5"],
        False,
    ),
    # ===== ADENOSINE =====
    (
        re.compile(
            rf"\b(a2b{CUSTOM_SP}receptor|adenosine{CUSTOM_SP}2b{CUSTOM_SP}receptor)\b"
        ),
        ["a2b", "adora2b"],
        False,
    ),
    (
        re.compile(
            rf"\b(a2a{CUSTOM_SP}receptor|adenosine{CUSTOM_SP}2a{CUSTOM_SP}receptor)\b"
        ),
        ["a2a", "adora2a"],
        False,
    ),
    (
        re.compile(
            rf"\b(adenosine{CUSTOM_SP}a3a{CUSTOM_SP}receptor|a3a{CUSTOM_SP}receptor|a3{CUSTOM_SP}receptor)\b"
        ),
        ["a3", "adora3"],
        False,
    ),
    (
        re.compile(
            rf"\b(adenosine{CUSTOM_SP}a{CUSTOM_SP}receptor|adenosine{CUSTOM_SP}receptor)\b"
        ),
        ["adora1", "adora2a", "adora2b", "adora3"],
        True,
    ),
    # ===== MELATONIN =====
    (
        re.compile(
            rf"\b(mt2(\s|$)|mt2{CUSTOM_SP}melatonin{CUSTOM_SP}receptor|melatonin{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}2|melatonin{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}2|melatonin{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}1b|melatonin{CUSTOM_SP}mt2{CUSTOM_SP}receptor)\b"
        ),
        ["mt2", "mtnr1b"],
        False,
    ),
    (
        re.compile(
            rf"\b(mt1(\s|$)|mt1a{CUSTOM_SP}receptor|ml1a{CUSTOM_SP}receptor|mt1{CUSTOM_SP}melatonin{CUSTOM_SP}receptor|melatonin{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}1|melatonin{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}1a)\b"
        ),
        ["mt1", "mtnr1a"],
        False,
    ),
    # ===== TRP =====
    (re.compile(rf"\btrmp8{CUSTOM_SP}receptor\b"), ["trpm8"], False),
    (
        re.compile(
            rf"\b(transient{CUSTOM_SP}receptor{CUSTOM_SP}potential{CUSTOM_SP}vanilloid(\s*type)?{CUSTOM_SP}1|vanilloid{CUSTOM_SP}receptor(\s*subtype)?{CUSTOM_SP}vr?1|vr1{CUSTOM_SP}receptor|trpv1{CUSTOM_SP}receptor)\b"
        ),
        ["trpv1", "vr1"],
        False,
    ),
    # ===== SMO =====
    (
        re.compile(rf"\b(smo{CUSTOM_SP}receptor|smoothened{CUSTOM_SP}receptor)\b"),
        ["smo"],
        False,
    ),
    # ===== MGLU / iGluR / NMDA / AMPA / KA =====
    (
        re.compile(rf"\bmglu1(\b|{CUSTOM_SP}receptor|{CUSTOM_SP}receptors)\b"),
        ["grm1"],
        False,
    ),
    (
        re.compile(rf"\bmglu2(\b|{CUSTOM_SP}receptor|{CUSTOM_SP}receptors)\b"),
        ["grm2"],
        False,
    ),
    (
        re.compile(rf"\bmglu3(\b|{CUSTOM_SP}receptor|{CUSTOM_SP}receptors)\b"),
        ["grm3"],
        False,
    ),
    (
        re.compile(
            rf"\bmglu5a{CUSTOM_SP}receptor\b|\bmglu5{CUSTOM_SP}receptor\b|\bglutamate{CUSTOM_SP}receptor{CUSTOM_SP}5\b"
        ),
        ["grm5"],
        False,
    ),
    (re.compile(rf"\bmglu7{CUSTOM_SP}receptor\b"), ["grm7"], False),
    (re.compile(rf"\biglur5{CUSTOM_SP}receptor\b"), ["grik1"], False),
    (
        re.compile(
            rf"\biglur6{CUSTOM_SP}receptor|iontropic{CUSTOM_SP}glutamate{CUSTOM_SP}receptor{CUSTOM_SP}6\b"
        ),
        ["grik2"],
        False,
    ),
    (re.compile(rf"\biglur7{CUSTOM_SP}receptor\b"), ["grik3"], False),
    (re.compile(rf"\bglun2a{CUSTOM_SP}receptor\b"), ["grin2a"], False),
    (re.compile(rf"\bglun2b{CUSTOM_SP}receptor\b"), ["grin2b"], False),
    (re.compile(rf"\beaat1{CUSTOM_SP}receptor\b"), ["slc1a3"], False),
    # ===== SEROTONIN =====
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine|5{CUSTOM_SP}?-?ht){CUSTOM_SP}6{CUSTOM_SP}receptor|ht6{CUSTOM_SP}receptor\b"
        ),
        ["htr6"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}5a{CUSTOM_SP}receptor\b"
        ),
        ["htr5a"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}2a{CUSTOM_SP}receptor\b"
        ),
        ["htr2a"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}2b{CUSTOM_SP}receptor\b"
        ),
        ["htr2b"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}1a{CUSTOM_SP}receptor\b|hydroxytryptamine{CUSTOM_SP}1a{CUSTOM_SP}receptor\b"
        ),
        ["htr1a"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}1b{CUSTOM_SP}receptor\b"
        ),
        ["htr1b"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}1d{CUSTOM_SP}receptor\b"
        ),
        ["htr1d"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}1f{CUSTOM_SP}receptor\b"
        ),
        ["htr1f"],
        False,
    ),
    (
        re.compile(
            rf"\b(5{CUSTOM_SP}?-?ht|5{CUSTOM_SP}hydroxy{CUSTOM_SP}?tryptamine){CUSTOM_SP}7{CUSTOM_SP}receptor\b"
        ),
        ["htr7"],
        False,
    ),
    (re.compile(rf"\b(ht3a{CUSTOM_SP}receptor)\b"), ["htr3a"], False),
    # ===== NPY =====
    (re.compile(rf"\bnpyy5{CUSTOM_SP}receptor\b"), ["npy5r"], False),
    (re.compile(rf"\bnpyy2{CUSTOM_SP}receptor\b"), ["npy2r"], False),
    (re.compile(rf"\bny4{CUSTOM_SP}receptor\b"), ["npy4r"], False),
    # ===== VASOPRESSIN / OXYTOCIN =====
    (
        re.compile(rf"\bvasopressin{CUSTOM_SP}receptor\b"),
        ["avpr1a", "avpr1b", "avpr2"],
        True,
    ),
    (
        re.compile(
            rf"\bvasopressin{CUSTOM_SP}v1a{CUSTOM_SP}receptor|v1a{CUSTOM_SP}receptor\b"
        ),
        ["avpr1a"],
        False,
    ),
    (
        re.compile(
            rf"\bvasopressin{CUSTOM_SP}v1b{CUSTOM_SP}receptor|v1b{CUSTOM_SP}receptor\b"
        ),
        ["avpr1b"],
        False,
    ),
    (
        re.compile(
            rf"\bvasopressin{CUSTOM_SP}v2{CUSTOM_SP}receptor|v2{CUSTOM_SP}receptor|v2{CUSTOM_SP}vasopressin{CUSTOM_SP}receptor\b"
        ),
        ["avpr2"],
        False,
    ),
    (
        re.compile(
            rf"\boxytocin{CUSTOM_SP}receptor|ot(r|{CUSTOM_SP}receptor)|oxr{CUSTOM_SP}receptor\b"
        ),
        ["oxtr"],
        False,
    ),
    # ===== ANGIOTENSIN =====
    (
        re.compile(
            rf"\b(angiotensin{CUSTOM_SP}ii{CUSTOM_SP}receptor(\s*,)?(\s*type)?{CUSTOM_SP}1|at-?1{CUSTOM_SP}receptor|angiotensin{CUSTOM_SP}2{CUSTOM_SP}type{CUSTOM_SP}-?{CUSTOM_SP}1{CUSTOM_SP}receptor|angiotensin{CUSTOM_SP}2{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}1|angiotensin{CUSTOM_SP}2{CUSTOM_SP}at1{CUSTOM_SP}receptor|angiotensin{CUSTOM_SP}1{CUSTOM_SP}receptor)\b"
        ),
        ["agtr1"],
        False,
    ),
    (
        re.compile(
            rf"\b(type{CUSTOM_SP}-?{CUSTOM_SP}2{CUSTOM_SP}angiotensin{CUSTOM_SP}-?{CUSTOM_SP}2{CUSTOM_SP}receptor|at2{CUSTOM_SP}receptor|angiotensin{CUSTOM_SP}ii{CUSTOM_SP}receptor(\s*,)?(\s*type)?{CUSTOM_SP}2)\b"
        ),
        ["agtr2"],
        False,
    ),
    # ===== CHEMOKINE =====
    (
        re.compile(
            rf"\bchemokine{CUSTOM_SP}receptor{CUSTOM_SP}5\b|\bc-?c{CUSTOM_SP}chemokine{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}5\b"
        ),
        ["ccr5"],
        False,
    ),
    (
        re.compile(
            rf"\bc-?c{CUSTOM_SP}chemokine{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}4\b"
        ),
        ["ccr4"],
        False,
    ),
    (
        re.compile(
            rf"\bc-?c{CUSTOM_SP}chemokine{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}3\b"
        ),
        ["ccr3"],
        False,
    ),
    (
        re.compile(
            rf"\bcxc{CUSTOM_SP}chemokine{CUSTOM_SP}receptor{CUSTOM_SP}1\b|interleukin{CUSTOM_SP}-?{CUSTOM_SP}8{CUSTOM_SP}receptor\b"
        ),
        ["cxcr1"],
        False,
    ),
    (
        re.compile(rf"\bcxc{CUSTOM_SP}chemokine{CUSTOM_SP}receptor{CUSTOM_SP}2\b"),
        ["cxcr2"],
        False,
    ),
    (
        re.compile(rf"\bcx3c{CUSTOM_SP}chemokine{CUSTOM_SP}receptor{CUSTOM_SP}[35]\b"),
        ["cx3cr1"],
        True,
    ),
    # ===== MELANOCORTIN =====
    (re.compile(r"\bmc3r(eceptor)?\b"), ["mc3r"], False),
    # ===== CNR / CANNABINOID =====
    (
        re.compile(
            rf"\b(cb1{CUSTOM_SP}receptor|cannabinoid{CUSTOM_SP}-?{CUSTOM_SP}1{CUSTOM_SP}receptor|cannabinoid{CUSTOM_SP}receptor{CUSTOM_SP}1|cannabinoid{CUSTOM_SP}cb1{CUSTOM_SP}receptor)\b"
        ),
        ["cb1", "cnr1"],
        False,
    ),
    (
        re.compile(
            rf"\b(cb2{CUSTOM_SP}receptor|cannabinoid{CUSTOM_SP}receptors?{CUSTOM_SP}2|cannabinoid{CUSTOM_SP}cb2{CUSTOM_SP}receptor|cannabinoid{CUSTOM_SP}receptor{CUSTOM_SP}2)\b"
        ),
        ["cb2", "cnr2"],
        False,
    ),
    # ===== DP/EP/IP/TP (PROSTANOID) =====
    (
        re.compile(
            rf"\bprostaglandin{CUSTOM_SP}i2{CUSTOM_SP}receptor|pgi2{CUSTOM_SP}receptor\b"
        ),
        ["ip", "ptgir"],
        False,
    ),
    (
        re.compile(rf"\bprostaglandin{CUSTOM_SP}e{CUSTOM_SP}receptor\b"),
        ["ptger1", "ptger2", "ptger3", "ptger4"],
        True,
    ),
    (
        re.compile(
            rf"\bep3(alpha|c)?{CUSTOM_SP}receptor|prostaglandin{CUSTOM_SP}ep3alpha{CUSTOM_SP}receptor\b"
        ),
        ["ep3", "ptger3"],
        False,
    ),
    # ===== ADRENERGIC =====
    (
        re.compile(
            rf"\bbeta2(\s*adrenoreceptor|\s*receptor|\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|{CUSTOM_SP}-?receptor)\b"
        ),
        ["adrb2"],
        False,
    ),
    (
        re.compile(
            rf"\bbeta1{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|beta-?{CUSTOM_SP}1{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|beta1{CUSTOM_SP}receptor\b"
        ),
        ["adrb1"],
        False,
    ),
    (
        re.compile(
            rf"\balpha1a(\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|\s*receptor|{CUSTOM_SP}-?adrenoreceptor)\b"
        ),
        ["adra1a"],
        False,
    ),
    (
        re.compile(
            rf"\balpha1b(\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|\s*receptor|{CUSTOM_SP}-?adrenoreceptor)\b"
        ),
        ["adra1b"],
        False,
    ),
    (
        re.compile(
            rf"\balpha1d(\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|\s*receptor|{CUSTOM_SP}-?adrenoreceptor)\b|\balpha1c{CUSTOM_SP}receptor\b"
        ),
        ["adra1d", "adra1a"],
        False,
    ),
    (
        re.compile(
            rf"\balpha2a(\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|\s*receptor)\b"
        ),
        ["adra2a"],
        False,
    ),
    (
        re.compile(
            rf"\balpha2b(\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|\s*receptor)\b"
        ),
        ["adra2b"],
        False,
    ),
    (
        re.compile(
            rf"\balpha2c(\s*-?{CUSTOM_SP}adrenergic{CUSTOM_SP}receptor|\s*receptor)\b"
        ),
        ["adra2c"],
        False,
    ),
    # ===== DOPAMINE =====
    (
        re.compile(
            rf"\bdopaminergic{CUSTOM_SP}d2{CUSTOM_SP}receptor|drd2(\s|\b)|dopamine{CUSTOM_SP}receptor{CUSTOM_SP}d2(l|s)?|dopamine{CUSTOM_SP}d2(l|s)?{CUSTOM_SP}receptor|d2(l|s)?{CUSTOM_SP}receptor\b"
        ),
        ["drd2"],
        False,
    ),
    (
        re.compile(
            rf"\bdopamine{CUSTOM_SP}d3{CUSTOM_SP}receptor|drd3{CUSTOM_SP}receptor|d3{CUSTOM_SP}receptor\b"
        ),
        ["drd3"],
        False,
    ),
    (
        re.compile(
            rf"\bd5{CUSTOM_SP}dopamine{CUSTOM_SP}receptor|d5{CUSTOM_SP}receptor\b"
        ),
        ["drd5"],
        False,
    ),
    (
        re.compile(
            rf"\bd4(\.4)?{CUSTOM_SP}receptor|dopamine{CUSTOM_SP}receptor{CUSTOM_SP}4(\.4)?|d4{CUSTOM_SP}dopamine{CUSTOM_SP}receptor\b"
        ),
        ["drd4"],
        False,
    ),
    (
        re.compile(
            rf"\bd1{CUSTOM_SP}dopamine{CUSTOM_SP}receptor|d1{CUSTOM_SP}receptor\b"
        ),
        ["drd1"],
        False,
    ),
    (
        re.compile(rf"\bdopamine{CUSTOM_SP}receptor(s)?\b"),
        ["drd1", "drd2", "drd3", "drd4", "drd5"],
        True,
    ),
    # ===== NUCLEAR RECEPTORS / прочее НЕ-GPCR =====
    (
        re.compile(
            rf"\bliver{CUSTOM_SP}x{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}alpha|lxr{CUSTOM_SP}alpha{CUSTOM_SP}receptor\b"
        ),
        ["nr1h3"],
        False,
    ),
    (
        re.compile(
            rf"\bliver{CUSTOM_SP}x{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}beta|lxr{CUSTOM_SP}beta{CUSTOM_SP}receptor\b"
        ),
        ["nr1h2"],
        False,
    ),
    (
        re.compile(
            rf"\bppar(alpha|{CUSTOM_SP}alpha){CUSTOM_SP}receptor|peroxisome{CUSTOM_SP}proliferator{CUSTOM_SP}activated{CUSTOM_SP}receptor{CUSTOM_SP}alpha\b"
        ),
        ["ppara"],
        False,
    ),
    (
        re.compile(
            rf"\bppar(gamma|{CUSTOM_SP}-?{CUSTOM_SP}gamma){CUSTOM_SP}receptor|peroxisome{CUSTOM_SP}proliferator{CUSTOM_SP}-?{CUSTOM_SP}activated{CUSTOM_SP}receptor{CUSTOM_SP}gamma\b"
        ),
        ["pparg"],
        False,
    ),
    (
        re.compile(
            rf"\bppar{CUSTOM_SP}delta{CUSTOM_SP}receptor|ppardelta{CUSTOM_SP}receptor\b"
        ),
        ["ppard"],
        False,
    ),
    (
        re.compile(
            rf"\bretinoid{CUSTOM_SP}x{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}alpha|rxr{CUSTOM_SP}alpha\b"
        ),
        ["rxra"],
        False,
    ),
    (
        re.compile(
            rf"\bretinoid{CUSTOM_SP}x{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}beta|rxr{CUSTOM_SP}beta\b"
        ),
        ["rxrb"],
        False,
    ),
    (
        re.compile(
            rf"\bretinoid{CUSTOM_SP}x{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}gamma|rxr{CUSTOM_SP}gamma\b"
        ),
        ["rxrg"],
        False,
    ),
    (
        re.compile(
            rf"\bretinoic{CUSTOM_SP}acid{CUSTOM_SP}receptor(\b|{CUSTOM_SP})-?{CUSTOM_SP}gamma\b"
        ),
        ["rarg"],
        False,
    ),
    (
        re.compile(
            rf"\bestrogen{CUSTOM_SP}related{CUSTOM_SP}receptor{CUSTOM_SP}gamma\b"
        ),
        ["esrrg"],
        False,
    ),
    (
        re.compile(
            rf"\bthyroid{CUSTOM_SP}hormone{CUSTOM_SP}receptor{CUSTOM_SP}beta(\s*-?\s*1)?\b"
        ),
        ["thrb"],
        False,
    ),
    (
        re.compile(
            rf"\bthyroid{CUSTOM_SP}hormone{CUSTOM_SP}receptor{CUSTOM_SP}alpha(\s*-?\s*1)?\b"
        ),
        ["thra"],
        False,
    ),
    (re.compile(rf"\bglucagon{CUSTOM_SP}receptor\b"), ["gcgr"], False),
    (re.compile(rf"\bglp1{CUSTOM_SP}receptor\b"), ["glp1r"], False),
    (re.compile(rf"\bgipr{CUSTOM_SP}receptor\b"), ["gipr"], False),
    (re.compile(rf"\bcd69{CUSTOM_SP}receptor\b"), ["cd69"], False),
    (
        re.compile(
            rf"\b(alpha|nicotinic){CUSTOM_SP}(\-?\s*)7{CUSTOM_SP}(nicotinic{CUSTOM_SP})?acetylcholine{CUSTOM_SP}receptor|chrna7\b"
        ),
        ["chrna7", "nachr alpha7"],
        False,
    ),
    (re.compile(rf"\bc3a{CUSTOM_SP}receptor\b"), ["c3ar1"], False),
    (re.compile(rf"\bgp6{CUSTOM_SP}receptor\b"), ["gp6"], False),
    (
        re.compile(
            rf"\bppar(gamma|{CUSTOM_SP}-?{CUSTOM_SP}gamma){CUSTOM_SP}receptor\b"
        ),
        ["pparg"],
        False,
    ),
    (re.compile(rf"\bcalcium{CUSTOM_SP}sensing{CUSTOM_SP}receptor\b"), ["casr"], False),
    (re.compile(rf"\bip3{CUSTOM_SP}receptor\b"), ["itpr1", "itpr2", "itpr3"], True),
    (
        re.compile(
            rf"\bthromboxane{CUSTOM_SP}a2{CUSTOM_SP}receptor(\s*alpha|\s*beta)?|tp(alpha|beta)?{CUSTOM_SP}receptor\b"
        ),
        ["tbxa2r", "tpalpha", "tpbeta"],
        False,
    ),
    (re.compile(rf"\bc5a{CUSTOM_SP}receptor\b"), ["c5ar1", "c5ar2"], True),
    (re.compile(rf"\bebi2{CUSTOM_SP}receptor\b"), ["gpr183", "ebi2"], False),
    (re.compile(rf"\bapj{CUSTOM_SP}receptor\b"), ["aplnr"], False),
    (
        re.compile(
            rf"\bpar2{CUSTOM_SP}receptor|protease{CUSTOM_SP}activated{CUSTOM_SP}receptor{CUSTOM_SP}1\b"
        ),
        ["f2rl1", "f2r"],
        False,
    ),
    (
        re.compile(
            rf"\bpar1{CUSTOM_SP}receptor|protease{CUSTOM_SP}activated{CUSTOM_SP}receptor{CUSTOM_SP}1\b"
        ),
        ["f2r"],
        False,
    ),
    (
        re.compile(
            rf"\bendothelin{CUSTOM_SP}receptor{CUSTOM_SP}type{CUSTOM_SP}a|et-?a{CUSTOM_SP}receptor|eta{CUSTOM_SP}receptor\b"
        ),
        ["ednra"],
        False,
    ),
    (
        re.compile(
            rf"\bendothelin{CUSTOM_SP}b{CUSTOM_SP}receptor|et-?b{CUSTOM_SP}receptor|ednrb\b"
        ),
        ["ednrb"],
        False,
    ),
    (
        re.compile(rf"\bepha2(\s*-?{CUSTOM_SP}fc)?{CUSTOM_SP}receptor\b"),
        ["epha2"],
        False,
    ),
    (
        re.compile(r"\bvegf(\s*receptor\s*2|\s*receptor\s*2\s*\(kdr\))|kdr\b"),
        ["kdr"],
        False,
    ),
    (
        re.compile(
            rf"\begf{CUSTOM_SP}receptor|epidermal{CUSTOM_SP}growth{CUSTOM_SP}factor{CUSTOM_SP}receptor\b"
        ),
        ["egfr"],
        False,
    ),
    (
        re.compile(
            rf"\bpdgf{CUSTOM_SP}receptor\b|platelet{CUSTOM_SP}derived{CUSTOM_SP}growth{CUSTOM_SP}factor{CUSTOM_SP}receptor{CUSTOM_SP}beta\b"
        ),
        ["pdgfrb"],
        False,
    ),
    (re.compile(rf"\binsulin{CUSTOM_SP}receptor\b"), ["insr"], False),
    (re.compile(rf"\bigf-?1{CUSTOM_SP}receptor\b"), ["igf1r"], False),
    (re.compile(rf"\bdat{CUSTOM_SP}receptor\b"), ["slc6a3"], False),
    (
        re.compile(rf"\bcocaine{CUSTOM_SP}receptor\b"),
        ["slc6a3", "sigmar1"],
        True,
    ),
    (re.compile(rf"\bnet{CUSTOM_SP}receptor\b"), ["slc6a2"], False),
    (
        re.compile(rf"\b5htt{CUSTOM_SP}receptor|sert{CUSTOM_SP}receptor\b"),
        ["slc6a4"],
        False,
    ),
    (re.compile(rf"\bgat1{CUSTOM_SP}receptor\b"), ["slc6a1"], False),
    (
        re.compile(rf"\bherg{CUSTOM_SP}receptor|erg{CUSTOM_SP}receptor\b"),
        ["kcnh2"],
        False,
    ),
    (
        re.compile(rf"\bgabac{CUSTOM_SP}rho(\s*-?\s*)1{CUSTOM_SP}receptor\b"),
        ["gabrr1"],
        False,
    ),
    (
        re.compile(
            rf"\bgabaalpha2{CUSTOM_SP}receptor\b|alpha{CUSTOM_SP}5{CUSTOM_SP}gaba{CUSTOM_SP}-?{CUSTOM_SP}a{CUSTOM_SP}receptor\b"
        ),
        ["gabra2", "gabra5"],
        True,
    ),
    (
        re.compile(
            rf"\bglycine{CUSTOM_SP}receptor{CUSTOM_SP}subunit{CUSTOM_SP}alpha{CUSTOM_SP}-?\s*1|glycine{CUSTOM_SP}receptor{CUSTOM_SP}alpha{CUSTOM_SP}1\b"
        ),
        ["glra1"],
        False,
    ),
]

# ---------------------------------------------------------------------------
# Ion channel-specific candidate rules

ION_DOT = r"\s*\.\s*"  # require a literal dot between digits
ION_SP = r"\s*(?:-|\s|/)\s*"

RULES_ION_CHANNELS: Sequence[CandidateRule] = [
    # ===== Nav (SCN) sodium channels =====
    (
        re.compile(
            rf"\b(?:na\s*v|nav)\s*1{ION_DOT}7\b|\bnav{ION_SP}1{ION_DOT}7(?:{ION_SP}(?:sodium|ion){ION_SP}channel)?\b"
        ),
        ["scn9a", "nav1.7"],
        False,
    ),
    (
        re.compile(
            rf"\b(?:na\s*v|nav)\s*1{ION_DOT}8\b|\bvoltage{ION_SP}gated{ION_SP}na{ION_SP}channel{ION_SP}1{ION_DOT}8\b"
        ),
        ["scn10a", "nav1.8"],
        False,
    ),
    (re.compile(rf"\b(?:na\s*v|nav)\s*1{ION_DOT}5\b"), ["scn5a", "nav1.5"], False),
    (
        re.compile(
            rf"\b(?:na\s*v|nav)\s*1{ION_DOT}3\b|\bnav{ION_SP}1{ION_DOT}3{ION_SP}sodium{ION_SP}channel\b"
        ),
        ["scn3a", "nav1.3"],
        False,
    ),
    (
        re.compile(
            rf"\b(?:na\s*v|nav)\s*1{ION_DOT}2\b|\bsodium{ION_SP}channel{ION_SP}type{ION_SP}ii?a\b"
        ),
        ["scn2a", "nav1.2"],
        False,
    ),
    # ===== Cav (CACNA) calcium channels =====
    (
        re.compile(rf"\bcav{ION_SP}1{ION_DOT}2{ION_SP}channel\b"),
        ["cacna1c", "cav1.2"],
        False,
    ),
    (
        re.compile(
            rf"\bl{ION_SP}type{ION_SP}(?:\[?\s*ca2\+\s*\]?{ION_SP})?calcium{ION_SP}channel(?:{ION_SP}receptor)?\b"
        ),
        ["cacna1c", "cacna1d"],
        True,
    ),
    (
        re.compile(
            rf"\bcav{ION_SP}2{ION_DOT}2{ION_SP}channel\b|\bn{ION_SP}type{ION_SP}calcium{ION_SP}channels?\b"
        ),
        ["cacna1b", "cav2.2"],
        False,
    ),
    (
        re.compile(
            rf"\bcav{ION_SP}3{ION_DOT}1{ION_SP}channel\b|\balpha{ION_SP}-?1g{ION_SP}t{ION_SP}type{ION_SP}calcium{ION_SP}channel\b|\bcav{ION_SP}3{ION_DOT}1\b"
        ),
        ["cacna1g", "cav3.1"],
        False,
    ),
    (
        re.compile(
            rf"\bcav{ION_SP}3{ION_DOT}2{ION_SP}channel\b|\balpha{ION_SP}-?1h{ION_SP}t{ION_SP}type{ION_SP}calcium{ION_SP}channel\b|\bcav{ION_SP}3{ION_DOT}2\b"
        ),
        ["cacna1h", "cav3.2"],
        False,
    ),
    (
        re.compile(rf"\bcav{ION_SP}3{ION_DOT}3{ION_SP}channel\b"),
        ["cacna1i", "cav3.3"],
        False,
    ),
    (
        re.compile(rf"\bt{ION_SP}type{ION_SP}(?:\[?\s*ca2\+\s*\]?{ION_SP})?channel\b"),
        ["cacna1g", "cacna1h", "cacna1i"],
        True,
    ),
    (
        re.compile(
            rf"\balpha{ION_SP}-?1g{ION_SP}calcium{ION_SP}channel\b|\bt{ION_SP}type{ION_SP}calcium{ION_SP}channel{ION_SP}alpha1g\b"
        ),
        ["cacna1g"],
        False,
    ),
    (re.compile(rf"\bcav{ION_SP}3{ION_DOT}2\b"), ["cacna1h", "cav3.2"], False),
    (re.compile(rf"\bcav{ION_SP}3{ION_DOT}1\b"), ["cacna1g", "cav3.1"], False),
    # ===== Auxiliary Ca2+ subunits =====
    (
        re.compile(
            rf"\bcalcium{ION_SP}channel{ION_SP}alpha2delta\b|\balpha{ION_SP}-?2delta{ION_SP}calcium{ION_SP}channel\b|\balpha{ION_SP}-?2delta{ION_SP}containing{ION_SP}calcium{ION_SP}channel\b"
        ),
        ["cacna2d1", "cacna2d2", "cacna2d3", "cacna2d4"],
        False,
    ),
    # ===== ERG / hERG (KCNH2) & IKs =====
    (
        re.compile(
            rf"\bherg{ION_SP}channel\b|\berg{ION_SP}channel\b|\berg{ION_SP}potassium{ION_SP}channel\b|\bkv11{ION_DOT}1\b|\bikr{ION_SP}channel\b|\berg{ION_SP}k\+\s*channel\b"
        ),
        ["kcnh2", "kv11.1", "herg"],
        False,
    ),
    (
        re.compile(rf"\biks{ION_SP}potassium{ION_SP}channel\b"),
        ["kcnq1", "kcne1"],
        False,
    ),
    # ===== Kv family =====
    (
        re.compile(rf"\bkv1{ION_DOT}3(?:{ION_SP}ion)?{ION_SP}channel\b"),
        ["kcna3", "kv1.3"],
        False,
    ),
    (
        re.compile(rf"\bkv1{ION_DOT}5(?:{ION_SP}ion)?{ION_SP}channel\b"),
        ["kcna5", "kv1.5"],
        False,
    ),
    # ===== Kir / ROMK =====
    (
        re.compile(rf"\bromk1?{ION_SP}channel\b"),
        ["kcnj1", "romk"],
        False,
    ),
    # ===== Anoctamin (ANO1 / TMEM16A) =====
    (re.compile(rf"\bano1{ION_SP}channel\b"), ["tmem16a", "ano1"], False),
    # ===== Exchangers =====
    (
        re.compile(
            rf"\bna\+\s*{ION_SP}h\+\s*exchanger{ION_SP}1\b|\bsodium{ION_SP}\/?hydrogen{ION_SP}exchanger{ION_SP}1\b"
        ),
        ["slc9a1", "nhe1"],
        False,
    ),
    (re.compile(rf"\bsodium{ION_SP}\/?hydrogen{ION_SP}exchanger\b"), ["slc9a1"], True),
    # ===== Gardos channel (KCNN4) =====
    (re.compile(rf"\bgardos{ION_SP}channel\b"), ["kcnn4", "kca3.1"], False),
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
    soft_out: List[str] = []
    rules: Sequence[CandidateRule] = [
        *[(p, v, False) for p, v in list(RULES_GPCR) + list(RULES_GPCR_EXTRA)],
        *RULES_CUSTOM,
        *RULES_ION_CHANNELS,
    ]
    for pat, val, soft in rules:
        m = pat.search(s)
        if m:
            adds = val(m) if callable(val) else val
            if soft:
                soft_out.extend(adds)
            else:
                out.extend(adds)

    # dedupe preserving order with soft matches at the end
    seen: set[str] = set()
    res: List[str] = []
    for x in out + soft_out:
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
    "xi": "11",
    "xii": "12",
    "xiii": "13",
    "xiv": "14",
    "xv": "15",
    "xvi": "16",
    "xvii": "17",
    "xviii": "18",
    "xix": "19",
    "xx": "20",
}

CONTROL_CHARS_RE = re.compile(r"[\u0000-\u001F\u007F]")
MULTI_SPACE_RE = re.compile(r"[\s\t]+")
TYPO_QUOTES_RE = re.compile(r"[“”«»„]|’")
LONG_DASH_RE = re.compile(r"[–—]")
PAREN_RE = re.compile(r"\(([^(]*)\)|\[([^\]]*)\]|\{([^}]*)\)")
# Split on whitespace, punctuation, and dots/commas not between digits
TOKEN_SPLIT_RE = re.compile(r"(?:[\s\-_/:;]|^,|,(?!\d)|(?<=\D),|(?<!\d)\.(?!\d))+")
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
    re.compile(r"\b[A-Z][0-9]+(?:del|ins|dup)[A-Z]*\b", re.IGNORECASE),
    re.compile(r"\b[A-Z][0-9]+fs\*?\d*\b", re.IGNORECASE),
    re.compile(r"(?:Δ|delta)\s?[A-Z][0-9]+", re.IGNORECASE),
    re.compile(r"\b(mutant|variant|mut\.)\b", re.IGNORECASE),
)
#: Tokens resembling mutations but representing legitimate receptor or gene names.
#: Extend this set via the ``mutation_whitelist`` parameter of ``normalize_target_name``.
MUTATION_WHITELIST: set[str] = {
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
        translation = str.maketrans(
            {**(greek or GREEK_LETTERS), **(superscripts or SUPERSCRIPTS)}
        )
    return text.translate(translation)


def replace_roman_numerals(text: str) -> str:
    """Replace standalone Roman numerals with digits.

    The mapping covers numerals from ``II`` to ``XX`` while excluding
    single-letter numerals such as ``v`` or ``x`` to avoid altering valid
    gene symbols. The input should already be lower-cased.
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


def find_mutations(text: str, whitelist: Sequence[str] | None = None) -> List[str]:
    """Extract mutation-like substrings from text.

    Parameters
    ----------
    text:
        Input string in its original form.
    whitelist:
        Additional mutation-like tokens to ignore during extraction.

    Returns
    -------
    List[str]
        Unique mutation substrings in order of appearance.

    Notes
    -----
    If the optional :mod:`hgvs` package is available, its parser is used to
    recognize HGVS-formatted variants beyond the built-in regex patterns.
    """

    allowed = set(MUTATION_WHITELIST)
    if whitelist:
        allowed.update(w.lower() for w in whitelist)

    found: List[str] = []
    for pattern in MUTATION_PATTERNS:
        for match in pattern.finditer(text):
            # Skip letter-digit-letter where the letters are identical (e.g., A123A)
            if pattern is LETTER_DIGIT_LETTER_RE:
                if match.group(1).upper() == match.group(3).upper():
                    continue
            token = match.group(0)
            lower = token.lower()
            if lower in allowed:
                continue
            if any(lower in f.lower() for f in found):
                continue
            found = [f for f in found if f.lower() not in lower]
            if token not in found:
                found.append(token)
    if _HGVS_PARSER is not None:
        tokens = re.findall(r"\S+", text)
        for tok in tokens:
            lower = tok.lower()
            if lower in allowed or any(lower == f.lower() for f in found):
                continue
            try:
                _HGVS_PARSER.parse(tok)
            except HGVSParseError:
                continue
            found.append(tok)
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


def normalize_target_name(
    name: str,
    *,
    strip_mutations: bool = True,
    mutation_whitelist: Sequence[str] | None = None,
    detect_mutations: bool = True,
    taxon: int = 9606,
) -> NormalizationResult:
    """Normalize a single target name.

    Parameters
    ----------
    name:
        Raw target name.
    strip_mutations:
        Remove mutation-like tokens from the results. Set to ``False`` to
        retain mutations in the output tokens.
    mutation_whitelist:
        Additional mutation-like tokens to preserve. Tokens are compared in
        lowercase and extend :data:`MUTATION_WHITELIST`.
    detect_mutations:
        Skip mutation detection when ``False`` for faster processing.
    taxon:
        NCBI taxonomy identifier to store in :class:`NormalizationResult`.

    Returns
    -------
    NormalizationResult
        Structured normalization information.
    """
    raw = name
    whitelist = set(MUTATION_WHITELIST)
    if mutation_whitelist:
        whitelist.update(t.lower() for t in mutation_whitelist)
    stage = sanitize_text(name)
    mutations: List[str] = []
    if detect_mutations:
        mutations = find_mutations(stage, whitelist=list(whitelist))
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

    mutation_tokens = mutation_token_set(mutations) if detect_mutations else set()
    tokens_no_stop_orig = list(tokens_no_stop)
    tokens_alt_orig = list(tokens_alt)
    tokens_base_no_stop_orig = list(tokens_base_no_stop)
    tokens_base_alt_orig = list(tokens_base_alt)

    if detect_mutations and strip_mutations:
        tokens_no_stop = [
            t for t in tokens_no_stop if t not in mutation_tokens or t in whitelist
        ]
        tokens_alt = [
            t for t in tokens_alt if t not in mutation_tokens or t in whitelist
        ]
        tokens_base_no_stop = [
            t for t in tokens_base_no_stop if t not in mutation_tokens or t in whitelist
        ]
        tokens_base_alt = [
            t for t in tokens_base_alt if t not in mutation_tokens or t in whitelist
        ]

    clean_tokens_check = [
        t for t in tokens_no_stop if not re.fullmatch(r"[a-z]$|\d+$", t)
    ]
    hints_mutations_only = False
    if detect_mutations and not clean_tokens_check:
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
        "mutations": mutations if detect_mutations else [],
    }
    if detect_mutations and hints_mutations_only:
        hints["mutations_only"] = True
    return NormalizationResult(
        raw=raw,
        clean_text=clean_text,
        clean_text_alt=clean_text_alt,
        query_tokens=tokens_no_stop,
        gene_like_candidates=final_cleanup(candidates),
        hint_taxon=taxon,
        hints=hints,
        rules_applied=rules_applied,
    )
