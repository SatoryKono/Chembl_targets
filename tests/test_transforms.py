from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))

from mylib.io_utils import read_target_names
from mylib.transforms import (
    apply_receptor_rules,
    normalize_target_name,
    replace_specials,
    sanitize_text,
)


def test_sanitize_and_replace():
    text = "\u0000β receptor"
    clean = sanitize_text(text)
    clean = replace_specials(clean)
    assert clean == "beta receptor"


def test_receptor_rules():
    text = "dopamine d2 receptor"
    new_text, candidates, rules = apply_receptor_rules(text)
    assert new_text == "dopamine d2"
    assert "drd2" in candidates
    assert rules


def test_full_normalization():
    sample = Path("tests/data/sample.csv")
    df = read_target_names(sample)
    result = normalize_target_name(df.loc[0, "target_name"])
    assert result.clean_text.startswith("beta2 adrenergic")
    assert "adrb2" in result.gene_like_candidates


def test_parenthetical_short_token_retained():
    result = normalize_target_name("histamine receptor (h3)")
    assert "h3" in result.query_tokens
    assert result.clean_text.endswith("h3")
    assert result.hints["parenthetical"] == ["h3"]
    assert "hrh3" in result.gene_like_candidates
