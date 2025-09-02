from pathlib import Path
import sys

import pytest

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
    assert result.clean_text_alt.startswith("beta2 adrenergic")
    assert "adrb2" in result.gene_like_candidates


def test_parenthetical_short_token_retained():
    result = normalize_target_name("histamine receptor (h3)")
    assert "h3" in result.query_tokens
    assert result.clean_text.endswith("h3")
    assert result.hints["parenthetical"] == ["h3"]
    assert "hrh3" in result.gene_like_candidates


def test_clean_text_alt_retains_stopwords():
    result = normalize_target_name("histamine receptor channel")
    assert result.clean_text == "histamine"
    assert result.clean_text_alt == "histamine receptor channel"
    assert "receptor" in result.hints["dropped"]
    assert "channel" in result.hints["dropped"]


def test_hyphen_variants_present():
    result = normalize_target_name("β2-adrenergic receptor")
    assert "beta2-adrenergic" in result.query_tokens
    assert "beta2adrenergic" in result.query_tokens
    assert "beta2-adrenergic" in result.clean_text.split()
    assert "beta2adrenergic" in result.clean_text.split()


def test_letter_digit_space_variants():
    result = normalize_target_name("h 3 receptor")
    assert "h" in result.query_tokens
    assert "3" in result.query_tokens
    assert "h3" in result.query_tokens
    assert "h-3" in result.query_tokens
    assert result.clean_text.split() == ["h3", "h-3"]


def test_parenthetical_complex_indices():
    res1 = normalize_target_name("p2x receptor (p2x7)")
    assert "p2x7" in res1.query_tokens
    assert "p2x7" in res1.clean_text.split()
    assert res1.hints["parenthetical"] == ["p2x7"]

    res2 = normalize_target_name("serotonin receptor (5-ht1a)")
    assert "5-ht1a" in res2.query_tokens
    assert "5ht1a" in res2.query_tokens
    assert "5-ht1a" in res2.clean_text.split()
    assert "5ht1a" in res2.clean_text.split()
    assert res2.hints["parenthetical"] == ["5-ht1a"]


@pytest.mark.parametrize(
    "name,expected",
    [
        ("histamine h4", "hrh4"),
        ("dopamine d3", "drd3"),
        ("adrenergic beta1", "adrb1"),
        ("p2x3", "p2rx3"),
        ("5-ht1b", "htr1b"),
        ("gaba a alpha2", "gabra2"),
    ],
)
def test_regex_gene_like_candidates(name: str, expected: str) -> None:
    result = normalize_target_name(name)
    assert expected in result.gene_like_candidates
