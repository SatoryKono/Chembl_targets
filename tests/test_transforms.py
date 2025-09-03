from pathlib import Path
import sys
import json

import pandas as pd
import pytest
from typing import List, cast

sys.path.append(str(Path(__file__).resolve().parents[1]))

from main import normalize_dataframe
from mylib.io_utils import read_target_names, write_with_new_columns
from mylib.transforms import (
    apply_receptor_rules,
    normalize_target_name,
    replace_specials,
    replace_roman_numerals,
    sanitize_text,
    normalize_unicode,
)


def test_sanitize_and_replace():
    text = "\u0000β receptor"
    clean = sanitize_text(text)
    clean = replace_specials(clean)
    assert clean == "beta receptor"


def test_translate_micro_and_superscript():
    text = normalize_unicode("µp² receptor")
    clean = replace_specials(text)
    assert clean == "mup2 receptor"


def test_replace_roman_numerals():
    text = "type viii receptor"
    replaced = replace_roman_numerals(text)
    assert replaced == "type 8 receptor"


def test_read_target_names_missing_column(tmp_path: Path) -> None:
    csv_path = tmp_path / "x.csv"
    csv_path.write_text("foo\n1\n2\n")
    with pytest.raises(KeyError) as exc:
        read_target_names(csv_path, column="bar")
    msg = str(exc.value)
    assert "Available columns" in msg
    assert "foo" in msg
    assert "1" in msg


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
    assert result.clean_text_alt.startswith("beta2 adrenergic receptor")
    assert "adrb2" in result.gene_like_candidates


def test_parenthetical_short_token_retained():
    result = normalize_target_name("histamine receptor (h3)")
    assert "h3" in result.query_tokens
    assert "histamine h3" in result.clean_text.split("|")
    assert "h3" in result.clean_text.split("|")
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
    assert "beta2 adrenergic" in result.clean_text.split("|")
    assert "beta2-adrenergic" in result.clean_text.split("|")
    assert "beta2adrenergic" in result.clean_text.split("|")


def test_letter_digit_space_variants():
    result = normalize_target_name("h 3 receptor")
    assert "h" in result.query_tokens
    assert "3" in result.query_tokens
    assert "h3" in result.query_tokens
    assert "h-3" in result.query_tokens
    assert result.clean_text.split("|") == ["h3", "h-3"]


def test_parenthetical_complex_indices():
    res1 = normalize_target_name("p2x receptor (p2x7)")
    assert "p2x7" in res1.query_tokens
    assert "p2x p2x7" in res1.clean_text.split("|")
    assert "p2x7" in res1.clean_text.split("|")
    assert res1.hints["parenthetical"] == ["p2x7"]

    res2 = normalize_target_name("serotonin receptor (5-ht1a)")
    assert "5-ht1a" in res2.query_tokens
    assert "5ht1a" in res2.query_tokens
    assert "serotonin 5-ht1a" in res2.clean_text.split("|")
    assert "serotonin 5ht1a" in res2.clean_text.split("|")
    assert "5-ht1a" in res2.clean_text.split("|")
    assert "5ht1a" in res2.clean_text.split("|")
    assert res2.hints["parenthetical"] == ["5-ht1a"]


def test_dataframe_hints_and_rules(tmp_path: Path) -> None:
    sample = Path("tests/data/sample.csv")
    df = read_target_names(sample)
    df_norm = normalize_dataframe(df, "target_name")
    assert isinstance(df_norm.loc[0, "hints"], dict)
    assert isinstance(df_norm.loc[0, "rules_applied"], list)
    out = tmp_path / "out.csv"
    write_with_new_columns(df_norm, out)
    saved = pd.read_csv(out)
    loaded_hints = json.loads(saved.loc[0, "hints"])
    loaded_rules = json.loads(saved.loc[0, "rules_applied"])
    assert isinstance(loaded_hints, dict)
    assert isinstance(loaded_rules, list)


@pytest.mark.parametrize(
    "name,expected",
    [
        ("histamine h4", "hrh4"),
        ("dopamine d3", "drd3"),
        ("adrenergic beta1", "adrb1"),
        ("p2x3", "p2rx3"),
        ("5-ht1b", "htr1b"),
        ("gaba a alpha2", "gabra2"),
        ("trp v 1 channel", "trpv1"),
        ("ampa glua1 receptor", "gria1"),
        ("nmda nr2b receptor", "grin2b"),
        ("kainate gluk3 receptor", "grik3"),
        ("mglur5 receptor", "grm5"),
        ("chemokine cc receptor 5", "ccr5"),
        ("cxcr4 receptor", "cxcr4"),
        ("adenosine a2a receptor", "adora2a"),
        ("nociceptin receptor", "oprl1"),
        ("npy y1 receptor", "npy1r"),
        ("melanocortin 4 receptor", "mc4r"),
        ("prostaglandin ep3 receptor", "ptger3"),
    ],
)
def test_regex_gene_like_candidates(name: str, expected: str) -> None:
    result = normalize_target_name(name)
    assert expected in result.gene_like_candidates


def test_family_level_candidates() -> None:
    res = normalize_target_name("ampa receptor")
    assert {"gria1", "gria2", "gria3", "gria4"} <= set(res.gene_like_candidates)
    res = normalize_target_name("nmda receptor")
    assert {
        "grin1",
        "grin2a",
        "grin2b",
        "grin2c",
        "grin2d",
        "grin3a",
        "grin3b",
    } <= set(res.gene_like_candidates)
    res = normalize_target_name("metabotropic glutamate receptor")
    assert {f"grm{i}" for i in range(1, 9)} <= set(res.gene_like_candidates)


def test_alias_candidates() -> None:
    res = normalize_target_name("sdf-1")
    assert "cxcr4" in res.gene_like_candidates
    res = normalize_target_name("il8")
    assert {"cxcr1", "cxcr2"} <= set(res.gene_like_candidates)
    res = normalize_target_name("rantes")
    assert {"ccr1", "ccr3", "ccr5"} <= set(res.gene_like_candidates)
    res = normalize_target_name("fractalkine")
    assert "cx3cr1" in res.gene_like_candidates


def test_gpcr_family_and_alias_order() -> None:
    res = normalize_target_name("adenosine receptor")
    assert {"adora1", "adora2a", "adora2b", "adora3"} <= set(res.gene_like_candidates)

    res = normalize_target_name("adenosine a2a receptor")
    assert res.gene_like_candidates[:2] == ["a2a", "adora2a"]

    res = normalize_target_name("nociceptin receptor")
    assert res.gene_like_candidates == ["nop", "orl1", "oprl1"]

    res = normalize_target_name("neuropeptide y receptor")
    assert {"npy1r", "npy2r", "npy4r", "npy5r"} <= set(res.gene_like_candidates)

    res = normalize_target_name("npy y1 receptor")
    assert res.gene_like_candidates[:2] == ["y1", "npy1r"]

    res = normalize_target_name("melanocortin 4 receptor")
    assert res.gene_like_candidates == ["mc4r"]

    res = normalize_target_name("prostaglandin ep3 receptor")
    assert res.gene_like_candidates[:2] == ["ep3", "ptger3"]


def test_gpcr_extra_rules() -> None:
    res = normalize_target_name("kisspeptin receptor")
    assert res.gene_like_candidates[:2] == ["kiss1r", "gpr54"]

    res = normalize_target_name("gpr120")
    assert res.gene_like_candidates[:2] == ["ffar4", "gpr120"]

    res = normalize_target_name("alx")
    assert res.gene_like_candidates == ["fpr2"]

    res = normalize_target_name("tgr5")
    assert res.gene_like_candidates[:2] == ["gpbar1", "tgr5"]


def test_mutation_extraction_and_removal() -> None:
    res = normalize_target_name("hiv1 protease I84V mutant")
    assert res.clean_text == "hiv1 protease"
    assert res.hints["mutations"] == ["I84V", "mutant"]

    res = normalize_target_name("BRAF V600E")
    assert res.clean_text == "braf"
    assert res.hints["mutations"] == ["V600E"]

    res = normalize_target_name("p.Gly12Asp KRAS")
    assert res.clean_text == "kras"
    assert res.hints["mutations"] == ["p.Gly12Asp"]

    res = normalize_target_name("CFTR ΔF508")
    assert res.clean_text == "cftr"
    assert res.hints["mutations"] == ["ΔF508"]


def test_letter_digit_letter_same_not_mutation() -> None:
    """Sequences like A123A should not be treated as mutations."""
    res = normalize_target_name("AKT1 E17E")
    assert res.clean_text == "akt1 e17e"
    assert res.hints["mutations"] == []


def test_mutation_whitelist() -> None:
    res = normalize_target_name("muscarinic (m2) receptor")
    assert "m2" in res.clean_text.split("|")
    assert res.hints["mutations"] == []

    text = (
        "BRAF p.Q123* p.Q123Ter p.Q123_L125del p.T78_S79insA "
        "p.A100dup p.R97fs*5 p.Met1? p.*100Y p.K45delinsST"
    )
    res = normalize_target_name(text)
    assert res.clean_text == "braf"
    mutations = cast(List[str], res.hints["mutations"])
    assert set(mutations) == {
        "p.Q123*",
        "p.Q123Ter",
        "p.Q123_L125del",
        "p.T78_S79insA",
        "p.A100dup",
        "p.R97fs*5",
        "p.Met1?",
        "p.*100Y",
        "p.K45delinsST",
    }
