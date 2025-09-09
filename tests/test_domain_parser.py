"""Tests for domain parsing utilities."""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library.domain_parser import parse_domain


def _extract(result: Any) -> tuple[Any, list[str], Any]:
    """Helper to extract fields of interest."""
    return result.domain_type, result.domain_loc, result.domain_index


@pytest.mark.parametrize(
    "text,expected_type,expected_loc,expected_index",
    [
        ("ace c-domain", None, ["C_TERM"], None),
        ("ace n-domain", None, ["N_TERM"], None),
        ("alk cytoplasmic domain", None, ["ICD"], None),
        ("alk kinase domain", "KD", [], None),
        ("androgen receptor lbd domain", "LBD", [], None),
        ("ar ligand binding domain", "LBD", [], None),
        ("brd4 bromo domain 1", "BROMO_BD1", [], 1),
        ("brd4 bromodomain2", "BROMO_BD2", [], 2),
        ("brd4 tandem domain", "BROMO_TANDEM", [], [1, 2]),
        (
            "egfr cytoplasmic domain l858r / t790m mutant",
            None,
            ["ICD"],
            None,
        ),
    ],
)
def test_parse_domain(
    text: str, expected_type: Any, expected_loc: list[str], expected_index: Any
) -> None:
    result = parse_domain(text)
    domain_type, domain_loc, domain_index = _extract(result)
    assert domain_type == expected_type
    assert domain_loc == expected_loc
    assert domain_index == expected_index
