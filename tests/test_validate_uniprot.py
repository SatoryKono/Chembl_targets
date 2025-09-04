from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library.validate import validate_uniprot_dataframe
from main import main as cli_main


def fake_record(_uniprot_id: str) -> dict:
    return {
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": "Hemoglobin subunit beta"}}
        },
        "genes": [{"geneName": {"value": "HBB"}, "synonyms": []}],
    }


def test_validate_uniprot_dataframe(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "library.validate.fetch_uniprot_record", lambda uid: fake_record(uid)
    )
    sample = Path("tests/data/uniprot_sample.csv")
    df = pd.read_csv(sample)
    result = validate_uniprot_dataframe(df, "uniprot_id", "protein_name")
    assert bool(result.loc[0, "uniprot_match"]) is True
    assert bool(result.loc[1, "uniprot_match"]) is False


def test_cli_validation(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "library.validate.fetch_uniprot_record", lambda uid: fake_record(uid)
    )
    sample = Path("tests/data/uniprot_sample.csv")
    out = tmp_path / "out.csv"
    argv = [
        "main.py",
        "--input",
        str(sample),
        "--output",
        str(out),
        "--column",
        "protein_name",
        "--id-column",
        "uniprot_id",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    cli_main()
    df = pd.read_csv(out)
    assert "uniprot_match" in df.columns
