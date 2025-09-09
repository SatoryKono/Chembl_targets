from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from main import add_domain_columns


def test_add_domain_columns() -> None:
    df = pd.DataFrame({"domain": ["ace c-domain", "alk kinase domain"]})
    out = add_domain_columns(df, "domain")
    assert list(out["domain_loc"]) == ["C_TERM", None]
    assert list(out["domain_type"]) == [None, "KD"]
