import httpx
import pytest

from library.uniprot import _CACHE, fetch_uniprot_records


@pytest.mark.asyncio
async def test_fetch_uniprot_records(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    async def fake_get(self, url: str, timeout: int = 10) -> httpx.Response:  # type: ignore[override]
        calls.append(url)
        uid = url.split("/")[-1].split(".")[0]
        data = {
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": f"Protein {uid}"}}
            },
            "genes": [],
        }
        req = httpx.Request("GET", url)
        return httpx.Response(200, json=data, request=req)

    monkeypatch.setattr(httpx.AsyncClient, "get", fake_get)
    result = await fetch_uniprot_records(["P1", "P2", "P1"])
    assert set(result) == {"P1", "P2"}
    assert len(calls) == 2  # duplicate ID should not trigger extra request
    # cache should be populated
    assert "P1" in _CACHE and "P2" in _CACHE
