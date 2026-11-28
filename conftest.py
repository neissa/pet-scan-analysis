import pytest

import store


@pytest.fixture(autouse=True)
def _data_dir(tmp_path, monkeypatch):
    """Les tests n'écrivent jamais dans le vrai dossier de données de l'utilisateur."""
    monkeypatch.setattr(store, "DIR", tmp_path / "netscan")
    monkeypatch.setattr(store, "LEGACY", [tmp_path / "legacy.json"])
