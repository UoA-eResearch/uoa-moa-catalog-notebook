import os
import tempfile
from pathlib import Path

import matplotlib

# Must happen before moa_catalog_client is imported: keep the test run's debug log
# out of lib/logs/, and draw figures without a display.
os.environ.setdefault("MOA_CLIENT_LOG", str(Path(tempfile.mkdtemp()) / "client.log"))
matplotlib.use("Agg")

import pytest  # noqa: E402

import moa_catalog_client as client  # noqa: E402


class FakeResponse:
    """Minimal stand-in for requests.Response."""

    def __init__(self, status_code=200, json_body=None, text=""):
        self.status_code = status_code
        self._json = json_body
        self.text = text
        self.reason = "Fake"

    @property
    def ok(self):
        return self.status_code < 400

    def json(self):
        if self._json is None:
            raise ValueError("No JSON body")
        return self._json

    def raise_for_status(self):
        if not self.ok:
            import requests

            raise requests.exceptions.HTTPError(f"HTTP {self.status_code}")


@pytest.fixture
def token(monkeypatch):
    monkeypatch.setattr(client, "_access_token", "test-token")


@pytest.fixture
def no_sleep(monkeypatch):
    monkeypatch.setattr(client.time, "sleep", lambda _seconds: None)
