"""Tests for the HTTP flows, with requests replaced by fakes - no network or token needed."""

import base64
import importlib
import io
import logging
import os

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from conftest import FakeResponse

import moa_catalog_client as client


def _fake_requests(monkeypatch, responses):
    """Route requests.request through a {(method, path): [responses...]} table, recording each call."""
    calls = []

    def fake_request(method, url, json=None, timeout=None, headers=None):
        path = url.removeprefix(client.API_BASE_URL)
        calls.append((method, path, json, headers))
        queue = responses[(method, path)]
        return queue.pop(0) if len(queue) > 1 else queue[0]

    monkeypatch.setattr(client.requests, "request", fake_request)
    return calls


def test_set_access_token_strips_whitespace(monkeypatch):
    monkeypatch.setattr(client, "_access_token", None)
    client.set_access_token("  abc.def.ghi \n")
    assert client._access_token == "abc.def.ghi"


def test_request_without_token_raises(monkeypatch):
    monkeypatch.setattr(client, "_access_token", None)
    with pytest.raises(client.APIError, match="No access token"):
        client._authed_request("GET", "/photometry_jobs/abc")


def test_bearer_token_is_sent(monkeypatch, token):
    calls = _fake_requests(monkeypatch, {("GET", "/ping"): [FakeResponse(200, {})]})
    client._authed_request("GET", "/ping")
    assert calls[0][3] == {"Authorization": "Bearer test-token", "User-Agent": client.USER_AGENT}


def _client_log_handlers(logger):
    log_path = os.path.abspath(client.DEBUG_LOG)
    return [h for h in logger.handlers if getattr(h, "baseFilename", None) == log_path]


def test_reload_does_not_duplicate_log_handlers():
    importlib.reload(client)
    importlib.reload(client)
    assert len(_client_log_handlers(logging.getLogger(client.__name__))) == 1
    assert _client_log_handlers(logging.getLogger()) == []


def _fits_cutout_b64(size=21):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [11.842125, -73.115006]
    wcs.wcs.crpix = [size / 2 + 0.5, size / 2 + 0.5]
    wcs.wcs.cdelt = [-0.58 / 3600, 0.58 / 3600]
    buffer = io.BytesIO()
    fits.PrimaryHDU(data=np.arange(size * size, dtype=float).reshape(size, size), header=wcs.to_header()).writeto(
        buffer
    )
    return base64.b64encode(buffer.getvalue()).decode()


class TestCheckInMoa:
    target = client.parse_coordinates("00 47 22.11", "-73 06 54.02", "fk5")

    def test_observed(self, monkeypatch, token):
        body = {
            "found": True,
            "message": "Observed in smc1-2",
            "cutout": {"data": _fits_cutout_b64()},
            "header_data": {"_cards": [{"keyword": "FIELD", "value": "smc1"}]},
            "source_coordinates": {"ra": 11.8422, "dec": -73.1150, "frame": "fk5"},
        }
        calls = _fake_requests(monkeypatch, {("POST", "/check_in_moa"): [FakeResponse(200, body)]})

        cutout, source, header = client.check_in_moa(self.target)

        assert calls[0][2]["coord_frame"] == "fk5"
        assert cutout.data.shape == (21, 21)
        assert cutout.wcs.wcs.ctype[0] == "RA---TAN"
        assert source.ra.degree == pytest.approx(11.8422)
        assert header["FIELD"] == "smc1"

    def test_observed_but_no_source_warns(self, monkeypatch, token, capsys):
        _fake_requests(monkeypatch, {("POST", "/check_in_moa"): [FakeResponse(200, {"found": True})]})
        assert client.check_in_moa(self.target) == (None, None, None)
        assert "[WARNING]" in capsys.readouterr().out

    def test_not_observed(self, monkeypatch, token):
        body = {"found": False, "message": "Not observed"}
        _fake_requests(monkeypatch, {("POST", "/check_in_moa"): [FakeResponse(200, body)]})
        assert client.check_in_moa(self.target) == (None, None, None)

    def test_expired_token(self, monkeypatch, token, capsys):
        _fake_requests(monkeypatch, {("POST", "/check_in_moa"): [FakeResponse(401, {"message": "Unauthorized"})]})
        assert client.check_in_moa(self.target) == (None, None, None)
        assert "401" in capsys.readouterr().out


class TestPhotometryJobs:
    job_id = "job-123"
    status_path = f"/photometry_jobs/{job_id}"
    result_path = f"/photometry_jobs/{job_id}/result"
    result_body = {"data_format": "csv", "photometry_data": "2450000.5,100.0,1.5\n2450001.5,101.0,1.6"}

    def _fake_s3(self, monkeypatch):
        downloads = []

        def fake_get(url, timeout=None):
            downloads.append(url)
            return FakeResponse(200, self.result_body)

        monkeypatch.setattr(client.requests, "get", fake_get)
        return downloads

    def test_submit_poll_and_download(self, monkeypatch, token, no_sleep):
        calls = _fake_requests(
            monkeypatch,
            {
                ("POST", "/photometry_jobs"): [FakeResponse(202, {"job_id": self.job_id})],
                ("GET", self.status_path): [
                    FakeResponse(200, {"status": "PENDING"}),
                    FakeResponse(200, {"status": "RUNNING"}),
                    FakeResponse(200, {"status": "SUCCEEDED"}),
                ],
                ("GET", self.result_path): [FakeResponse(200, {"result_url": "https://s3.example/result.json"})],
            },
        )
        downloads = self._fake_s3(monkeypatch)

        source = client.parse_coordinates("11.8422", "-73.1150", "fk5")
        photometry = client.retrieve_photometry(source, star_radius=1.0, bkg_in=1.5, bkg_out=2.0, poll_interval=0)

        assert photometry.shape == (2, 3)
        submitted = calls[0][2]
        assert submitted["percentage_of_data"] == 100
        assert (submitted["star_radius"], submitted["bkg_in"], submitted["bkg_out"]) == (1.0, 1.5, 2.0)
        assert [c[1] for c in calls].count(self.status_path) == 3
        assert downloads == ["https://s3.example/result.json"]

    def test_invalid_aperture_makes_no_request(self, monkeypatch, token, capsys):
        calls = _fake_requests(monkeypatch, {})
        source = client.parse_coordinates("11.8422", "-73.1150", "fk5")
        assert client.retrieve_photometry(source, star_radius=0.5, bkg_in=1.5, bkg_out=2.0) is None
        assert calls == []
        assert "star radius" in capsys.readouterr().out

    def test_rate_limited_submission(self, monkeypatch, token, capsys):
        _fake_requests(monkeypatch, {("POST", "/photometry_jobs"): [FakeResponse(429, {"message": "Too Many"})]})
        source = client.parse_coordinates("11.8422", "-73.1150", "fk5")
        assert client.retrieve_photometry(source, star_radius=1.0, bkg_in=1.5, bkg_out=2.0) is None
        assert "Rate limit" in capsys.readouterr().out

    def test_failed_job(self, monkeypatch, token, no_sleep, capsys):
        _fake_requests(
            monkeypatch,
            {("GET", self.status_path): [FakeResponse(200, {"status": "FAILED", "error": "Empty or corrupt FITS"})]},
        )
        assert client.get_photometry_job(self.job_id, poll_interval=0) is None
        assert "Empty or corrupt FITS" in capsys.readouterr().out

    def test_transient_poll_failures_give_up_with_resume_hint(self, monkeypatch, token, no_sleep, capsys):
        unavailable = FakeResponse(503, text="Service Unavailable")
        calls = _fake_requests(monkeypatch, {("GET", self.status_path): [unavailable]})
        assert client.get_photometry_job(self.job_id, poll_interval=0) is None
        assert len(calls) == client.JOB_MAX_POLL_FAILURES
        assert f"get_photometry_job('{self.job_id}')" in capsys.readouterr().out

    def test_transient_failure_then_recovery(self, monkeypatch, token, no_sleep):
        _fake_requests(
            monkeypatch,
            {
                ("GET", self.status_path): [
                    FakeResponse(502, text="Bad Gateway"),
                    FakeResponse(200, {"status": "SUCCEEDED"}),
                ],
                ("GET", self.result_path): [FakeResponse(200, {"result_url": "https://s3.example/result.json"})],
            },
        )
        self._fake_s3(monkeypatch)
        assert client.get_photometry_job(self.job_id, poll_interval=0).shape == (2, 3)

    def test_max_wait_stops_waiting(self, monkeypatch, token, no_sleep, capsys):
        _fake_requests(monkeypatch, {("GET", self.status_path): [FakeResponse(200, {"status": "RUNNING"})]})
        assert client.get_photometry_job(self.job_id, poll_interval=0, max_wait=0) is None
        out = capsys.readouterr().out
        assert "stopped waiting" in out
        assert "may still be running server-side" in out

    def test_result_download_omits_bearer_token(self, monkeypatch, token, no_sleep):
        _fake_requests(
            monkeypatch,
            {
                ("GET", self.status_path): [FakeResponse(200, {"status": "SUCCEEDED"})],
                ("GET", self.result_path): [FakeResponse(200, {"result_url": "https://s3.example/result.json"})],
            },
        )
        seen = {}

        def fake_get(url, timeout=None, **kwargs):
            seen.update(kwargs)
            return FakeResponse(200, self.result_body)

        monkeypatch.setattr(client.requests, "get", fake_get)
        client.get_photometry_job(self.job_id, poll_interval=0)
        assert "headers" not in seen
