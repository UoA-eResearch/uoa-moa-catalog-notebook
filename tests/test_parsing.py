import numpy as np
import pytest
from astropy.io import fits

import moa_catalog_client as client


class TestParseCoordinates:
    def test_decimal_degrees(self):
        coords = client.parse_coordinates("268.1204134350", "-28.8896223900", "fk5")
        assert coords.ra.degree == pytest.approx(268.1204134350)
        assert coords.dec.degree == pytest.approx(-28.8896223900)
        assert coords.frame.name == "fk5"

    def test_sexagesimal_ra_is_hourangle(self):
        coords = client.parse_coordinates("00 47 22.11", "-73 06 54.02", "icrs")
        assert coords.ra.degree == pytest.approx((47 * 60 + 22.11) / 3600 * 15)
        assert coords.dec.degree == pytest.approx(-(73 + 6 / 60 + 54.02 / 3600))
        assert coords.frame.name == "icrs"


class TestParsePhotometryData:
    def test_csv_rows_with_nan(self):
        data = client.parse_photometry_data("2450000.5,100.0,1.5\n2450001.5,nan,nan", "csv")
        assert data.shape == (2, 3)
        assert data[0].tolist() == [2450000.5, 100.0, 1.5]
        assert np.isnan(data[1, 1:]).all()

    def test_single_row_stays_2d(self):
        assert client.parse_photometry_data("2450000.5,100.0,1.5", "csv").shape == (1, 3)

    def test_empty(self):
        assert client.parse_photometry_data("  \n", "csv").shape == (0, 3)

    def test_unsupported_format(self):
        with pytest.raises(client.APIError, match="data_format"):
            client.parse_photometry_data("1,2,3", "json")


class TestLoadPhotometryCsv:
    def test_with_header(self):
        content = b"jd_mid,aperture_sum,uncertainty\n2450000.5,100.0,1.5\n2450001.5,101.0,1.6\n"
        assert client.load_photometry_csv(content).shape == (2, 3)

    def test_without_header(self):
        assert client.load_photometry_csv(b"2450000.5,100.0,1.5\n").shape == (1, 3)

    def test_old_two_column_file(self):
        content = b"jd_mid,aperture_sum\n2450000.5,100.0\n2450001.5,101.0\n"
        assert client.load_photometry_csv(content).shape == (2, 2)

    def test_doubled_windows_line_endings(self):
        content = b"jd_mid,aperture_sum,uncertainty\r\r\n2450000.5,100.0,1.5\r\r\n2450001.5,101.0,1.6\r\r\n"
        data = client.load_photometry_csv(content)
        assert data.shape == (2, 3)
        assert data[1, 1] == 101.0


class TestValidateAperture:
    @pytest.mark.parametrize("radii", [(1.0, 1.5, 2.0), (5.0, 6.0, 7.0), (2.0, 3.0, 3.0)])
    def test_valid(self, radii):
        assert client._validate_aperture(*radii) is None

    @pytest.mark.parametrize(
        "radii, message",
        [
            ((0.5, 1.5, 2.0), "star radius"),
            ((5.5, 6.0, 7.0), "star radius"),
            ((2.0, 4.0, 3.0), "bkg_in must be less than bkg_out"),
            ((2.0, 2.0, 3.0), "bkg_in should be greater"),
            ((2.0, 3.0, client.MAX_BKG_OUT + 1), "bkg_out should be"),
        ],
    )
    def test_invalid(self, radii, message):
        assert message in client._validate_aperture(*radii)


class TestDictToFitsHeader:
    def test_none(self):
        assert client.dict_to_fits_header(None) is None

    def test_ordered_cards(self):
        header = client.dict_to_fits_header(
            {
                "_cards": [
                    {"keyword": "FIELD", "value": "gb5", "comment": "MOA field"},
                    {"keyword": "CHIP", "value": 2},
                ]
            }
        )
        assert isinstance(header, fits.Header)
        assert list(header.keys()) == ["FIELD", "CHIP"]
        assert header["CHIP"] == 2
        assert header.comments["FIELD"] == "MOA field"

    def test_legacy_flat_dict(self):
        header = client.dict_to_fits_header({"FIELD": "gb5", "HISTORY": ["one", "two"]})
        assert header["FIELD"] == "gb5"
        assert list(header["HISTORY"]) == ["one", "two"]


class TestValidatePhotometry:
    def test_valid(self):
        client.validate_photometry(np.array([[2450000.5, 100.0, 1.5], [2450001.5, np.nan, np.nan]]))

    @pytest.mark.parametrize(
        "photometry, message",
        [
            (None, "retrieve_photometry failed"),
            (np.array([1.0, 2.0]), "Expected an"),
            (np.empty((0, 3)), "empty"),
            (np.array([[np.nan, 100.0, 1.5]]), "non-finite Julian dates"),
        ],
    )
    def test_invalid(self, photometry, message):
        with pytest.raises(ValueError, match=message):
            client.validate_photometry(photometry)


class TestParsePhotometryResult:
    def test_missing_photometry(self):
        with pytest.raises(client.APIError, match="No photometry"):
            client._parse_photometry_result({"data_format": "csv"})

    def test_lists_unreadable_files_capped(self, capsys):
        result = client._parse_photometry_result(
            {
                "data_format": "csv",
                "photometry_data": "2450000.5,100.0,1.5",
                "num_observations": 1,
                "unreadable_files": [f"/archive/file{i}.fits" for i in range(12)],
            }
        )
        assert result.shape == (1, 3)
        out = capsys.readouterr().out
        assert "12 observation file(s) were unreadable" in out
        assert "/archive/file9.fits" in out
        assert "/archive/file10.fits" not in out
        assert "... and 2 more" in out


def test_plot_lightcurve_uses_days_since_first_observation():
    ax = client.plot_lightcurve(np.array([[2450001.5, 101.0, 1.0], [2450000.5, 100.0, 1.0]]))
    assert "JD 2450000.5000" in ax.get_xlabel()
    assert sorted(ax.collections[0].get_offsets()[:, 0].tolist()) == [0.0, 1.0]
