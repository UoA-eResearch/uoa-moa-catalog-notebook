"""
HTTP client for the MOA Catalog API: coverage check (check_in_moa), async
photometry jobs (retrieve_photometry / get_photometry_job), and the client-side
helpers the notebook uses to parse, validate, plot and reload the results.
"""

import base64
import io
import logging
import time
import warnings

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.exceptions import AstropyUserWarning, AstropyWarning
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS

from logging_config import DEBUG_LOG

plt.rcParams['figure.figsize'] = (1.5,1.5)

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=AstropyUserWarning)
warnings.filterwarnings("ignore", message="invalid value encountered in subtract")

# Mirrors moaPhotometryUtils.MAX_BKG_OUT. Duplicated rather than imported so this
# client stays free of the backend's database/archive dependencies - keep the two
# in step if the patch size in generate_photometry ever changes.
MAX_BKG_OUT = 7.0

# Module logger writing to DEBUG_LOG. Scoped to this module rather than the root
# logger so other libraries' debug output stays out of the file, and guarded so
# %autoreload re-executing this module doesn't stack up duplicate handlers.
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
if not logger.handlers:
    fh = logging.FileHandler(DEBUG_LOG)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(fh)


# --- Configuration ---
# The MOA Catalog API, and the auth app that issues access tokens for it.
API_BASE_URL = "https://apis.auckland.ac.nz/vera-rubin/v1"
AUTH_APP_URL = "https://moacatalog.auckland.ac.nz"

# The process API only accepts percentage_of_data = 100 (api-vera-rubin-v1's
# retrieve-photometry-request.ts) - photometry always runs over the full archive,
# same as app-vera-rubin-auth's FIXED_PERCENTAGE_OF_DATA.
PERCENTAGE_OF_DATA = 100

# Photometry runs as an async job (POST /photometry_jobs, then poll its status),
# since a full-archive run takes minutes - far longer than one gateway request
# should be held open. Status polls are rate-limited at 60/minute.
JOB_POLL_INTERVAL = 60  # seconds between status polls
JOB_MAX_WAIT = 3600  # seconds to wait before giving up (the job keeps running server-side)
# (connect, read) timeout for each individual job-API / result-download request.
JOB_REQUEST_TIMEOUT = (10, 60)
# A status poll failing this many times in a row (network blip, 5xx) ends the wait.
JOB_MAX_POLL_FAILURES = 3

# Identifies this notebook client to the process API, whose audit trail classifies
# a request as jupyter_notebook by the "vera-rubin-notebook/" User-Agent prefix
# (api-vera-rubin-v1's context/RSM-5725-audit-log.md) - without it, calls look like
# any other python-requests script. Bump CLIENT_VERSION on notable client changes.
CLIENT_VERSION = "1.0.0"
USER_AGENT = f"vera-rubin-notebook/{CLIENT_VERSION} {requests.utils.default_user_agent()}"


class APIError(Exception):
    """Raised when API request fails"""

# --- Authentication ---
# The process API requires a Cognito bearer token with
# the selected environment's api-vera-rubin-v1 scope. Get one by opening that
# environment's auth app (AUTH_APP_URL), signing in, and copying the access token
# shown on its dashboard - see the notebook's Authentication cell and this repo's
# README. Held in memory only for the life of the kernel process, never written
# to disk.
_access_token: str | None = None


def set_access_token(token: str) -> None:
    """Store the Cognito access token (from the auth app at AUTH_APP_URL) used to
    authenticate calls to the process API. Call this before check_in_moa /
    retrieve_photometry."""
    global _access_token
    _access_token = token.strip()


def _authed_request(method: str, path: str, payload: dict | None = None,
                    timeout=JOB_REQUEST_TIMEOUT) -> requests.Response:
    """
    Send a request to API_BASE_URL + path with the access token set via
    set_access_token(). A 401 means the token is missing/expired/wrongly-scoped -
    surfaced with a message telling the user how to get a fresh one, rather than
    retried automatically (there's no automated login to retry with; a human has
    to go back to the auth app at AUTH_APP_URL).
    """
    if not _access_token:
        raise APIError(
            f"No access token set. Sign in via {AUTH_APP_URL}, copy the access token from its "
            "dashboard, and call set_access_token(<token>) before making API calls."
        )
    response = requests.request(
        method, f"{API_BASE_URL}{path}", json=payload, timeout=timeout,
        headers={"Authorization": f"Bearer {_access_token}", "User-Agent": USER_AGENT}
    )
    if response.status_code == 401:
        print(
            "[ERROR] API rejected the current access token (401) - it may be expired. Get a fresh one "
            f"from {AUTH_APP_URL} and call set_access_token(<token>) again."
        )
    return response

def check_in_moa(given_coords: SkyCoord) -> tuple[Cutout2D | None, SkyCoord | None, fits.Header | None]:
    """
    Check if the given coordinates have been observed by MOA

    Returns:
        Tuple of (image_cutout, source_position, cutout_header)
    """

    # Extract RA/Dec from SkyCoord
    ra_deg = given_coords.ra.degree
    dec_deg = given_coords.dec.degree
    frame = given_coords.frame.name.lower()

    # Prepare request - send directly as JSON
    payload = {
        "ra": str(ra_deg),
        "dec": str(dec_deg),
        "coord_frame": frame
    }

    print(f"Checking MOA observations for RA={ra_deg:.6f}, Dec={dec_deg:.6f}, Frame={frame} ...")
    print(f"Sending request with payload: {payload}")

    try:
        # Call API endpoint
        response = _authed_request("POST", "/check_in_moa", payload, timeout=300)
        response.raise_for_status()

        data = response.json()

        message = data.get("message", "")

        # found=true with no source_coordinates means the coordinates WERE observed and
        # the cutout loaded fine, but DAOStarFinder found no star clearing the 5-sigma
        # threshold. The return value is all-None either way, so this printed line is the
        # only place that distinction is visible to the caller — flag it as a warning
        # instead of letting it blend into the ordinary status line, and record it in
        # DEBUG_LOG for after-the-fact diagnosis.
        if data.get("found") and not data.get("source_coordinates"):
            detail = message or (
                "Coordinates were observed by MOA, but no source bright enough to clear "
                "the 5-sigma detection threshold was found in the 15 arcsec cutout."
            )
            print(f"[WARNING] MOA API: {detail}")
            logger.warning(f"check_in_moa: no source centroid at ra={ra_deg}, dec={dec_deg} — {detail}")
        elif message:
            print(f"MOA API: {message}")

        # Do not continue unless the API reports found=true AND provides a source
        # coordinate — without a centroid there is nothing usable to hand back.
        if not data.get("found") or not data.get("source_coordinates"):
            return None, None, None

        # Decode cutout FITS data - THIS IS THE KEY!
        cutout_data = data.get("cutout")
        if cutout_data and cutout_data.get("data"):
            # Get the base64 string from the cutout object
            fits_base64 = cutout_data["data"]

            # Decode base64 to bytes
            fits_bytes = base64.b64decode(fits_base64)

            # Open FITS from bytes - THIS CONTAINS THE COMPLETE HEADER!
            hdul = fits.open(io.BytesIO(fits_bytes))
            image_data = hdul[0].data
            full_header = hdul[0].header  # This is the COMPLETE fits.Header object

            # Create WCS from the full header
            wcs_obj = WCS(full_header)

            # Create Cutout2D object
            cutout = Cutout2D(
                data=image_data,
                position=(image_data.shape[1] // 2, image_data.shape[0] // 2),
                size=image_data.shape,
                wcs=wcs_obj
            )
            cutout.header = full_header
            hdul.close()

            header_dict = data.get("header_data")
            header_data = dict_to_fits_header(header_dict)
        else:
            cutout = None
            header_data = None

        # Reconstruct source coordinates from JSON (guaranteed present by the guard above)
        src_coords_data = data["source_coordinates"]
        source_position = SkyCoord(
            ra=src_coords_data["ra"] * u.deg,
            dec=src_coords_data["dec"] * u.deg,
            frame=src_coords_data["frame"]
        )

        return cutout, source_position, header_data  # Returns proper fits.Header!

    except (requests.exceptions.RequestException, APIError) as e:
        print(f"Error calling MOA API: {e}")
        return None, None, None

def dict_to_fits_header(header_dict):
    """Convert dictionary back to FITS Header object"""
    if header_dict is None:
        return None

    # Check if this is the ordered format
    if '_cards' in header_dict:
        # Reconstruct from ordered cards
        cards = []
        for card_dict in header_dict['_cards']:
            card = fits.Card(
                keyword=card_dict['keyword'],
                value=card_dict['value'],
                comment=card_dict.get('comment', '')
            )
            cards.append(card)
        return fits.Header(cards)

    # Fallback to old method for backward compatibility
    cards = []
    for key, value in header_dict.items():
        if isinstance(value, list):
            for item in value:
                cards.append((key, str(item)))
        else:
            cards.append((key, value))
    return fits.Header(cards)

def parse_coordinates(ra, dec, coord_frame="fk5"):
    """
    Parse coordinate strings and return SkyCoord object.
    This function is identical to the backend version (no API call needed).

    Parameters:
        ra: Right ascension (string)(h m s; h:m:s; float)
        dec: Declination (string)(h m s; h:m:s; float)
        coord_frame: Coordinate system (icrs, fk4, or fk5; default fk5)

    Return:
        Astropy SkyCoord object
    """
    logger.debug(f"parse_coordinates called with ra={ra}, dec={dec}, coord_frame={coord_frame}")

    coord_unit = u.deg  # assume float

    # Check if it's a float; if not, assume hourangle format
    try:
        float(ra)
        float(dec)
        logger.debug("Input parsed as decimal degrees")
    except ValueError:
        coord_unit = (u.hourangle, u.deg)
        logger.debug("Input parsed as hourangle/degrees format")

    coord_vals = SkyCoord(ra, dec, unit=coord_unit, frame=coord_frame)
    logger.debug(f"Created SkyCoord: RA={coord_vals.ra.degree:.4f}°, Dec={coord_vals.dec.degree:.4f}°, "
                 f"Frame={coord_vals.frame.name}")

    return coord_vals


def display_cutout(img_cutout, position_tuple):
    """
    Display 15 arcsec region around given coordinate.
    This function should be implemented locally on the client
    since visualization is client-side.

    Parameters:
        img_cutout: Cutout2D object of the select region
        position_tuple: central pixel coordinates of the object to highlight
    """

    inverval = ZScaleInterval()
    zmin, zmax = inverval.get_limits(img_cutout.data)

    plt.subplot(projection=img_cutout.wcs)
    ax = plt.gca()

    ax.coords['ra'].set_ticklabel_position('l')
    ax.coords['ra'].set_format_unit('degree')
    ax.coords['dec'].set_ticklabel_position('b')
    ax.coords['dec'].set_format_unit('degree')
    ax.set_ylabel('RA')
    ax.set_xlabel('Dec')
    ax.axis('off')
    ax.set_aspect('equal')

    circular_patch = mpatches.Circle(position_tuple, radius=5, edgecolor='red', facecolor='none')
    ax.add_patch(circular_patch)

    plt.imshow(img_cutout.data, origin='lower', vmin=zmin, vmax=zmax)


def define_annulus(img_cutout, center_tuple, star_radius, bkg_in, bkg_out):
    """
    Display aperture and background annuli.
    This function should be implemented locally on the client
    since visualization is client-side.

    Parameters:
        img_cutout: Cutout2D object of the select region
        center_tuple: central pixel coordinates of the object
        star_radius: aperture radius in pixels
        bkg_in: inner background annulus radius
        bkg_out: outer background annulus radius
    """

    inverval = ZScaleInterval()
    zmin, zmax = inverval.get_limits(img_cutout.data)

    _, ax = plt.subplots()

    ax.axis('off')

    star_patch = mpatches.Circle(center_tuple, radius=star_radius, edgecolor='red', facecolor='none')
    ax.add_patch(star_patch)
    bkg_in_patch = mpatches.Circle(center_tuple, radius=bkg_in, edgecolor='green', facecolor='none')
    ax.add_patch(bkg_in_patch)
    bkg_out_patch = mpatches.Circle(center_tuple, radius=bkg_out, edgecolor='green', facecolor='none')
    ax.add_patch(bkg_out_patch)

    plt.imshow(img_cutout.data, origin='lower', vmin=zmin, vmax=zmax)


def _report_failed_response(response: requests.Response) -> None:
    """Print a non-200 API response. Gateway errors (e.g. a 504 from Kong when
    the request runs too long) come back with a non-JSON body, so don't assume
    .json() will parse."""
    try:
        detail = response.json()
    except ValueError:
        detail = response.text.strip() or response.reason
    print(f"[ERROR] API request failed ({response.status_code}): {detail}")
    if response.status_code == 504:
        print("[HINT] The request timed out at the gateway - try again shortly.")
    elif response.status_code == 429:
        print("[HINT] Rate limit reached (job submissions: 3/minute, 15/hour) - wait and try again.")


def parse_photometry_data(photometry_data: str, data_format: str | None) -> np.ndarray:
    """
    Parse the photometry_data field of a /retrieve_photometry response into an
    (N, 3) array of [jd, flux, uncertainty] rows.

    With data_format "csv" (currently the only format the API returns),
    photometry_data is "jd,flux,uncertainty" rows joined by newlines, no header
    row. A failed-aperture observation's values are the literal "nan", which
    np.loadtxt parses back to NaN.

    Parameters:
        photometry_data: the response's photometry_data string
        data_format: the response's data_format field
    """
    if data_format != "csv":
        raise APIError(f"Unsupported photometry data_format: {data_format!r} (expected 'csv')")
    if not photometry_data.strip():
        return np.empty((0, 3))
    return np.loadtxt(io.StringIO(photometry_data), delimiter=',', ndmin=2)


def _validate_aperture(star_radius, bkg_in, bkg_out) -> str | None:
    """Return a message describing the first out-of-range aperture radius, or None if all are valid."""
    if star_radius > 5.0 or star_radius < 1.0:
        return "Please enter a star radius between 1.0 and 5.0"
    if bkg_in > bkg_out:
        return "bkg_in must be less than bkg_out"
    if bkg_in <= star_radius or bkg_in > 10.0:
        return "bkg_in should be greater than star radius and less than 10"
    if bkg_out <= star_radius or bkg_out > MAX_BKG_OUT:
        return f"bkg_out should be greater than the star radius and at most {MAX_BKG_OUT}"
    return None


def retrieve_photometry(given_coords, star_radius, bkg_in, bkg_out,
                        poll_interval=JOB_POLL_INTERVAL, max_wait=JOB_MAX_WAIT):
    """
        Run aperture photometry for the source nearest the given coordinates as a
        photometry job on the process API, and wait for its result.

        Mirrors check_in_moa's coordinate handling: the backend locates the
        reference header and detects the source itself from ra/dec/coord_frame,
        so no prior check_in_moa cutout or header_data is needed here.

        Submits POST /photometry_jobs, prints the job_id, then hands off to
        get_photometry_job() to poll and download the result - the same flow as
        app-vera-rubin-auth's job pages. Always runs over 100% of the archive
        (PERCENTAGE_OF_DATA), the only value the API accepts.

        If the wait ends early (max_wait reached, an expired token, or a kernel
        interrupt) the job keeps running server-side: call get_photometry_job(job_id)
        to pick it up rather than submitting a new one.

        Returns an (N, 3) numpy array, one row per observation: [mid-exposure JD,
        aperture sum (flux), uncertainty], or None on failure.

        Parameters:
            given_coords: SkyCoord with the sky coordinates to search
            star_radius: aperture radius in pixels; must be between 1.0 and 5.0
            bkg_in: inner background annulus radius; must be greater than star_radius and at most 10.0
            bkg_out: outer background annulus radius; must be at least bkg_in, greater than
                star_radius, and at most MAX_BKG_OUT (7.0)
            poll_interval: seconds between status polls; defaults to JOB_POLL_INTERVAL
            max_wait: seconds to wait for the job to finish; defaults to JOB_MAX_WAIT
    """
    problem = _validate_aperture(star_radius, bkg_in, bkg_out)
    if problem:
        print(problem)
        return None

    payload = {
        "ra": str(given_coords.ra.degree),
        "dec": str(given_coords.dec.degree),
        "coord_frame": given_coords.frame.name.lower(),
        "percentage_of_data": PERCENTAGE_OF_DATA,
        "star_radius": star_radius,
        "bkg_in": bkg_in,
        "bkg_out": bkg_out
    }

    try:
        response = _authed_request("POST", "/photometry_jobs", payload)
        if not response.ok:
            _report_failed_response(response)
            return None
        job_id = response.json()["job_id"]
    except (requests.exceptions.RequestException, APIError, ValueError, KeyError) as e:
        print(f"[ERROR] Failed to submit the photometry job: {e}")
        return None

    print(f"Submitted photometry job {job_id}. It can take several minutes "
          f"(~30 s per 1000 observations).")
    print(f"You can also follow this job at {AUTH_APP_URL}/jobs/{job_id}")
    return get_photometry_job(job_id, poll_interval=poll_interval, max_wait=max_wait)


def get_photometry_job(job_id, poll_interval=JOB_POLL_INTERVAL, max_wait=JOB_MAX_WAIT):
    """
    Wait for an already-submitted photometry job to finish and return its result.

    Use it to resume a job whose wait in retrieve_photometry ended early (max_wait
    reached, an expired token - call set_access_token() with a fresh one first - or
    a kernel interrupt) without submitting a new one. Polls
    GET /photometry_jobs/{job_id} every poll_interval seconds until the job is
    SUCCEEDED or FAILED, then downloads the result.

    Returns the same (N, 3) array as retrieve_photometry, or None on failure.

    Parameters:
        job_id: the job_id printed by retrieve_photometry
        poll_interval: seconds between status polls; defaults to JOB_POLL_INTERVAL
        max_wait: seconds to wait for the job to finish; defaults to JOB_MAX_WAIT
    """
    try:
        job = _wait_for_job(job_id, poll_interval, max_wait)
        if job is None:
            _print_resume_hint(job_id)
            return None
        if job["status"] == "FAILED":
            print(f"[ERROR] Photometry job {job_id} failed: {job.get('error') or 'no error message'}")
            return None
        return _download_job_result(job_id)
    except KeyboardInterrupt:
        print()
        _print_resume_hint(job_id)
        raise
    except (requests.exceptions.RequestException, APIError, ValueError, KeyError) as e:
        print(f"[ERROR] Failed to retrieve photometry for job {job_id}: {e}")
        _print_resume_hint(job_id)
        return None


def _print_resume_hint(job_id) -> None:
    print(f"[INFO] Job {job_id} may still be running server-side - resume with "
          f"get_photometry_job('{job_id}') instead of submitting a new job.")


class _TransientPollError(Exception):
    """A status poll failure worth retrying: a network error, 429 or 5xx."""


def _get_job_status(job_id) -> dict:
    """
    One GET /photometry_jobs/{job_id}. Raises _TransientPollError for failures
    worth retrying, and APIError for ones that aren't (401, 403, 404, ...).
    """
    try:
        response = _authed_request("GET", f"/photometry_jobs/{job_id}")
    except requests.exceptions.RequestException as e:
        raise _TransientPollError(str(e)) from e
    if response.ok:
        return response.json()
    if response.status_code == 429 or response.status_code >= 500:
        raise _TransientPollError(f"HTTP {response.status_code}")
    _report_failed_response(response)
    raise APIError(f"Status poll for job {job_id} failed ({response.status_code})")


def _wait_for_job(job_id, poll_interval, max_wait) -> dict | None:
    """
    Poll the job's status until it is SUCCEEDED or FAILED and return the final
    status body, or None (after printing why) if max_wait passes or polling fails
    JOB_MAX_POLL_FAILURES times in a row. Prints each status change, and a
    still-running line about once a minute so a long wait doesn't look hung.
    """
    start = time.monotonic()
    last_status, last_print, failures = None, start, 0
    while True:
        try:
            job = _get_job_status(job_id)
            failures = 0
        except _TransientPollError as e:
            failures += 1
            print(f"[WARN] Status poll failed ({e}) - attempt {failures} of {JOB_MAX_POLL_FAILURES}.")
            if failures >= JOB_MAX_POLL_FAILURES:
                return None
            job = None

        now = time.monotonic()
        if job is not None:
            if job["status"] != last_status or now - last_print >= 60:
                last_status, last_print = job["status"], now
                print(f"[{int(now - start) // 60}m{int(now - start) % 60:02d}s] Job {job_id}: {last_status}")
            if last_status in ("SUCCEEDED", "FAILED"):
                return job

        if now - start >= max_wait:
            print(f"[WARN] Job {job_id} has not finished after {max_wait} s - stopped waiting.")
            return None
        time.sleep(poll_interval)


def _download_job_result(job_id) -> np.ndarray:
    """Get a SUCCEEDED job's presigned result URL, download the result JSON and parse it."""
    response = _authed_request("GET", f"/photometry_jobs/{job_id}/result")
    if not response.ok:
        _report_failed_response(response)
        raise APIError(f"Could not get the result URL for job {job_id} ({response.status_code})")
    result_url = response.json()["result_url"]

    # A presigned S3 URL, not a process API endpoint - fetched without the bearer
    # token, since the URL's own signature is the authorisation (S3 rejects a
    # request carrying a second credential). Same as app-vera-rubin-auth's job
    # detail page.
    result = requests.get(result_url, timeout=JOB_REQUEST_TIMEOUT)
    result.raise_for_status()
    return _parse_photometry_result(result.json())


def _parse_photometry_result(resp_json: dict) -> np.ndarray:
    """
    Parse a photometry result body ({data_format, photometry_data,
    num_observations, unreadable_files}) into the (N, 3) array, printing a
    summary and any unreadable archive files.
    """
    if resp_json.get("photometry_data") is None:
        raise APIError("No photometry in the job result")

    photometry_data = parse_photometry_data(resp_json["photometry_data"], resp_json.get("data_format"))

    # Print observations summary
    print(f"{photometry_data.shape[0]} observations collated "
          f"(num_observations reported by the API: {resp_json.get('num_observations')}).")

    # Archive files the server had to skip because the FITS was empty/corrupt.
    # Listed (capped) so a bad file can be identified and fixed on the share.
    unreadable = resp_json.get("unreadable_files") or []
    if unreadable:
        print(f"[WARN] {len(unreadable)} observation file(s) were unreadable and were skipped:")
        for path in unreadable[:10]:
            print(f"  {path}")
        if len(unreadable) > 10:
            print(f"  ... and {len(unreadable) - 10} more")

    return photometry_data


def validate_photometry(photometry) -> None:
    """
    Sanity-check the array returned by retrieve_photometry before it is
    plotted or saved. Raises ValueError describing the first problem found.

    Parameters:
        photometry: 2D numpy array, column 0 = mid-exposure JD, column 1 = flux,
            column 2 (if present) = flux uncertainty
    """
    if photometry is None:
        raise ValueError("No photometry - retrieve_photometry failed; see the messages above.")
    photometry = np.asarray(photometry)
    if photometry.ndim != 2 or photometry.shape[1] < 2:
        raise ValueError(f"Expected an (N, 2) array of (JD, flux), got shape {photometry.shape}.")
    if photometry.shape[0] == 0:
        raise ValueError("Photometry array is empty - no observations were returned.")
    if not np.isfinite(photometry[:, 0]).all():
        raise ValueError("Photometry contains non-finite Julian dates.")


def plot_lightcurve(photometry, ax=None):
    """
    Scatter-plot flux against days since the first observation.

    Parameters:
        photometry: 2D numpy array, column 0 = mid-exposure JD, column 1 = flux
        ax: optional matplotlib Axes to draw on; a new 10x5 figure is created if omitted
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))
    t0 = np.min(photometry[:, 0])
    ax.scatter(photometry[:, 0] - t0, photometry[:, 1], marker='d')
    ax.set_ylabel('Aperture Sum')
    ax.set_xlabel(f'Time (Days since JD {t0:.4f})')
    return ax


def load_photometry_csv(content: bytes):
    """
    Parse a fluxes CSV previously written by the notebook back into a 2D array.

    Accepts files with or without a "jd_mid,aperture_sum,uncertainty" header row
    (older files may have "jd_mid,aperture_sum" or "col_2"), and
    normalises line endings: csv.writer on Windows text-mode files can produce
    doubled "\\r\\r\\n" terminators that trip up np.loadtxt.

    Parameters:
        content: raw bytes of the CSV file
    """
    text = bytes(content).decode('utf-8').replace('\r\n', '\n').replace('\r', '\n')
    lines = [line for line in text.split('\n') if line.strip()]
    if lines and not lines[0].lstrip()[:1].isdigit():
        lines = lines[1:]  # skip header row
    return np.loadtxt(io.StringIO('\n'.join(lines)), delimiter=',', ndmin=2)
