# uoa-moa-catalog-notebook

Jupyter notebook + HTTP client library for searching the MOA (Microlensing
Observations in Astrophysics) image catalog by sky coordinate, viewing a
FITS cutout, and running aperture photometry, all via HTTP against the
**MOA Catalog API**. No archive or database access is needed locally; only
network access to the API and a Cognito access token.

## Layout

```
lib/
├── moa_catalog_client.py   # HTTP client: check_in_moa, retrieve_photometry, parse_coordinates, ...
└── logging_config.py       # Logging setup (writes to lib/logs/, gitignored)
tests/                      # pytest suite for lib/ (no network or token needed)
moa_catalog_search.ipynb
requirements.txt            # Dependencies, pinned
requirements-dev.txt        # + pytest, ruff, nbstripout, pre-commit
pyproject.toml              # ruff and pytest configuration
outputs/                    # Saved FITS cutouts / flux CSVs (gitignored)
```

## Setup

Requires **Python 3.13** (or newer). Follow the step-by-step guide for your
operating system:

| OS | Guide |
|---|---|
| Windows (PowerShell) | [setup-windows.md](setup-windows.md) |
| macOS (Terminal) | [setup-mac.md](setup-mac.md) |
| Linux | [setup-linux.md](setup-linux.md) |

## Run the notebook

Open `moa_catalog_search.ipynb` and select the
**uoa-moa-catalog-notebook (Python 3.13)** kernel (VS Code: kernel picker in the
top-right; classic Jupyter: Kernel → Change Kernel). If numpy/astropy
imports fail with `ModuleNotFoundError`, the wrong kernel is selected.
Pick the one registered during setup, not a system/global Python.

Everything you'd normally change (target RA/Dec, aperture radii, output
directory, a job id to resume) lives in the single
**Parameters** cell near the top. The notebook is meant to run cleanly with
*Restart Kernel and Run All*: if the target isn't observed (or the request
fails) it stops at the coverage cell with a clear error.

## MOA Catalog API

The notebook uses `https://apis.auckland.ac.nz/vera-rubin/v1`. Get your access
token from <https://moacatalog.auckland.ac.nz>. No environment configuration is
needed. API documentation: <https://moacatalog.auckland.ac.nz/docs>.

Photometry runs as an async job: `retrieve_photometry` submits
`POST /photometry_jobs`, polls the job's status every minute, then downloads the
result. It always covers 100% of the archive, the only value the API accepts.
If the wait ends early (after 1 hour, an expired token, or a kernel
interrupt) the job keeps running server-side. Collect it later with
`get_photometry_job(job_id)`, or set `RESUME_JOB_ID` in the notebook's
Parameters cell. You can also follow a submitted job in the MOA Catalog app at
`https://moacatalog.auckland.ac.nz/jobs/<job_id>`.

## Authentication

Calls to the MOA Catalog API require a Cognito **access token**.

1. Open <https://moacatalog.auckland.ac.nz> and sign in.
2. On its dashboard, copy the **access token** (there's a Copy button), not
   the ID token, which carries no `scope` claim and will be rejected.
3. In the notebook's "Authentication" cell, run it and paste the token at the
   prompt. It's held in memory for the life of the kernel session only and is
   never written to disk.

Access tokens are short-lived. If a call fails with a 401, repeat the steps
above to get a fresh token and re-run the Authentication cell.

For headless runs, set the `VERA_RUBIN_ACCESS_TOKEN` environment variable
instead; the Authentication cell uses it when present and only prompts
otherwise. Never paste the token into a code cell.

## Development

Install the dev tools into the same `.venv` and enable the git hook once per
clone:

```bash
python -m pip install -r requirements-dev.txt
pre-commit install
```

The hook strips notebook outputs (`nbstripout`) and lints `lib/`, `tests/` and
the notebook (`ruff`) on every commit. Run the checks by hand with:

```bash
pre-commit run --all-files
pytest
```

CI (`.github/workflows/ci.yml`) runs the same lint, output check and tests on
every push to `main` and on pull requests.
