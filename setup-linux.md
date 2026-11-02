# Setup guide: Linux (bash)

Step-by-step instructions for getting `uoa-moa-catalog-notebook` running on Linux.
It takes about 15 minutes. For what the repo is and how it's laid out, see
[README.md](README.md). All commands below are run in a **bash** terminal.

Other platforms: [Windows](setup-windows.md) · [macOS](setup-mac.md)

---

## 1. Prerequisites

| Tool | Version | Check with |
|---|---|---|
| Python | **3.13+** | `python3.13 --version` |
| Git (optional) | any recent. Only needed to clone; you can download a ZIP instead | `git --version` |
| VS Code (recommended) | with the **Python** and **Jupyter** extensions | — |

> Don't have Python 3.13? Install it with your package manager (e.g. on Ubuntu,
> `sudo apt install python3.13 python3.13-venv` from the deadsnakes PPA).

> **Don't use an older Python** (e.g. a system 3.7). It will still create a
> venv, but the pinned dependencies are tested against 3.13 only.

## 2. Get the code

**Option A: Download the ZIP (no Git needed).** On the
[repo page](https://github.com/dsin867/uoa-moa-catalog-notebook) choose
**Code → Download ZIP**, unzip it, and open a terminal in the unzipped
`uoa-moa-catalog-notebook` folder.

**Option B: Clone with Git.**

```bash
git clone https://github.com/dsin867/uoa-moa-catalog-notebook.git
cd uoa-moa-catalog-notebook
```

## 3. Create the virtual environment

The venv **must** be called `.venv` and live at the repo root.
VS Code expects it there.

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python --version                  # should print 3.13.x
```

## 4. Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 5. Register the Jupyter kernel

```bash
python -m ipykernel install --user --name uoa-moa-catalog-notebook \
    --display-name "uoa-moa-catalog-notebook (Python 3.13)"
```

Check it with `jupyter kernelspec list`.

## 6. Get an access token

The notebook talks to the **MOA Catalog API** (`https://apis.auckland.ac.nz/vera-rubin/v1`).
No environment configuration is needed; you only need an access token from
<https://moacatalog.auckland.ac.nz>.

### Getting a token

1. Open <https://moacatalog.auckland.ac.nz> and sign in.
2. Copy the **access token** from the dashboard. Don't use the ID token: it has
   no `scope` claim, so the API rejects it.
3. Either paste it when the notebook's *Authentication* cell prompts you, or set
   it before starting Jupyter:
   ```bash
   export VERA_RUBIN_ACCESS_TOKEN=<token>
   ```

Tokens are short-lived. If you get a **401**, fetch a new token and re-run the
Authentication cell. Never paste a token into a code cell, since it would end up
in the notebook's source.

## 7. Run the notebook

1. Open the repo folder in VS Code and open `moa_catalog_search.ipynb`.
2. In the kernel picker (top-right), choose **uoa-moa-catalog-notebook (Python 3.13)**
   or `.venv`. Don't pick a global or system Python.
3. **Restart Kernel → Run All.**

With the defaults it should:

- accept your token in the *Authentication* cell;
- find the default target (in the SMC) observed in MOA (the coverage check can take up
  to ~3 min);
- show a cutout of the closest source, then submit a photometry job over the full
  archive and poll it every minute until it finishes (roughly 30 s per 1000
  observations; the default target is a small field, so it should finish
  much sooner than a dense bulge field) before drawing the
  light curve;
- write `outputs/my_source_<ra>_<dec>.fits` and
  `outputs/fluxes_<ra>_<dec>_moa.csv`.

To change the target or aperture, edit only the **Parameters** cell. If the
photometry wait stops early, the job keeps running server-side: set
`RESUME_JOB_ID` to the printed job id and re-run section 2b of the notebook.

Generated files go to `outputs/` (FITS cutouts and flux CSVs), and the HTTP
client's debug log to `lib/logs/client.log`.

---

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ModuleNotFoundError: numpy` / `astropy` in the notebook | Wrong kernel selected. Pick **uoa-moa-catalog-notebook (Python 3.13)** or `.venv`. |
| `ModuleNotFoundError: moa_catalog_client` | The bootstrap cell wasn't run, or the notebook was opened outside the repo. Run all cells from the top. |
| `401 Unauthorized` | The token expired or is an ID token. Get a fresh **access** token from <https://moacatalog.auckland.ac.nz> (step 6). |
| `429` on submitting photometry | Job submissions are rate-limited (3/minute, 15/hour). Wait, or resume an existing job with `RESUME_JOB_ID`. |
| Photometry stops waiting after 1 hour | The job is still running server-side. Set `RESUME_JOB_ID` to the printed job id and re-run section 2b of the notebook. |
| Target reports "not observed" or FITS errors | Check the *Example targets* table in the notebook. Some coordinates are known to be too faint or to have corrupt FITS files. |
| VS Code Pylance can't resolve `lib` imports | Open the **repo root** as the workspace. |
