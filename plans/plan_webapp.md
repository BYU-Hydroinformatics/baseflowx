# Baseflow Explorer — Web App Plan

## Names

- **Package:** `baseflowx` (renamed from `pybaseflow`)
- **Web app:** **Baseflow Explorer**
- **Public URL target:** `baseflow-explorer.streamlit.app`

Web app name alternates considered: *BaseflowLab*, *HydroSep*, *SepCompare*. Going with Baseflow Explorer — clear to hydrologists and students, hints at the interactive "pick a gage, see what happens" UX.

## Goal

A free, public web app that lets anyone:
1. Enter a USGS stream gage ID (or pick from a short curated list).
2. Pull daily streamflow from NWIS for a chosen date range.
3. Run every separation algorithm in `baseflowx` and plot the results together.
4. Download results as CSV.

The app is a living demo — when `baseflowx` ships a new method on PyPI, a redeploy picks it up automatically.

## Hosting

**Streamlit Community Cloud.** Free, unlimited public apps, auto-redeploys on `git push`, custom `*.streamlit.app` subdomain. Apps sleep after ~7 days idle (30s cold start on wake) — acceptable for a docs-linked demo.

Fallback if sleep becomes a problem: Hugging Face Spaces (Streamlit SDK, no sleep, 16GB RAM free).

## Architecture

```
USGS NWIS  ──dataretrieval──>  Streamlit app  ──baseflowx──>  Plotly chart + CSV
                                      │
                                      └── @st.cache_data on NWIS pulls
```

- **Data source:** `dataretrieval.nwis.get_dv` (daily values service). No API key needed.
- **Compute:** in-process call to `baseflowx`. For daily series (<30 years = ~11k points) every method runs well under a second even without numba.
- **Caching:** `@st.cache_data` on the NWIS fetch keyed by `(site_id, start, end)` so the same gage doesn't re-download.
- **Plotting:** Plotly (interactive zoom/pan matters for hydrographs) over Matplotlib.

---

## Prerequisite 1: rename `pybaseflow` → `baseflowx`

PyPI does **not** allow renaming an existing project — you publish under the new name and retire the old one. Steps:

### GitHub repo

1. Decide the new repo slug: `baseflowx` (matches package, clean URL).
2. On GitHub: **Settings → General → Rename repository** to `baseflowx`. GitHub sets up automatic redirects from `njones61/pybaseflow` → `njones61/baseflowx` for clones, issues, and PRs — safe.
3. Update local remote: `git remote set-url origin git@github.com:njones61/baseflowx.git`.
4. Rename the working directory locally: `/Users/njones/python_projects/pybaseflow` → `/Users/njones/python_projects/baseflowx`. (Optional but tidy.)

### Python package

1. Rename the package directory: `pybaseflow/` → `baseflowx/`.
2. Update `pyproject.toml`:
   - `name = "baseflowx"`
   - `[tool.setuptools.packages.find] include = ["baseflowx*"]`
   - `[tool.setuptools.package-data] baseflowx = ["data/*.csv", "data/*.md"]`
   - `[tool.setuptools.dynamic] version = {attr = "baseflowx.__version__"}`
   - Bump version to `0.2.0` (name change is a breaking API change — import path differs).
3. Grep the entire repo for `pybaseflow` and replace. Touches at minimum:
   - `baseflowx/__init__.py` (if it self-references)
   - any internal `from pybaseflow...` imports
   - `README.md` (install snippet, import examples, badges)
   - `docs/**/*.md` (all MkDocs pages)
   - `mkdocs.yml` (`site_name`, `repo_url`, `repo_name`, any edit_uri)
   - `pybaseflow.egg-info/` — just delete, it regenerates
   - test files if any reference the package name
4. Rebuild and run the test suite.

### PyPI

1. Register and publish `baseflowx 0.2.0` (fresh project under the new name — new PyPI project, you're the sole owner).
2. On the old `pybaseflow` PyPI project, publish **one** final version `pybaseflow 0.2.0` as a thin redirect shim:
   - Description prominently says "Renamed to baseflowx. Install `baseflowx` instead."
   - `pyproject.toml` lists `baseflowx>=0.2.0` as its only runtime dependency.
   - `pybaseflow/__init__.py` does `from baseflowx import *` and emits a `DeprecationWarning` on import.
   - This keeps every existing `pip install pybaseflow` working while nudging users off it.
3. Do **not** yank older `pybaseflow` releases — breaks anyone pinning a specific version. Leave them up with the deprecation notice on the project page.

### Docs / site

1. Update `mkdocs.yml` site name, repo URL, and any hard-coded install commands.
2. If the MkDocs site is deployed to GitHub Pages on `njones61.github.io/pybaseflow/`, GitHub's repo redirect does **not** cover Pages — the new Pages URL becomes `njones61.github.io/baseflowx/`. Update any external links (from your personal site, the PyPI description, etc.) after the rename.
3. Update README badges (PyPI, CI, docs link) to point to the new project.

### Verification checklist

- [ ] `pip install baseflowx` works in a fresh venv; `import baseflowx` succeeds.
- [ ] `pip install pybaseflow` still works and imports without error (deprecation warning visible).
- [ ] MkDocs site builds and deploys under the new URL.
- [ ] `git clone` of the old repo URL still succeeds (via GitHub redirect).
- [ ] All tests pass.

## Prerequisite 2: make numba optional

Blocker for clean cloud deploys and for an eventual Pyodide port. Small change:

- `baseflowx/separation.py`, `estimate.py`, `utils.py`: replace `from numba import njit, prange` with a try/except shim — `njit` becomes an identity decorator and `prange = range` when numba is absent.
- `pyproject.toml`: move `numba` from `dependencies` to `[project.optional-dependencies].speed`.
- Document `pip install baseflowx[speed]` in README for users who want the JIT.
- Ship as `baseflowx 0.2.1`.

Verify: existing test suite passes with `numba` uninstalled.

---

## UI sketch

Sidebar:
- Gage selector: text input for site ID **plus** a dropdown of ~10 curated gages with known-interesting baseflow behavior (e.g. one spring-fed, one flashy, one snowmelt-dominated). BN77 from the docs example should be in there.
- Date range picker (default: last 5 water years).
- Method multi-select (default: all). Expose method-specific params behind an expander.

Main pane:
- Summary card: site name, drainage area, record period, mean Q.
- Hydrograph: log-scale Q on y-axis, one line per selected method + the raw streamflow.
- BFI table: baseflow index per method.
- Download buttons: CSV of all separations, PNG of plot.

## Phases

1. **Rename `pybaseflow` → `baseflowx`** (prerequisite 1 above). Ship `baseflowx 0.2.0` + deprecation shim `pybaseflow 0.2.0`.
2. **Numba-optional refactor** (prerequisite 2 above). Ship `baseflowx 0.2.1`.
3. **Scaffold Streamlit app** in a `/webapp` folder of the `baseflowx` repo.
   - `streamlit_app.py`, `requirements.txt` pinning `baseflowx>=0.2.1`.
   - Minimal working version: one gage, one method, one chart.
4. **Fill in methods + params** — wire up the full method set with per-method param controls.
5. **Polish** — curated gage list, site metadata lookup, download buttons, BFI table.
6. **Deploy to Streamlit Cloud** — connect repo, claim `baseflow-explorer.streamlit.app`, add badge in docs.
7. **Link from docs** — "Try it" callout on the docs landing page and each method page.
8. **(Later) Pyodide demo** — embed a lightweight version directly in MkDocs pages so readers can run an example without leaving the docs.

## Open questions

- **Repo layout:** `/webapp` subfolder in `baseflowx` vs. separate `baseflow-explorer` repo. Subfolder keeps the demo in sync with the library and is simpler. Lean toward subfolder.
- **Data fetch library:** `dataretrieval` (official USGS client, friendlier) vs. raw NWIS JSON via `requests` (fewer deps, faster cold start). Lean toward `dataretrieval`.
- **Curated gage list:** which 10? Need a spread of hydrologic regimes to showcase where each method shines or fails. Draft list TBD.
- **Parameter UX:** defaults only, or expose method-specific sliders? Start with defaults, add sliders in phase 4.
- **Transitional shim lifespan:** how long to maintain `pybaseflow 0.2.x` on PyPI? Suggest one year, then stop publishing updates to it (the last version stays installable indefinitely).
