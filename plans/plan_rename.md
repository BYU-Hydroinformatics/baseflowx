# Rename `pybaseflow` → `baseflowx`

PyPI does **not** allow renaming an existing project — you publish under the new name and retire the old one. Steps:

## GitHub repo

1. Decide the new repo slug: `baseflowx` (matches package, clean URL).
2. On GitHub: **Settings → General → Rename repository** to `baseflowx`. GitHub sets up automatic redirects from `njones61/pybaseflow` → `njones61/baseflowx` for clones, issues, and PRs — safe.
3. Update local remote: `git remote set-url origin git@github.com:njones61/baseflowx.git`.
4. Rename the working directory locally: `/Users/njones/python_projects/pybaseflow` → `/Users/njones/python_projects/baseflowx`. (Optional but tidy.)

## Python package

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

## PyPI

1. Register and publish `baseflowx 0.2.0` (fresh project under the new name — new PyPI project, you're the sole owner).
2. On the old `pybaseflow` PyPI project, publish **one** final version `pybaseflow 0.2.0` as a thin redirect shim:
   - Description prominently says "Renamed to baseflowx. Install `baseflowx` instead."
   - `pyproject.toml` lists `baseflowx>=0.2.0` as its only runtime dependency.
   - `pybaseflow/__init__.py` does `from baseflowx import *` and emits a `DeprecationWarning` on import.
   - This keeps every existing `pip install pybaseflow` working while nudging users off it.
3. Do **not** yank older `pybaseflow` releases — breaks anyone pinning a specific version. Leave them up with the deprecation notice on the project page.

## Docs / site

1. Update `mkdocs.yml` site name, repo URL, and any hard-coded install commands.
2. Docs are served on ReadTheDocs. The RTD project is tied to a specific repo and slug — after the GitHub rename, check that the RTD project's repo URL still resolves (GitHub redirects help) and decide whether to rename the RTD project from `pybaseflow` to `baseflowx` so the served URL matches `mkdocs.yml`'s `site_url`.
3. Update README badges (PyPI, CI, docs link) to point to the new project.

## Verification checklist

- [ ] `pip install baseflowx` works in a fresh venv; `import baseflowx` succeeds.
- [ ] `pip install pybaseflow` still works and imports without error (deprecation warning visible).
- [ ] MkDocs site builds and deploys under the new URL.
- [ ] `git clone` of the old repo URL still succeeds (via GitHub redirect).
- [ ] All tests pass.

## Open questions

- **Transitional shim lifespan:** how long to maintain `pybaseflow 0.2.x` on PyPI? Suggest one year, then stop publishing updates to it (the last version stays installable indefinitely).
