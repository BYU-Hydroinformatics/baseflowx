# pybaseflow Refactoring Plan

## Goal

Strip pybaseflow down to a clean, focused Python package for running baseflow separation algorithms on a hydrograph. Remove everything tied to the original paper's large-scale study (multi-station batch processing, geospatial ice lookups, KGE comparisons, bundled data files).

---

## Phase 1: Remove Study-Specific Code

### Files to delete entirely

| File | Reason |
|---|---|
| `baseflow/read_npz.py` | Only reads `thawed.npz`, a study artifact |
| `baseflow/thawed.npz` | Global thaw/freeze raster data used in the study |
| `baseflow/example.csv` | Study-specific example stations; replace later with a minimal example if needed |
| `docs/` (entire folder) | Will be rewritten from scratch |
| `mkdocs.yml` | Old docs config |
| `readthedocs.yml` | Old docs config |

### Functions/code to remove

**`separation.py`**
- `separation()` — batch multi-station separator with `tqdm` loop, `df_sta` geospatial lookups, KGE return logic. This is the study's batch runner, not a separation algorithm.
- `single()` — orchestrator that couples all methods together with KGE evaluation, ice handling. This is study comparison glue.
- The `if __name__ == "__main__"` block at the bottom — test scaffold for `bn77`.
- Duplicate `import numpy as np` and `from numba import njit, prange` statements mid-file (lines 429 and 796-797).

**`utils.py`**
- `geo2imagexy()` — converts geographic coords to raster pixel coords for the thawed.npz lookup. Study-specific.
- `exist_ice()` — ice period detection for the global study.
- `kge()` — evaluation metric for comparing methods in the study.
- `format_method()` — convenience for the `single()`/`separation()` orchestrators.
- `flow_duration_curve()` — utility that doesn't belong in a baseflow separation package (or move to a separate analysis/plotting package).

**`estimate.py`**
- `maxmium_BFI()` — estimates BFImax from annual data, used inside `single()` (currently commented out there). Remove unless Eckhardt filter needs a default BFImax estimator.

**`baseflow/plotting`** (note: this is a file, not a directory, and has no `.py` extension)
- Imports `plotly` and `baseflow` in broken ways. Plotting can be reintroduced properly later.

### Functions to keep and clean up

**`separation.py`**
- `strict_baseflow()` — useful for a variety of applications. Required by `recession_coefficient()`, which is used by 6 filter methods. Clean up: remove the `ice` parameter coupling to the thawed.npz data (already optional).
- `bn77()` and its helpers (`_estimate_recession_slope`, `_identify_recession_episodes`, `_eliminate_points`) — legitimate Brutsaert-Nieber recession analysis method with broad utility. Fix bug: `bn77()` passes 5 args to `_eliminate_points()` which expects 7 (missing `S` and `quantile` args). Rewrite to work correctly.

**`__init__.py`**
- Remove `example = Path(__file__).parent / 'example.csv'` line.
- Clean up wildcard imports after the module contents change.

**`setup.py`**
- Remove `package_data` entry for `example.csv` and `thawed.npz`.
- Update package name to `pybaseflow`, author, URL, description.
- Remove the `UploadCommand` class (use `twine` directly or switch to `pyproject.toml`).
- Remove `tqdm` from `REQUIRED` (no longer needed after removing batch `separation()`).

---

## Phase 2: Clean Up the Separation Methods

### What stays (the core algorithms)

These are all legitimate, published baseflow separation methods:

| Method | Function | Type |
|---|---|---|
| Lyne-Hollick | `lh()` | Recursive digital filter (1-pass, 2-pass) |
| Lyne-Hollick multi-pass | `lh_multi()` | Recursive digital filter (n-pass) |
| Chapman | `chapman()` | Recursive digital filter |
| Chapman-Maxwell | `chapman_maxwell()` | Recursive digital filter |
| Boughton | `boughton()` | Two-parameter recursive digital filter |
| Furey-Gupta | `furey()` | Physically-based digital filter |
| Eckhardt | `eckhardt()` | Two-parameter recursive digital filter |
| EWMA | `ewma()` | Exponential weighted moving average filter |
| Willems | `willems()` | Digital filter |
| UKIH | `ukih()` | Graphical (smoothed minima) |
| Local minimum | `local()` | Graphical (HYSEP) |
| Fixed interval | `fixed()` | Graphical (HYSEP) |
| Sliding interval | `slide()` | Graphical (HYSEP) |
| WHAT | `what()` | Eckhardt variant |
| HydRun | `hyd_run()` | Multi-pass digital filter |

### Cleanup tasks for kept methods

1. **Deduplicate the `initial_method` boilerplate** — the same 15-line block for setting `b[0]` is copy-pasted in 7+ functions. Extract a helper: `_init_baseflow(Q, initial_method)`.

2. **Deduplicate the `return_exceed` pattern** — similarly repeated. Consider whether this feature is needed at all for end-user usage. If kept, standardize via a wrapper or the init helper.

3. **Standardize function signatures** — all filter methods should accept `(Q, ...)` as a plain numpy array and return a numpy array of baseflow. Remove pandas coupling from individual methods.

4. **`what()` vs `eckhardt()`** — `what()` is mathematically identical to `eckhardt()`. Either remove `what()` and note it in the Eckhardt docstring, or keep it as an alias.

5. **`hyd_run()`** — appears to be a standalone multi-pass digital filter added separately. Review whether it duplicates `lh_multi()`. If distinct, keep and document; if redundant, remove.

6. **Helper functions** — `hysep_interval()`, `linear_interpolation()`, `local_turn()`, `ukih_turn()`, `slide_interpolation()`, `fixed_interpolation()` are internal helpers used by the graphical methods. Keep but prefix with `_` to signal they're private.

7. **`numba` dependency** — several functions use `@njit`. Consider whether the performance gain justifies the dependency for typical use cases (single hydrograph). Could offer numba as optional. Low priority.

---

## Phase 3: Clean Up Supporting Modules

### `estimate.py` — keep and simplify

- `recession_coefficient()` — needed by Chapman, CM, Boughton, Furey, Eckhardt, Willems. Keep.
- `param_calibrate()` / `param_calibrate_jit()` — auto-calibration for Boughton, Furey, Eckhardt, EWMA, Willems. Keep. These rely on `recession_period()`, `moving_average()`, `multi_arange()`, and `backward()` which should also stay.
- `maxmium_BFI()` — optional BFImax estimator for Eckhardt. Keep if useful as a standalone utility.

### `utils.py` — slim down

After removals, what remains:
- `clean_streamflow()` — useful input cleaning. Keep.
- `moving_average()` — used by `recession_period()`. Keep.
- `multi_arange()` — used by `recession_period()`. Keep.
- `backward()` — used by `maxmium_BFI()`. Keep.

---

## Phase 4: Modernize Packaging

1. Replace `setup.py` with `pyproject.toml`.
2. Rename package to `pybaseflow` (both the directory and the pip install name).
3. Pin minimum Python version to 3.9+.
4. Add a proper `__version__` in `__init__.py`.
5. Add a `py.typed` marker if type hints are added.

---

## Phase 5: Add a Simple User-Facing API

Create a clean top-level function that replaces the removed `single()`/`separation()`:

```python
import pybaseflow

# Run a single method
baseflow = pybaseflow.eckhardt(Q, a=0.98, BFImax=0.8)

# Or use individual functions directly
from pybaseflow.separation import lh, chapman, eckhardt
```

The individual method functions should be the primary API. A convenience wrapper could be added later but is not essential.

---

## Phase 6: Baseflow Separation Methods to Add

The following are well-established methods not currently in the package:

### High priority

| Method | Reference | Notes |
|---|---|---|
| **Nathan-McMahon (1990)** | Nathan & McMahon, 1990 | One of the most widely used single-parameter recursive digital filters. Often confused with Lyne-Hollick but has its own parameterization. Verify `lh()` doesn't already cover this. |
| **Sweetwater (SWT) / RDF** | Eckhardt, 2008 | Generalized recursive digital filter framework. Eckhardt showed most filters are special cases of a general form. Implementing the general form would unify several existing methods. |
| **PART** | Rutledge, 1998 (USGS) | USGS program for estimating baseflow from daily streamflow using the antecedent recession requirement. Very widely used in practice, especially in the US. |
| **BFlow** | Arnold & Allen, 1999 | Automated baseflow separation, widely used with SWAT. Based on recursive digital filter with specific calibration approach. |
| **Conductivity Mass Balance (CMB)** | Stewart et al., 2007 | Uses specific conductance as a tracer to separate baseflow. Different paradigm (tracer-based) but very useful when conductivity data is available. |

### Medium priority

| Method | Reference | Notes |
|---|---|---|
| **Brutsaert-Nieber (1977)** | Brutsaert & Nieber, 1977 | Already in codebase as `bn77()` but has bugs — fix rather than rewrite. See Phase 1 notes. |
| **HYSEP Sliding (improved)** | Sloto & Crouse, 1996 | The current `slide()` implementation may not perfectly match the USGS HYSEP specification. Worth verifying and correcting. |
| **Constant-k** | Blume et al., 2007 | Simple exponential recession filter with a constant recession constant. Easy to implement, commonly used as a baseline. |
| **Jakeman-Hornberger (IHACRES)** | Jakeman & Hornberger, 1993 | Two-component recursive filter from the IHACRES rainfall-runoff model. Separates quick flow and slow flow. |

### Lower priority / future

| Method | Reference | Notes |
|---|---|---|
| **Generalized Conductivity Mass Balance** | Miller et al., 2014 | Extension of CMB that accounts for multiple end-members. |
| **Wavelet-based separation** | Various | Uses wavelet decomposition to separate high-frequency (quick flow) from low-frequency (baseflow) signals. Research-grade, less standardized. |
| **Kalinlin-type** | Various Russian literature | Recession-based methods popular in Eastern European hydrology. Niche but would broaden the package's scope. |

---

## Suggested Order of Work

1. **Delete** study-specific files and code (Phase 1)
2. **Clean up** remaining methods and helpers (Phase 2-3)
3. **Modernize** packaging (Phase 4)
4. **Add** PART and BFlow as the first new methods (Phase 6, high priority)
5. **Rewrite** documentation
6. **Add** remaining new methods incrementally
