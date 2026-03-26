# pybaseflow Refactoring Plan

## Goal

Strip pybaseflow down to a clean, focused Python package for running baseflow separation algorithms on a hydrograph. Remove everything tied to the original paper's large-scale study (multi-station batch processing, geospatial ice lookups, KGE comparisons, bundled data files). The refactored package will serve as the basis for a journal article describing pybaseflow as a comprehensive, well-documented Python toolkit for baseflow separation.

---

## Phase 1: Remove Study-Specific Code

### Files to delete entirely

| File | Reason |
|---|---|
| `baseflow/read_npz.py` | Only reads `thawed.npz`, a study artifact |
| `baseflow/thawed.npz` | Global thaw/freeze raster data used in the study |
| `baseflow/example.csv` | Study-specific example stations; replaced with new sample data in Phase 4 |
| `docs/` (entire folder) | Will be rebuilt from scratch in Phase 5 |
| `mkdocs.yml` | Old docs config; replaced in Phase 5 |
| `readthedocs.yml` | Old docs config; replaced in Phase 5 |

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

## Phase 2: Clean Up Separation Methods and Supporting Modules

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

## Phase 3: Modernize Packaging

1. Replace `setup.py` with `pyproject.toml`.
2. Rename package to `pybaseflow` (both the directory and the pip install name).
3. Pin minimum Python version to 3.9+.
4. Add a proper `__version__` in `__init__.py`.
5. Add a `py.typed` marker if type hints are added.
6. **Bundle sample data** in `pybaseflow/data/`:
   - A small daily streamflow CSV (1-3 years) from a USGS gauge with public-domain data (e.g., NWIS). Include date, discharge, and ideally specific conductance for CMB examples.
   - A convenience loader: `pybaseflow.load_sample_data()` returning a numpy array (or dict with Q, dates, and optionally SC).
   - Keep the data file small (<500 KB) and include it in `package_data`.
   - Document the gauge ID, period, and source in a `data/README.md`.
7. **USGS data retrieval** — add `pybaseflow.fetch_usgs(site_id, start_date, end_date, parameter='discharge')` that pulls daily values from the USGS NWIS Water Services API:
   - Use the NWIS daily values REST endpoint (`https://waterservices.usgs.gov/nwis/dv/`)
   - Support parameter codes: `'00060'` (discharge) and `'00095'` (specific conductance)
   - Accept a convenience `parameter` kwarg mapping `'discharge'` → `'00060'`, `'sc'` → `'00095'`
   - Return a dict with `dates` (array of datetime), `Q` (numpy array in cms or cfs with a `units` field), and optionally `SC`
   - Use only stdlib (`urllib`) to avoid adding a requests dependency; or make it an optional dependency
   - Handle common issues: missing days, provisional data flags, ice-affected flags, zero-flow padding
   - Place in `pybaseflow/io.py` as a standalone module
   - Example usage:
     ```python
     from pybaseflow.io import fetch_usgs
     data = fetch_usgs('01013500', '2015-01-01', '2020-12-31')
     b = pybaseflow.eckhardt(data['Q'], a=0.98, BFImax=0.8)
     ```
8. **User-facing API** — individual method functions are the primary API, accessible from the top-level namespace via wildcard imports:
   ```python
   import pybaseflow
   baseflow = pybaseflow.eckhardt(Q, a=0.98, BFImax=0.8)

   # Or import directly
   from pybaseflow.separation import lh, chapman, eckhardt
   ```

---

## Phase 4: Add New Methods

Detailed research into each candidate method, with findings that inform implementation priority and approach.

### Existing Filter Taxonomy

Before adding new methods, it is important to understand the structural taxonomy of the existing recursive digital filters. Eckhardt (2005, 2008) showed that most filters follow the general form:

```
b(t) = alpha * b(t-1) + beta * (Q(t) + gamma * Q(t-1))
```

The filters split into **two structural families** based on gamma:

**Family 1: gamma = 0 (current-timestep-only, linear reservoir derived)**

| Filter | alpha | beta | gamma |
|---|---|---|---|
| Chapman-Maxwell (1996) | `k / (2-k)` | `(1-k) / (2-k)` | 0 |
| Boughton (1993) | `k / (1+C)` | `C / (1+C)` | 0 |
| Eckhardt (2005) | `(1-BFImax)*k / (1-k*BFImax)` | `(1-k)*BFImax / (1-k*BFImax)` | 0 |
| EWMA / Tularam-Ilahee (2008) | `a` | `1-a` | 0 |

Key relationships within this family:
- Chapman-Maxwell is Eckhardt with `BFImax = 0.5`
- Boughton is Eckhardt with `C = (1-k)*BFImax / (1-BFImax)` (a reparameterization)
- EWMA corresponds to the degenerate limit `BFImax = 1`

**Family 2: gamma = 1 (two-timestep, signal-processing derived)**

| Filter | alpha | beta | gamma |
|---|---|---|---|
| Lyne-Hollick (1979) | `a` | `(1-a)/2` | 1 |
| Chapman (1991) | `(3k-1)/(3-k)` | `(1-k)/(3-k)` | 1 |
| Willems (2009) | `(a-v)/(1+v)` | `v/(1+v)` | 1 |

Lyne-Hollick, Chapman (1991), and Willems are **not** special cases of the Eckhardt equation — they include `Q(t-1)`, while Eckhardt's form uses only `Q(t)`. The gamma=0 family derives from linear reservoir theory (physically based); the gamma=1 family originates from signal processing (low-pass filtering).

This taxonomy directly informs the Phase 2 refactoring: implement a private `_recursive_digital_filter(Q, alpha, beta, gamma, b0)` as the shared computational core, then keep all named public functions as thin wrappers that compute their specific alpha/beta/gamma from published parameters. This eliminates the boilerplate duplication while preserving the user-facing API and makes adding new filters trivial.

References:
- Eckhardt, K. (2005). How to construct recursive digital filters for baseflow separation. *Hydrological Processes*, 19, 507-515.
- Eckhardt, K. (2008). A comparison of baseflow indices. *Journal of Hydrology*, 352, 168-173.
- Eckhardt, K. (2023). How physically based is hydrograph separation by recursive digital filtering? *HESS*, 27, 495-499.

---

### 4.1 Nathan-McMahon (1990) — ALREADY IMPLEMENTED

**Reference:** Nathan, R.J. and McMahon, T.A. (1990). Evaluation of automated techniques for base flow and recession analyses. *Water Resources Research*, 26(7), 1465-1473.

**Finding:** Nathan-McMahon is **not a distinct algorithm**. It is the Lyne-Hollick (1979) filter with a specific recommended parameterization and application protocol. The equation is identical:

```
b(t) = alpha * b(t-1) + (1-alpha)/2 * (Q(t) + Q(t-1))
```

This matches the existing `lh()` implementation exactly (separation.py line 535). Nathan & McMahon's contributions were:
1. Standardized `alpha = 0.925` (the existing default in `lh()`)
2. Validated the filter against manual separation across 186 Australian catchments
3. Recommended a 3-pass protocol (forward-backward-forward), whereas the current `lh()` does 2 passes

**Implementation:** No new function needed. Document the Nathan-McMahon (1990) reference in the `lh()` docstring. Note that `lh_multi(Q, beta=0.925, num_pass=3)` reproduces the Nathan-McMahon recommended protocol exactly.

---

### 4.2 Generalized Recursive Digital Filter Refactoring

**Reference:** Eckhardt, K. (2005, 2008) — see taxonomy section above.

**Implementation:** This is an internal refactoring task, not a new user-facing method. The approach:

1. Create `_recursive_digital_filter(Q, alpha, beta, gamma, b0)` as the shared computational core
2. All existing filter functions become thin wrappers that:
   - Accept their published parameters (e.g., `eckhardt(Q, a, BFImax)`)
   - Compute alpha/beta/gamma from those parameters
   - Delegate to the core function
3. This also makes adding new filters (e.g., IHACRES) trivial — just define the parameter mapping

The core function handles:
- Initial value setting (`b[0] = b0`)
- Forward iteration with the general equation
- Streamflow capping (`b[t] = min(b[t], Q[t])`)
- The `return_exceed` pattern (if retained)

This should be done as part of Phase 2 cleanup, not Phase 6.

---

### 4.3 PART (Rutledge, 1998) — NEW, HIGH PRIORITY

**Reference:** Rutledge, A.T. (1998). Computer programs for describing the recession of ground-water discharge and for estimating mean ground-water recharge and discharge from streamflow records — Update. *USGS Water-Resources Investigations Report 98-4148*.

**Why this matters:** PART is one of the most widely used baseflow separation methods in US groundwater studies. It represents a fundamentally different paradigm from digital filters — a recession-based graphical method that identifies days where streamflow equals baseflow, then interpolates. No Python implementation currently exists. This would be a significant and unique addition to the package.

**Algorithm (from Figure 12, USGS WRI 98-4148):**

**Step 1 — Compute the recession duration N:**
```
N = A^0.2
```
where A = drainage area in square miles. The existing `hysep_interval()` already computes this (with km²-to-sq-mi conversion). N is generally not an integer.

**Step 2 — Identify qualifying (groundwater discharge) days:**
A day qualifies if streamflow has been continuously declining (or equal) for the preceding N consecutive days: `Q[i] >= Q[i+1]` for all i in the N-day window before the day in question.

**Step 3 — Apply the 0.1 log-cycle check:**
For each qualifying day, check if it is followed by a daily decline exceeding 0.1 log cycles:
```
if log10(Q[day]) - log10(Q[day+1]) > 0.1  →  disqualify the day
```
Rationale (Barnes, 1939): A decline > 0.1 log cycle/day indicates surface runoff is still present. The 0.1 threshold is the default but can be user-adjustable.

**Step 4 — Set baseflow on qualifying days:**
On all remaining qualifying days, set `baseflow = streamflow`.

**Step 5 — Log-linear interpolation:**
For all non-qualifying days, interpolate baseflow linearly between the log values of the nearest preceding and following qualifying days:
```
log10(b[t]) = linear_interp(log10(b[t_prev]), log10(b[t_next]))
b[t] = 10^(interpolated value)
```
This produces exponential-decay curves between anchor points.

**Step 6 — Iterative correction for baseflow > streamflow:**
If any interpolated baseflow exceeds streamflow:
1. Find the day in the offending interval with the largest `log(baseflow) - log(streamflow)`
2. Make that day a new anchor point (baseflow = streamflow)
3. Return to Step 5 and re-interpolate
4. Repeat until no violations remain

**Step 7 — Three-pass combination:**
Run the entire procedure three times with:
- Pass 1: N = floor(A^0.2) (minimum 1)
- Pass 2: N = ceil(A^0.2) (minimum 2)
- Pass 3: N = ceil(A^0.2) + 1

Combine the three mean baseflow estimates using second-order polynomial (curvilinear) interpolation as a function of N, evaluated at the exact real-valued N = A^0.2. (Note: the USGS-R DVstats implementation simplifies this to linear interpolation between floor and ceil results.)

**Key parameters:**

| Parameter | Default | Notes |
|---|---|---|
| `area` | Required | Drainage area (km² or sq mi, with conversion) |
| `log_cycle_threshold` | 0.1 | Daily log10 decline threshold; user-adjustable |

**Edge cases:**
- Zero flows: convert to a small epsilon (e.g., 1e-99) to avoid log(0); set baseflow < 1e-6 to 0
- Record edges: extrapolate using the nearest qualifying day's value (constant extrapolation)
- Drainage area limits: recommended 1-500 sq mi; below 1, N < 1 day; above ~500, nonuniform weather

**Implementation notes:**
- Reuse `hysep_interval()` for N computation (already converts km² to sq mi)
- This is a graphical/recession method — no recursive filter, no numba needed
- The iterative correction loop (Step 6) converges quickly in practice
- Function signature: `part(Q, area, log_cycle_threshold=0.1)`
- Returns a numpy array of daily baseflow values

---

### 4.4 BFlow (Arnold & Allen, 1999) — FILTER EXISTS, ADD RECESSION ANALYSIS

**References:**
- Arnold, J.G. and Allen, P.M. (1999). Automated methods for estimating baseflow and ground water recharge from streamflow records. *JAWRA*, 35(2), 411-424.
- Arnold, J.G., Allen, P.M., Muttiah, R., and Bernhardt, G. (1995). Automated base flow separation and recession analysis techniques. *Groundwater*, 33(6), 1010-1018.

**Finding:** The BFlow filter equation is identical to Lyne-Hollick / Nathan-McMahon:
```
q(t) = beta * q(t-1) + ((1+beta)/2) * (Q(t) - Q(t-1))    [quickflow form]
b(t) = Q(t) - q(t)                                          [baseflow by subtraction]
```
with `beta = 0.925` and 3 passes (forward-backward-forward). This is already covered by `lh_multi(Q, beta=0.925, num_pass=3)`.

**What BFlow adds beyond the filter:** An automated recession analysis that computes the **alpha factor** (baseflow recession constant), which maps directly to SWAT's `ALPHA_BF` parameter. This is the novel contribution:

1. Identify baseflow recession segments from the separated baseflow time series
2. Select segments between NDMIN (default 10) and NDMAX (default 300) days
3. Assemble into a master recession curve
4. Fit the exponential recession: `Q(t) = Q(0) * exp(-alpha * t)`
5. Compute **baseflow days** (BFD): days for recession to decline one log cycle, where `alpha = 2.3 / BFD`

**Outputs per pass:**
- Baseflow fraction (Fr1, Fr2, Fr3 from each pass)
- Alpha factor (recession constant)
- Number of recession segments used (NPR)
- Baseflow days

**Implementation:**
- The filter itself requires no new code — use `lh_multi()` with `num_pass=3`
- Implement `bflow_recession_analysis(Q, b, ndmin=10, ndmax=300)` as a new function in `estimate.py` that takes the separated baseflow and returns the alpha factor, baseflow days, and number of recession segments
- Provide a convenience wrapper `bflow(Q, beta=0.925)` that runs both the filter and recession analysis, returning a results object with all BFlow outputs
- This gives SWAT users a direct path to `ALPHA_BF` calibration

---

### 4.5 Conductivity Mass Balance — CMB (Stewart et al., 2007) — NEW, HIGH PRIORITY

**Reference:** Stewart, M.T., Cimino, J., and Ross, M. (2007). Calibration of base flow separation methods with streamflow conductivity. *Groundwater*, 45(1), 17-27.

**Why this matters:** CMB is a completely different paradigm — tracer-based rather than mathematical. It provides a physically grounded baseflow estimate that is commonly used as a calibration reference for digital filters (e.g., calibrating Eckhardt's BFImax). No Python implementation currently exists. Adding CMB would make pybaseflow the first Python package to offer both filter-based and tracer-based separation methods.

**Core equation (two-component mixing model):**

Water balance: `Q(t) = Q_b(t) + Q_r(t)`
Solute mass balance: `Q(t) * SC(t) = Q_b(t) * SC_BF + Q_r(t) * SC_RO`

Solving for baseflow:
```
Q_b(t) = Q(t) * (SC(t) - SC_RO) / (SC_BF - SC_RO)
```

Where:
- `Q(t)` = total streamflow at time t
- `SC(t)` = measured specific conductance (uS/cm) at time t
- `SC_BF` = baseflow end-member conductivity (constant)
- `SC_RO` = surface runoff end-member conductivity (constant)

**Data requirements:**
- Streamflow time series `Q(t)` (daily or sub-daily)
- Concurrent specific conductance time series `SC(t)`
- Minimum 6 months of paired data recommended; multi-year records ideal

**End-member estimation (the critical step):**
- `SC_RO`: Use the **1st percentile** of the SC record (not the absolute minimum, which is noise-sensitive)
- `SC_BF`: Three approaches:
  1. Maximum SC — simple but sensitive to outliers (not recommended)
  2. 99th percentile of entire SC record — more robust
  3. Dynamic yearly 99th percentile with interpolation (recommended) — accommodates temporal changes in groundwater chemistry

Sensitivity: BFI is more sensitive to SC_BF than SC_RO. A 10% error in SC_BF causes ~13% error in BFI; the same error in SC_RO causes only ~6% error.

**Physical constraints:**
- If `Q_b(t) > Q(t)` (SC(t) > SC_BF): cap at `Q_b(t) = Q(t)`
- If `Q_b(t) < 0` (SC(t) < SC_RO): set `Q_b(t) = 0`

**When CMB works well:**
- Strong inverse correlation between Q and SC (r <= -0.5)
- Headwater/tributary streams, steep terrain
- Large contrast between groundwater and runoff conductivity
- Minimal anthropogenic influence

**When CMB fails:**
- No consistent inverse Q-SC correlation (r > -0.3)
- Anthropogenic interference (irrigation, mining, road salt)
- Reservoir influence (evaporation increases SC)
- Large catchments (>34,000 km² — only ~11% show strong inverse correlation)
- Non-conservative SC behavior

**Common integrated workflow in the literature:**
1. Collect 6+ months of paired Q and SC data
2. Run CMB to get a reference BFI
3. Use that BFI to calibrate Eckhardt's BFImax parameter
4. Apply calibrated Eckhardt to the full discharge record where SC data is unavailable

Note: Zhang et al. (2021) found that while cumulative baseflow from calibrated Eckhardt matches CMB well (NSE 0.91-1.00), daily baseflow series are inconsistent (mean NSE -0.30), suggesting the two methods measure fundamentally different things at short timescales.

**Implementation:**
- Function signature: `cmb(Q, SC, SC_BF=None, SC_RO=None, sc_bf_percentile=99, sc_ro_percentile=1)`
- If end-members are not provided, estimate from the SC record using percentiles
- Add `estimate_endmembers(SC, bf_percentile=99, ro_percentile=1)` utility
- Returns numpy array of daily baseflow values
- Add `calibrate_eckhardt_from_cmb(Q, SC)` convenience function that runs CMB, estimates BFI, and returns a calibrated BFImax for use with `eckhardt()`
- Place in a new module `pybaseflow/tracer.py` to separate tracer-based methods from digital filters

**Additional references:**
- Li, Q., et al. (2020). Key challenges facing the application of the CMB method. *HESS*, 24, 6075-6090.
- Zhang, J., et al. (2021). Can the two-parameter recursive digital filter really be calibrated by CMB? *HESS*, 25, 1747-1760.
- Cartwright, I., et al. (2022). Implications of variations in stream specific conductivity. *HESS*, 26, 183-195.
- Mei, Y., et al. (2024). Optimal baseflow separation through chemical mass balance. *Water Resources Research*, 60, e2023WR036386.

---

### 4.6 Jakeman-Hornberger / IHACRES (1993) — NEW, MEDIUM PRIORITY

**Reference:** Jakeman, A.J. and Hornberger, G.M. (1993). How much complexity is warranted in a rainfall-runoff model? *Water Resources Research*, 29(8), 2637-2649.

**Finding:** The IHACRES baseflow filter is a 3-parameter extension of the existing Boughton filter. The equation:

```
b(t) = [a / (1+C)] * b(t-1) + [C / (1+C)] * (Q(t) + alpha_s * Q(t-1))
```

When `alpha_s = 0`, this reduces exactly to `boughton()` (separation.py line 182). The extra `alpha_s * Q(t-1)` term adds dependence on the previous timestep's total flow, accounting for delayed baseflow response.

**Relationship to the full IHACRES model:**
The IHACRES rainfall-runoff model routes effective rainfall through two parallel linear stores:
```
X_q(t) = alpha_q * X_q(t-1) + beta_q * U(t)    [quick flow store]
X_s(t) = alpha_s * X_s(t-1) + beta_s * U(t)    [slow flow store]
Q(t) = X_q(t) + X_s(t)
```
The baseflow filter form above is the inverse problem: extracting the slow flow component from observed total streamflow without requiring rainfall data.

**Key parameters:**

| Parameter | Meaning | Typical Range | Constraint |
|---|---|---|---|
| `a` | Recession coefficient for slow flow | 0.9-0.999 | Related to `tau_s`: `a = exp(-1/tau_s)` |
| `C` | Partitioning parameter | 0.1-1.0 | Controls fraction attributed to baseflow |
| `alpha_s` | Previous-timestep dependence | -0.99 to 0 | `alpha_s = -exp(-1/k)` where k is recession constant; 0 reduces to Boughton |

**In the generalized filter taxonomy:**
```
alpha = a / (1 + C)
beta  = C / (1 + C)
gamma = alpha_s    (NOT 0 or 1, but a calibrated value)
```
This makes IHACRES unique — it is neither purely gamma=0 nor gamma=1, but occupies the continuum between the two families.

**Implementation:**
- Extend the Boughton implementation pattern with one additional parameter
- Function signature: `ihacres(Q, a, C, alpha_s, initial_method='Q0', return_exceed=False)`
- Calibration: extend `param_calibrate()` to handle the 3-parameter space, or provide alpha_s estimation from recession analysis
- With the generalized `_recursive_digital_filter()` core, this is trivial to add

---

### 4.7 Brutsaert-Nieber (1977) — BUG FIX, MEDIUM PRIORITY

**Reference:** Brutsaert, W. and Nieber, J.W. (1977). Regionalized drought flow hydrographs from a mature glaciated plateau. *Water Resources Research*, 13(3), 637-643.

**Finding:** The existing `bn77()` implementation has a confirmed critical bug.

**The bug (separation.py line 821):**
```python
drought_flow_points = _eliminate_points(recession_episodes, L_min,
    snow_freeze_period, observational_precision, Q, quantile)
```

**The function signature (line 879):**
```python
def _eliminate_points(recession_episodes, L_min, snow_freeze_period,
    observational_precision, Q, S, quantile):
```

`_eliminate_points()` expects 7 arguments but is called with 6 — the recession slope array `S` is missing. `S` is computed at line 815 via `S = _estimate_recession_slope(Q)` but never passed to the function. This causes a runtime crash (Python will bind `quantile` to the `S` parameter and raise TypeError for missing `quantile`).

**Fix:** Add `S` to the call:
```python
drought_flow_points = _eliminate_points(recession_episodes, L_min,
    snow_freeze_period, observational_precision, Q, S, quantile)
```

After fixing the bug, review the full bn77 pipeline for correctness:
- `_estimate_recession_slope(Q)` — computes dQ/dt
- `_identify_recession_episodes(Q)` — finds continuous recession periods
- `_eliminate_points()` — filters out non-drought points using S and quantile thresholds

---

### 4.8 HYSEP Sliding Verification — LOW PRIORITY

**Reference:** Sloto, R.A. and Crouse, M.Y. (1996). HYSEP: A computer program for streamflow hydrograph separation and analysis. *USGS Water-Resources Investigations Report 96-4040*.

**Finding:** The current `slide()` implementation (separation.py line 592) uses a sliding window minimum approach with the correct interval calculation (`N = A^0.2`, odd integer between 3-11). Edge handling uses the minimum of the edge segment, which is reasonable but may differ slightly from the USGS spec in some interpretations. This is a minor verification task, not a rewrite.

---

### 4.9 Constant-k (Blume et al., 2007) — LOW PRIORITY, LIKELY SKIP

**Reference:** Blume, T., Zehe, E., and Bronstert, A. (2007). Rainfall-runoff response, event-based runoff coefficients and hydrograph separation. *Hydrological Sciences Journal*, 52(5), 843-862.

**Finding:** The constant-k method is an **event-based** recession technique, not a continuous filter. During recession: `b(t) = k * b(t-1)`. This is what every linear-reservoir-based filter reduces to when `b(t) = Q(t)` (no direct runoff). The recession constant `k` is already computed by `recession_coefficient()` in estimate.py.

Implementing constant-k as a standalone method would require:
- Storm event identification as a preprocessing step
- Extrapolation of the pre-event recession beneath the event hydrograph
- Event boundary detection using the instantaneous recession ratio

The event-detection complexity adds significant scope for marginal value — the recession behavior is already implicit in all the existing filters. **Recommend skipping** unless there is a specific use case for event-based separation.

---

### 4.10 Lower Priority / Future

| Method | Reference | Notes |
|---|---|---|
| **Generalized Conductivity Mass Balance** | Miller et al., 2014 | Extension of CMB for multiple end-members. Implement after CMB is stable. |
| **Wavelet-based separation** | Various | Wavelet decomposition to separate high/low frequency signals. Research-grade, less standardized. |
| **Kalinlin-type** | Various Russian literature | Recession-based methods popular in Eastern European hydrology. Niche but broadens scope. |

---

## Summary: Method Coverage After Refactoring

The refactored pybaseflow package will offer methods spanning four distinct paradigms:

### Recursive Digital Filters (gamma=0 family — linear reservoir based)
| Method | Function | Parameters | Reference |
|---|---|---|---|
| Chapman-Maxwell | `chapman_maxwell()` | k | Chapman & Maxwell, 1996 |
| Boughton | `boughton()` | k, C | Boughton, 1993 |
| Furey-Gupta | `furey()` | a, A | Furey & Gupta, 2001 |
| Eckhardt | `eckhardt()` | a, BFImax | Eckhardt, 2005 |
| EWMA | `ewma()` | e | Tularam & Ilahee, 2008 |
| WHAT | `what()` (alias) | a, BFImax | Lim et al., 2005 |

### Recursive Digital Filters (gamma=1 family — signal processing based)
| Method | Function | Parameters | Reference |
|---|---|---|---|
| Lyne-Hollick | `lh()` / `lh_multi()` | beta, num_pass | Lyne & Hollick, 1979; Nathan & McMahon, 1990 |
| Chapman | `chapman()` | k | Chapman, 1991 |
| Willems | `willems()` | a, w | Willems, 2009 |

### Recursive Digital Filters (variable gamma — hybrid)
| Method | Function | Parameters | Reference |
|---|---|---|---|
| Jakeman-Hornberger | `ihacres()` | a, C, alpha_s | Jakeman & Hornberger, 1993 |

### Graphical / Recession-Based Methods
| Method | Function | Parameters | Reference |
|---|---|---|---|
| UKIH | `ukih()` | area | Piggott et al., 2005 |
| Local minimum | `local()` | area | Sloto & Crouse, 1996 |
| Fixed interval | `fixed()` | area | Sloto & Crouse, 1996 |
| Sliding interval | `slide()` | area | Sloto & Crouse, 1996 |
| PART | `part()` | area | Rutledge, 1998 |
| Brutsaert-Nieber | `bn77()` | L_min, quantile | Brutsaert & Nieber, 1977 |

### Tracer-Based Methods
| Method | Function | Parameters | Reference |
|---|---|---|---|
| Conductivity Mass Balance | `cmb()` | SC, SC_BF, SC_RO | Stewart et al., 2007 |

### Supporting Functions
| Function | Module | Purpose |
|---|---|---|
| `recession_coefficient()` | estimate.py | Compute recession constant from discharge |
| `param_calibrate()` | estimate.py | Auto-calibrate filter parameters |
| `maxmium_BFI()` | estimate.py | Estimate BFImax from annual data |
| `bflow_recession_analysis()` | estimate.py | BFlow/SWAT alpha factor computation |
| `estimate_endmembers()` | tracer.py | Estimate SC_BF and SC_RO from conductivity data |
| `calibrate_eckhardt_from_cmb()` | tracer.py | Calibrate Eckhardt BFImax using CMB reference |
| `clean_streamflow()` | utils.py | Input validation and cleaning |

---

## Progress and Order of Work

1. **Phase 1 — Delete study-specific code** ✓
2. **Phase 2 — Refactor existing methods** ✓
   - Implemented `_recursive_digital_filter()` core (generalized RDF)
   - Converted existing filters to thin wrappers
   - Extracted `_init_baseflow()` to deduplicate boilerplate
   - Fixed `bn77()` bug (missing `S` argument)
   - Removed `hyd_run()` (duplicate of `lh_multi()`)
   - Made `what()` an alias for `eckhardt()`
   - Prefixed internal helpers with `_`
   - Cleaned up estimate.py and utils.py
3. **Phase 3 — Modernize packaging** ✓
   - Replaced `setup.py` with `pyproject.toml`
   - Renamed package to `pybaseflow`
   - Added `__version__`, sample data, `fetch_usgs()`, `load_sample_data()`
4. **Phase 4 — Add new methods** ✓
   - **PART** — Rutledge 1998 recession-based graphical method
   - **CMB** + `calibrate_eckhardt_from_cmb()` — tracer-based separation in `tracer.py`
   - **BFlow** recession analysis — `bflow()` and `bflow_recession_analysis()` for SWAT interop
   - **IHACRES** — 3-parameter Boughton extension via `_recursive_digital_filter()`
5. **Phase 5 — Rebuild documentation from scratch**
   a. **Docs infrastructure** — set up MkDocs with Material theme, `mkdocstrings` for auto-generated API reference, `pyproject.toml` integration, and GitHub Pages deployment
   b. **API reference** — auto-generate from docstrings; ensure every public function has a complete docstring with parameters, return type, equation, and literature reference
   c. **Method guide** — a narrative page for each method category (digital filters, graphical/recession, tracer-based) explaining the theory, when to use each method, and how they relate to each other; include the two-family filter taxonomy (gamma=0 vs gamma=1)
   d. **Getting started** — installation, minimal example (load a hydrograph, run a filter, plot the result), and a quick comparison of 3-4 methods on the same data
   e. **Worked examples** — Jupyter notebooks or literate docs showing:
      - Single-method usage for each method
      - Multi-method comparison on a real hydrograph
      - Parameter calibration workflow (recession coefficient → filter parameters)
      - CMB-to-Eckhardt calibration pipeline
      - BFlow recession analysis for SWAT users
   f. **Parameter selection guide** — practical guidance on choosing filter parameters, recession coefficient estimation, BFImax estimation, and when default values are appropriate
   g. **Method comparison table** — a single reference page listing every method with its equation, parameters, type, reference, and recommended use case
   h. **Contributing guide** — how to add a new filter method using the `_recursive_digital_filter()` core
   i. **Unified notation** — consistent variable names across all docs (Q for streamflow, b for baseflow, a/k for recession coefficient, etc.)
6. **Phase 6 — Journal article**
   - Describing pybaseflow's scope, taxonomy, and implementation
   - Method comparison / benchmarking suite with real hydrograph data
   - Positioning relative to existing tools (USGS HYSEP, PART, BFlow, R packages)
