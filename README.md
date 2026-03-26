# pybaseflow

A comprehensive Python toolkit for baseflow separation from streamflow hydrographs.

pybaseflow implements 17 baseflow separation methods spanning four paradigms: recursive digital filters, graphical/interval methods, recession-based methods, and tracer-based methods. It also provides parameter estimation, USGS data retrieval, and a CMB-to-Eckhardt calibration bridge.

This project builds on the [baseflow](https://github.com/xiejx5/baseflow) package by Xie et al. (2020), which implemented methods described in "Evaluation of typical methods for baseflow separation in the contiguous United States" (Journal of Hydrology, 583, 124628). pybaseflow extends that work with new methods (PART, CMB, BFlow, IHACRES), a unified filter architecture, modern packaging, and USGS data integration.

This project is funded by [CIROH](https://ciroh.ua.edu/).

## Install

```bash
pip install pybaseflow
```

## Quick Start

```python
import pybaseflow

# Load bundled sample data (USGS 01013500, Fish River, Maine)
data = pybaseflow.load_sample_data()
Q = data['Q']

# Run a single method
b = pybaseflow.eckhardt(Q, a=0.98, BFImax=0.8)

# Or fetch your own data from USGS NWIS
from pybaseflow.io import fetch_usgs
data = fetch_usgs('01013500', '2019-01-01', '2020-12-31')
```

## Methods

### Recursive Digital Filters

All recursive digital filters share a generalized core:
`b[t] = alpha * b[t-1] + beta * (Q[t] + gamma * Q[t-1])`

**gamma=0 family** (linear reservoir based):

| Method | Function | Reference |
|---|---|---|
| Boughton | `boughton(Q, a, C)` | Boughton, 1993 |
| Chapman-Maxwell | `chapman_maxwell(Q, a)` | Chapman & Maxwell, 1996 |
| Eckhardt | `eckhardt(Q, a, BFImax)` | Eckhardt, 2005 |
| EWMA | `ewma(Q, e)` | Tularam & Ilahee, 2008 |
| Furey-Gupta | `furey(Q, a, A)` | Furey & Gupta, 2001 |
| WHAT | `what(Q, BFImax, a)` | Lim et al., 2005 (alias for Eckhardt) |

**gamma=1 family** (signal processing based):

| Method | Function | Reference |
|---|---|---|
| Lyne-Hollick | `lh(Q)` / `lh_multi(Q, num_pass=3)` | Lyne & Hollick, 1979; Nathan & McMahon, 1990 |
| Chapman | `chapman(Q, a)` | Chapman, 1991 |
| Willems | `willems(Q, a, w)` | Willems, 2009 |

**Variable gamma** (hybrid):

| Method | Function | Reference |
|---|---|---|
| IHACRES | `ihacres(Q, a, C, alpha_s)` | Jakeman & Hornberger, 1993 |

### Graphical / Recession-Based Methods

| Method | Function | Reference |
|---|---|---|
| UKIH | `ukih(Q, b_LH)` | UKIH, 1980 |
| Local minimum | `local(Q, b_LH, area)` | Sloto & Crouse, 1996 |
| Fixed interval | `fixed(Q, area)` | Sloto & Crouse, 1996 |
| Sliding interval | `slide(Q, area)` | Sloto & Crouse, 1996 |
| PART | `part(Q, area)` | Rutledge, 1998 |

### Tracer-Based Methods

| Method | Function | Reference |
|---|---|---|
| Conductivity Mass Balance | `cmb(Q, SC)` | Stewart et al., 2007 |

### Recession Analysis

| Function | Description |
|---|---|
| `bflow(Q)` | BFlow 3-pass filter + recession analysis (Arnold & Allen, 1999) |
| `bn77(Q, ...)` | Brutsaert-Nieber drought flow identification (Cheng et al., 2016) |

## Parameter Estimation

```python
import pybaseflow

data = pybaseflow.load_sample_data()
Q = data['Q']

# Estimate recession coefficient from the hydrograph
strict = pybaseflow.strict_baseflow(Q)
a = pybaseflow.recession_coefficient(Q, strict)

# Use the estimated recession coefficient with any filter
b = pybaseflow.eckhardt(Q, a, BFImax=0.8)
b = pybaseflow.chapman_maxwell(Q, a)
b = pybaseflow.willems(Q, a, w=0.5)
```

## CMB Calibration Bridge

Use specific conductance data to calibrate Eckhardt's BFImax:

```python
from pybaseflow.tracer import calibrate_eckhardt_from_cmb

cal = calibrate_eckhardt_from_cmb(Q, SC)
b = pybaseflow.eckhardt(Q, cal['a'], cal['BFImax'])
```

## BFlow / SWAT Integration

```python
result = pybaseflow.bflow(Q)
print(f"ALPHA_BF = {result['alpha_factor']:.4f}")  # for SWAT calibration
print(f"BFI = {result['BFI']:.3f}")
print(f"Baseflow days = {result['baseflow_days']:.1f}")
```

## USGS Data Retrieval

```python
from pybaseflow.io import fetch_usgs

# Fetch daily discharge
data = fetch_usgs('01013500', '2015-01-01', '2020-12-31')
Q = data['values']

# Fetch specific conductance (for CMB)
sc_data = fetch_usgs('01013500', '2015-01-01', '2020-12-31', parameter='sc')
```

## License

MIT
