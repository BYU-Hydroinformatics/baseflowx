"""Shared fixtures for baseflowx's test suite."""

import csv
from pathlib import Path

import numpy as np
import pytest

_HERE = Path(__file__).parent
_PACKAGE_DATA = _HERE.parent / "baseflowx" / "data"

# (site_id, csv path, drainage area km^2, ecoregion)
VALIDATION_GAGES = [
    ("01013500", _PACKAGE_DATA / "fish_river.csv", 2253.0, "Northeast (snowmelt headwater)"),
    ("10343500", _HERE / "data" / "10343500.csv", 27.6, "Western Mountains (snowmelt)"),
    ("14138900", _HERE / "data" / "14138900.csv", 21.7, "Pacific maritime"),
    ("04027000", _HERE / "data" / "04027000.csv", 1619.8, "Great Lakes shield"),
    ("02315550", _HERE / "data" / "02315550.csv", 6622.2, "SE coastal plain / karst"),
    ("02204130", _HERE / "data" / "02204130.csv", 67.0, "SE Plains"),
]


def _load_discharge(csv_path):
    values = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            values.append(float(row["discharge_cfs"]))
    return np.array(values, dtype=np.float64)


@pytest.fixture(scope="session")
def validation_gages():
    """Load all validation gages as a list of (site_id, Q, area_km2, ecoregion)."""
    out = []
    for site_id, path, area, ecoregion in VALIDATION_GAGES:
        Q = _load_discharge(path)
        assert np.all(Q > 0), f"{site_id}: non-positive discharge present"
        out.append((site_id, Q, area, ecoregion))
    return out


@pytest.fixture(scope="session")
def fish_river_Q():
    return _load_discharge(_PACKAGE_DATA / "fish_river.csv")
