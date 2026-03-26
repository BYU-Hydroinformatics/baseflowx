"""pybaseflow — a comprehensive Python toolkit for baseflow separation."""

__version__ = "0.1.0"

import csv
from datetime import date
from pathlib import Path

import numpy as np

from pybaseflow.separation import *  # noqa: F401,F403
from pybaseflow.estimate import *  # noqa: F401,F403
from pybaseflow.utils import *  # noqa: F401,F403
from pybaseflow.tracer import cmb, estimate_endmembers, calibrate_eckhardt_from_cmb  # noqa: F401

_DATA_DIR = Path(__file__).parent / 'data'


def load_sample_data():
    """Load the bundled Fish River sample dataset.

    Returns daily discharge for USGS site 01013500 (Fish River near Fort Kent,
    Maine) for 2019-2020.

    Returns:
        dict: Dictionary with keys:
            - 'dates': numpy array of datetime.date objects
            - 'Q': numpy array of discharge in ft³/s
            - 'site_id': '01013500'
            - 'units': 'ft3/s'
    """
    path = _DATA_DIR / 'fish_river.csv'
    dates = []
    values = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            dates.append(date.fromisoformat(row['date']))
            val = row['discharge_cfs']
            values.append(float(val) if val else np.nan)
    return {
        'dates': np.array(dates),
        'Q': np.array(values, dtype=np.float64),
        'site_id': '01013500',
        'units': 'ft3/s',
    }
