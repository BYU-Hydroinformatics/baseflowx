"""Validate baseflowx against results published in the peer-reviewed and
USGS literature, on named USGS gages, using software with no code lineage
to baseflowx or to the Xie et al. (2020) package used in
test_reference_implementations.py.

Two sources:

1. Risser, Gburek & Folmar (2005), USGS Scientific Investigations Report
   2005-5038. Runs the *original* USGS PART (Rutledge, 1993/1998) and HYSEP
   (Sloto & Crouse, 1996) programs -- not a Python reimplementation of any
   kind -- on East Mahantango Creek near Dalmatia, PA (USGS 01555500,
   162 mi^2), 1968-2001. Table 5 reports mean-annual base flow in
   inches/year for PART and all three HYSEP variants.

2. Swanson et al. (2020), Hydrogeology Journal 29, 723-736. Runs the R
   `EcoHydRology` package's three-pass Lyne-Hollick filter (Arnold & Allen,
   1999 "BFLOW" method; filter parameter 0.9, 3 passes) on six USGS-gaged
   Colorado Plateau tributaries. Table 4 reports base flow as a percentage
   of total tributary discharge over each gage's period of record.

Both were fetched fresh from NWIS for this test (tests/data/), which is
~5 years newer than the source papers -- gage records can be revised or
extended, so exact reproduction of the papers' stated periods of record is
not always possible. Where it isn't, that is noted explicitly rather than
silently worked around; see the Kanab Creek and Dirty Devil River cases
below.
"""

from datetime import date
from pathlib import Path

import numpy as np
import pytest

from baseflowx import separation

_DATA = Path(__file__).parent / "data"


def _load(site):
    import csv
    rows = list(csv.DictReader(open(_DATA / f"{site}.csv")))
    dates = [date.fromisoformat(r["date"]) for r in rows]
    Q = np.array([float(r["discharge_cfs"]) for r in rows])
    return dates, Q


def _contiguous_segments(dates, values):
    """Split a record at any date gap so the filter is never run across
    missing days as if they were adjacent."""
    segs, start = [], 0
    for i in range(1, len(dates)):
        if (dates[i] - dates[i - 1]).days != 1:
            segs.append(values[start:i])
            start = i
    segs.append(values[start:])
    return [s for s in segs if len(s) >= 10]


# ---------------------------------------------------------------------------
# Risser, Gburek & Folmar (2005) -- PART and HYSEP vs. the original USGS
# programs, East Mahantango Creek near Dalmatia PA (01555500), 1968-2001.
# ---------------------------------------------------------------------------

# Table 5, mean-annual, 1968-2001, inches/year.
_RISSER_2005 = {
    "part": 12.9,
    "local": 10.8,
    "slide": 12.3,
    "fixed": 12.2,
}
_DALMATIA_AREA_MI2 = 162.0


def _cfs_to_inches_per_year(b, n_years, area_mi2):
    area_ft2 = area_mi2 * 5280**2
    total_ft3 = np.sum(b) * 86400.0
    return (total_ft3 / area_ft2) * 12.0 / n_years


@pytest.mark.parametrize("method,tolerance_pct", [
    ("part", 2), ("fixed", 2), ("slide", 2), ("local", 5),
])
def test_matches_risser_2005_dalmatia(method, tolerance_pct):
    dates, Q = _load("01555500")
    assert dates[0] == date(1968, 1, 1) and dates[-1] == date(2001, 12, 31)
    assert len(Q) == 12419, "expected a gap-free 1968-2001 daily record"

    area_km2 = _DALMATIA_AREA_MI2 * 2.58999
    b_lh = separation.lh(Q)
    b = {
        "part": separation.part(Q, area_km2),
        "fixed": separation.fixed(Q, area_km2),
        "slide": separation.slide(Q, area_km2),
        "local": separation.local(Q, b_lh, area_km2),
    }[method]

    got = _cfs_to_inches_per_year(b, n_years=34, area_mi2=_DALMATIA_AREA_MI2)
    published = _RISSER_2005[method]
    pct_diff = 100 * abs(got - published) / published
    assert pct_diff < tolerance_pct, (
        f"{method}: baseflowx={got:.2f} in/yr vs. Risser et al. (2005) "
        f"Table 5={published:.2f} in/yr, {pct_diff:.1f}% apart "
        f"(tolerance {tolerance_pct}%)"
    )


# ---------------------------------------------------------------------------
# Swanson et al. (2020) -- three-pass Lyne-Hollick (beta=0.9) vs. R's
# EcoHydRology package, six Colorado Plateau tributaries.
# ---------------------------------------------------------------------------

# site: (published % of tributary discharge from Table 4, tolerance in
# percentage points, note)
_SWANSON_2020 = {
    "09337500": (43.0, 5, "Escalante River"),
    "09404115": (93.0, 5, "Havasu Creek"),
    "09402300": (69.0, 5, "Little Colorado River"),
    "09382000": (41.0, 5, "Paria River"),
    "09333500": (56.0, 10, "Dirty Devil River -- NWIS record here starts "
                 "2001-05-02, not Jan 2001 as the paper's period implies; "
                 "wider tolerance reflects the period mismatch, not a "
                 "code disagreement"),
    "09403850": (38.0, None, "Kanab Creek -- NWIS currently holds only "
                 "2018-04-17 onward (1.7 yr) vs. the paper's stated "
                 "2016-2019 (4 yr); not a meaningful comparison, computed "
                 "for the record but not asserted"),
}


@pytest.mark.parametrize("site", list(_SWANSON_2020))
def test_matches_swanson_2020_colorado_plateau(site):
    published_pct, tolerance_pts, note = _SWANSON_2020[site]
    dates, Q = _load(site)
    segments = _contiguous_segments(dates, Q)

    total_b, total_q = 0.0, 0.0
    for seg in segments:
        b = separation.lh_multi(seg, beta=0.9, num_pass=3)
        total_b += np.sum(b)
        total_q += np.sum(seg)
    got_pct = 100 * total_b / total_q

    if tolerance_pts is None:
        pytest.skip(f"{note} (baseflowx computed {got_pct:.1f}%, "
                     f"published {published_pct:.1f}%, not compared)")

    diff = abs(got_pct - published_pct)
    assert diff < tolerance_pts, (
        f"{site}: baseflowx={got_pct:.1f}% vs. Swanson et al. (2020) "
        f"Table 4={published_pct:.1f}%, {diff:.1f} points apart "
        f"(tolerance {tolerance_pts} points). {note}"
    )
