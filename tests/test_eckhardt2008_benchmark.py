"""Validate seven baseflowx methods -- most importantly Eckhardt itself --
against Eckhardt (2008), J. Hydrol. 352, 168-173, "A comparison of baseflow
indices, which were calculated with seven different baseflow separation
methods."

This is the strongest external reference available to baseflowx for two
reasons. First, it is Eckhardt's *own* paper, run with his own custom
software ("a free computer program... provided by the author on request")
-- there is no code relationship whatsoever to baseflowx, Xie(2020),
EcoHydRology, or WHAT. Second, unlike the single- or six-gage comparisons
in test_literature_benchmarks.py, Table 1 reports per-method BFI for 65
individually-named USGS gages, which is enough gages to make a
*correlation* the right comparison metric rather than requiring exact
per-gage agreement.

Correlation is the right metric here for a specific reason: Eckhardt
derives the recession constant `a` via the Langbein (1938) correlation
method (fitting the upper envelope of a Q[k+1] vs Q[k] scatterplot during
recessions, with a 2% tolerance rule), which is a different technique from
baseflowx's own `recession_coefficient()` (5th-percentile of -dQ/Q at
`strict_baseflow()` points -- see docs/methods/recession-analysis.md,
"Percentile fit vs. log-log envelope fit"). A per-gage exact match is not
expected because the two `a` values will differ; what should hold is that
baseflowx's Eckhardt BFI, computed with baseflowx's *own* recession
constant, tracks the published Eckhardt BFI closely across many gages --
exactly the standard the paper itself uses (its Table 3 reports 0.85-0.98
correlation *among* the seven methods it compares, not exact agreement).

The other six methods in the paper (three HYSEP variants, PART, BFLOW --
recomputed by Eckhardt as a single-pass Lyne-Hollick filter, not the
three-pass version -- and UKIH) have no such `a`-estimation ambiguity, so
they serve as a sanity check on the correlation methodology itself: if
those six also show near-machine correlation, the Eckhardt result's
weaker (but still strong) correlation is credibly attributable to the
disclosed `a`-estimation difference, not to a bug in the Eckhardt
implementation.

Streamflow was fetched fresh from NWIS for each gage's own actual period
of record (tests/data/eckhardt2008/_best_windows.json), which generally
does not match whatever period the original 2005/2008 studies used (not
stated in either paper) -- many of these are historic Great Lakes gages
discontinued decades ago. BFI is a comparatively stable long-term
catchment characteristic, so this is expected to still correlate well; it
is one more reason exact agreement is not the bar being applied.
"""

import csv
import json
from datetime import date
from pathlib import Path

import numpy as np
import pytest

from baseflowx import estimate, separation

_DIR = Path(__file__).parent / "data" / "eckhardt2008"


def _load_published():
    return {r["site"]: r for r in csv.DictReader(open(_DIR / "published_bfi.csv"))}


def _load_areas():
    return json.load(open(_DIR / "_areas.json"))


def _load_gage(site):
    rows = list(csv.DictReader(open(_DIR / f"{site}.csv")))
    dates = [date.fromisoformat(r["date"]) for r in rows]
    Q = np.array([float(r["discharge_cfs"]) for r in rows])
    return dates, Q


def _segments(dates, Q, min_len=30):
    segs, start = [], 0
    for i in range(1, len(dates)):
        if (dates[i] - dates[i - 1]).days != 1:
            segs.append((dates[start:i], Q[start:i]))
            start = i
    segs.append((dates[start:], Q[start:]))
    return [(d, q) for d, q in segs if len(q) >= min_len]


def _bflow_1pass(Q, beta=0.925):
    """BFLOW as recomputed in Eckhardt (2008): single-pass Lyne-Hollick,
    not baseflowx's default 2-pass lh()."""
    n = len(Q)
    b = np.zeros(n)
    b[0] = Q[0]
    for i in range(n - 1):
        b[i + 1] = beta * b[i] + (1 - beta) / 2 * (Q[i] + Q[i + 1])
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
    return b


@pytest.fixture(scope="module")
def bfi_by_method():
    """Run all seven methods over all 65 gages; returns {method: {site: bfi}}."""
    published = _load_published()
    areas_mi2 = _load_areas()
    sites = list(published)

    out = {m: {} for m in
           ["HYSEP1", "HYSEP2", "HYSEP3", "PART", "BFLOW", "UKIH", "Eckhardt"]}

    for site in sites:
        dates, Q = _load_gage(site)
        area_km2 = areas_mi2[site] * 2.58999 if areas_mi2.get(site) else None
        segs = _segments(dates, Q)

        totals = {m: [0.0, 0.0] for m in out}
        for _, q in segs:
            b_lh = separation.lh(q)

            totals["BFLOW"][0] += np.sum(_bflow_1pass(q))
            totals["BFLOW"][1] += np.sum(q)
            totals["HYSEP1"][0] += np.sum(separation.fixed(q, area_km2))
            totals["HYSEP1"][1] += np.sum(q)
            totals["HYSEP2"][0] += np.sum(separation.slide(q, area_km2))
            totals["HYSEP2"][1] += np.sum(q)
            try:
                totals["HYSEP3"][0] += np.sum(separation.local(q, b_lh, area_km2))
                totals["HYSEP3"][1] += np.sum(q)
            except IndexError:
                pass
            try:
                totals["UKIH"][0] += np.sum(separation.ukih(q, b_lh))
                totals["UKIH"][1] += np.sum(q)
            except IndexError:
                pass
            if np.all(q > 0):
                totals["PART"][0] += np.sum(separation.part(q, area_km2))
                totals["PART"][1] += np.sum(q)
                strict = separation.strict_baseflow(q)
                if np.sum(strict) > 20:
                    a = estimate.recession_coefficient(q, strict)
                    if 0 < a < 1:
                        b_eck = separation.eckhardt(q, a, 0.8, initial_method="LH")
                        totals["Eckhardt"][0] += np.sum(b_eck)
                        totals["Eckhardt"][1] += np.sum(q)

        for m in out:
            sb, sq = totals[m]
            if sq > 0:
                out[m][site] = sb / sq

    return out


_PUBLISHED_COLUMN = {
    "HYSEP1": "HYSEP1", "HYSEP2": "HYSEP2", "HYSEP3": "HYSEP3",
    "PART": "PART", "BFLOW": "BFLOW", "UKIH": "UKIH", "Eckhardt": "Eckhardt",
}
# Minimum acceptable Pearson r and minimum gage count, per method.
_THRESHOLDS = {
    "HYSEP1": (0.97, 40), "HYSEP2": (0.97, 40), "HYSEP3": (0.95, 40),
    "PART": (0.97, 40), "BFLOW": (0.95, 40), "UKIH": (0.95, 40),
    # Looser bound for Eckhardt: baseflowx and Eckhardt (2008) estimate the
    # recession constant `a` by genuinely different methods (see module
    # docstring), so this is expected to correlate less tightly than the
    # six methods above, which have no such ambiguity.
    "Eckhardt": (0.80, 40),
}


@pytest.mark.parametrize("method", list(_THRESHOLDS))
def test_correlates_with_eckhardt_2008(bfi_by_method, method):
    published = _load_published()
    min_r, min_n = _THRESHOLDS[method]
    col = _PUBLISHED_COLUMN[method]

    xs, ys = [], []
    for site, bfx_bfi in bfi_by_method[method].items():
        xs.append(bfx_bfi)
        ys.append(float(published[site][col]))
    xs, ys = np.array(xs), np.array(ys)

    assert len(xs) >= min_n, (
        f"{method}: only {len(xs)} gages produced a result, expected >= {min_n}"
    )
    r = np.corrcoef(xs, ys)[0, 1]
    assert r >= min_r, (
        f"{method}: correlation with Eckhardt (2008) Table 1 is {r:.3f} "
        f"across {len(xs)} gages, below the {min_r} threshold"
    )
