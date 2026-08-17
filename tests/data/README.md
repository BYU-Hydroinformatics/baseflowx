# Validation gage data

Daily-mean discharge (NWIS parameter 00060), 2015-01-01 to 2020-12-31,
for five USGS gages used by `tests/test_reference_implementations.py`
to validate baseflowx against the Xie et al. (2020) `baseflow` package
across multiple hydroclimatic settings.

Gages are drawn from the Tier A conductivity-mass-balance reference pool
screened for the manuscript revision (`review1/analysis/cmb_candidate_gages_screened.csv`
in the paper repo), chosen to span distinct ecoregions:

| site | name | ecoregion | drainage area (km²) |
|---|---|---|---|
| 10343500 | Sagehen Ck nr Truckee, CA | Western Mountains (snowmelt) | 27.6 |
| 14138900 | N.F. Bull Run R nr Multnomah Falls, OR | Pacific maritime | 21.7 |
| 04027000 | Bad R nr Odanah, WI | Great Lakes shield | 1619.8 |
| 02315550 | Suwannee R at Suwannee Springs, FL | SE coastal plain / karst | 6622.2 |
| 02204130 | Honey Ck nr Conyers, GA | SE Plains | 67.0 |

Fetched via `baseflowx.io.fetch_usgs`; all five records are complete
(0% missing) over the period. Regenerate with:

```bash
python -c "
from baseflowx.io import fetch_usgs
import csv
for site in ['10343500', '14138900', '04027000', '02315550', '02204130']:
    d = fetch_usgs(site, '2015-01-01', '2020-12-31', 'discharge')
    with open(f'tests/data/{site}.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['date', 'discharge_cfs'])
        for dt, v in zip(d['dates'], d['values']):
            w.writerow([dt.isoformat(), v])
"
```

The bundled Fish River sample (`baseflowx/data/fish_river.csv`, USGS
01013500, 2019-2020) is used as the sixth validation gage and needs no
separate fetch.

**Note:** the original semi-arid-plains representative, 07301500 (N.F.
Red River nr Carter, OK), was dropped after fetching -- 2015-2020 daily
discharge there includes days of zero flow (stream goes intermittent).
`recession_coefficient()`, `strict_baseflow()`, and `PART` all divide by
`Q` or take `log(Q)` and are not defined at zero flow, so this gage
cannot be used for validation as-is. This is a genuine package
limitation (baseflowx assumes perennial flow) worth noting in the
manuscript's Limitations section, separate from the validation-suite
question of which reference gage to use here. 02204130 (Honey Ck nr
Conyers, GA) was substituted as the SE Plains representative.

## Literature-benchmark gages

Used by `tests/test_literature_benchmarks.py` to check baseflowx against
values published in the peer-reviewed/USGS literature (not against another
codebase -- see that file's docstring).

| site | name | source | period fetched |
|---|---|---|---|
| 01555500 | East Mahantango Ck nr Dalmatia, PA | Risser et al. (2005), USGS SIR 2005-5038 | 1968-01-01 to 2001-12-31 (gap-free, matches the paper's period exactly) |
| 09333500 | Dirty Devil River, UT | Swanson et al. (2020), *Hydrogeol. J.* 29 | 2001-01-01 to 2019-12-31 (NWIS record starts 2001-05-02; 2 small 2-day gaps in 2014) |
| 09337500 | Escalante River, UT | Swanson et al. (2020) | 2001-01-01 to 2019-12-31, gap-free |
| 09404115 | Havasu Creek, AZ | Swanson et al. (2020) | 2001-01-01 to 2019-12-31, with a 2009-03-05 to 2011-09-30 gap (NWIS's current gap is larger than the paper's stated 2010-only gap -- likely a record revision since 2020) |
| 09403850 | Kanab Creek, AZ | Swanson et al. (2020) | 2018-04-17 to 2019-12-31 only -- NWIS does not currently hold the paper's stated 2016-2019 range for this gage; comparison is computed but not asserted (see test) |
| 09402300 | Little Colorado River, AZ | Swanson et al. (2020) | 2001-01-01 to 2019-12-31, gap-free (record starts 2003-05-25) |
| 09382000 | Paria River, UT/AZ | Swanson et al. (2020) | 2001-01-01 to 2019-12-31, gap-free |

Regenerate with:

```bash
python -c "
from baseflowx.io import fetch_usgs
import csv
sites = ['01555500', '09333500', '09337500', '09404115', '09403850', '09402300', '09382000']
ranges = {'01555500': ('1968-01-01', '2001-12-31'),
          '09403850': ('2016-01-01', '2019-12-31')}
for site in sites:
    start, end = ranges.get(site, ('2001-01-01', '2019-12-31'))
    d = fetch_usgs(site, start, end, 'discharge')
    with open(f'tests/data/{site}.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['date', 'discharge_cfs'])
        for dt, v in zip(d['dates'], d['values']):
            w.writerow([dt.isoformat(), v])
"
```
