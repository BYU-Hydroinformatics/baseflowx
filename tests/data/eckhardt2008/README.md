# Eckhardt (2008) 65-gage benchmark

Used by `tests/test_eckhardt2008_benchmark.py`. See that file's docstring
for the full rationale (why correlation, not exact match, is the right
comparison metric here).

- `published_bfi.csv` — Table 1 of Eckhardt, K. (2008), "A comparison of
  baseflow indices, which were calculated with seven different baseflow
  separation methods," *J. Hydrol.* 352, 168–173. Transcribed by hand;
  65 named USGS gages x 7 methods (HYSEP1=Fixed, HYSEP2=Sliding,
  HYSEP3=Local Minimum, PART, BFLOW=single-pass Lyne-Hollick, UKIH,
  Eckhardt).
- `{site}.csv` — daily-mean discharge (NWIS 00060/00003) for each of the
  65 gages, fetched via `baseflowx.io.fetch_usgs` over each gage's own
  actual period of record (`_best_windows.json`), not necessarily the
  period Eckhardt (2008) used (not stated in the paper). Many of these
  are historic Great Lakes-basin gages discontinued decades ago —
  periods of record range from 1942–1964 (04017000) to as short as
  1952–1957 (04032500, ~5.4 years). 8 of the 65 records have an internal
  gap; the test suite splits at gaps rather than treating non-adjacent
  days as consecutive.
- `_best_windows.json` — per-site {begin_date, end_date, count} for the
  parameter-00060/stat-00003 series with the most observations, from
  NWIS's `seriesCatalogOutput`. Regenerate via the site service:
  `https://waterservices.usgs.gov/nwis/site/?format=rdb&sites=<comma-list>&seriesCatalogOutput=true&parameterCd=00060`
- `_areas.json` — drainage area (mi²) per gage from NWIS site metadata
  (`drain_area_va`), needed for the HYSEP/PART interval calculation.

51–63 of 65 gages produce a result per method (PART and Eckhardt need
`Q > 0` throughout and enough `strict_baseflow()` points to fit a
recession constant; a handful of gages don't qualify).
