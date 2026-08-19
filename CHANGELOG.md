# Changelog

## 0.2.2 — 2026-08-19

### Fixed — data correctness (affects any prior use of `fetch_usgs`)

- **`fetch_usgs()` returned daily MAXIMUM values instead of daily mean.** The
  NWIS request did not pin a statistic code, and the parser took the first
  series returned. NWIS publishes several statistic series per parameter and
  orders them maximum, minimum, mean, median, so any site with a maximum series
  received daily maxima.

  Specific conductance is affected far more often than discharge, because
  maximum and minimum series are routinely published for conductance while most
  gages publish only daily mean discharge. In a 31-catchment tracer study, 26
  gages received maximum conductance and 1 received maximum discharge.

  `fetch_usgs()` now takes a `statistic` argument defaulting to `'00003'`
  (daily mean), pins it in the request, echoes it in the returned dict, and
  raises rather than guessing if more than one series is ever returned.

  **Anyone who has fetched data with an earlier release should re-fetch.**

- **`local()` (HYSEP local minimum) interpolated linearly between turning
  points; it now interpolates in log space,** as the original USGS Fortran does.
  Against the published benchmark of Risser, Gburek & Folmar (2005), linear
  interpolation gave 11.12 in/yr against their 10.8 (3.0% error); log-space
  gives 10.77 (0.3%), in line with the other HYSEP and PART comparisons.

## 0.2.1

Initial release described in the accompanying manuscript.
