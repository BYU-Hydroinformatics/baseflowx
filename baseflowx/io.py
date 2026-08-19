"""I/O utilities for fetching streamflow data from USGS NWIS."""

import json
from datetime import datetime
from urllib.request import urlopen
from urllib.parse import urlencode

import numpy as np

# NWIS parameter code mapping
_PARAM_CODES = {
    'discharge': '00060',
    'streamflow': '00060',
    'sc': '00095',
    'specific_conductance': '00095',
}

_NWIS_DV_URL = 'https://waterservices.usgs.gov/nwis/dv/'


def fetch_usgs(site_id, start_date, end_date, parameter='discharge',
               statistic='00003'):
    """Fetch daily values from the USGS NWIS Water Services API.

    Args:
        site_id (str): USGS site number (e.g., '01013500').
        start_date (str): Start date as 'YYYY-MM-DD'.
        end_date (str): End date as 'YYYY-MM-DD'.
        parameter (str): Parameter to fetch. Accepts convenience names
            ('discharge', 'streamflow', 'sc', 'specific_conductance') or
            a raw NWIS parameter code (e.g., '00060'). Default 'discharge'.
        statistic (str): NWIS statistic code selecting which daily series to
            return. Default '00003' (daily mean), which is almost always what
            is wanted. Other codes include '00001' (daily maximum) and '00002'
            (daily minimum). NWIS publishes several statistic series for the
            same parameter at many sites, so this must be specified explicitly
            to obtain a well-defined series.

    Returns:
        dict: Dictionary with keys:
            - 'dates': numpy array of datetime.date objects
            - 'values': numpy array of float64 values (NaN for missing)
            - 'units': str describing the measurement units
            - 'site_id': str echoing the requested site
            - 'parameter': str echoing the parameter code used
            - 'statistic': str echoing the NWIS statistic code used
            - 'qualifiers': list of qualifier strings per timestep

    Raises:
        ValueError: If the parameter name is unrecognized or no data is found.
        ConnectionError: If the NWIS API request fails.

    Example:
        >>> from baseflowx.io import fetch_usgs
        >>> data = fetch_usgs('01013500', '2015-01-01', '2015-12-31')
        >>> data['values'][:5]
        array([...])
    """
    # Resolve parameter code
    param_code = _PARAM_CODES.get(parameter, parameter)
    if not param_code.isdigit():
        raise ValueError(
            f"Unrecognized parameter: {parameter!r}. Use 'discharge', 'sc', "
            f"or a numeric NWIS code like '00060'.")

    # Build request URL
    params = {
        'format': 'json',
        'sites': site_id,
        'startDT': start_date,
        'endDT': end_date,
        'parameterCd': param_code,
        'statCd': statistic,
        'siteStatus': 'all',
    }
    url = _NWIS_DV_URL + '?' + urlencode(params)

    # Fetch data
    try:
        with urlopen(url, timeout=60) as response:
            raw = json.loads(response.read().decode('utf-8'))
    except Exception as exc:
        raise ConnectionError(f"NWIS request failed: {exc}") from exc

    # Parse JSON response
    try:
        ts = raw['value']['timeSeries']
        if not ts:
            raise ValueError(f"No data returned for site {site_id}, "
                             f"parameter {param_code}, "
                             f"{start_date} to {end_date}.")
        if len(ts) != 1:
            # Should not happen once statCd is pinned, but never guess: an
            # ambiguous response previously caused the first series (daily
            # maximum) to be returned in place of the requested statistic.
            found = []
            for s in ts:
                opts = s['variable'].get('options', {}).get('option', [])
                found.append(next((o.get('optionCode') for o in opts
                                   if o.get('name') == 'Statistic'), '?'))
            raise ValueError(
                f"NWIS returned {len(ts)} series for site {site_id}, "
                f"parameter {param_code}, statistic {statistic} "
                f"(statistic codes: {found}). Expected exactly one; "
                f"specify `statistic` explicitly.")
        series = ts[0]
        values_list = series['values'][0]['value']
        unit = series['variable']['unit']['unitCode']
    except (KeyError, IndexError) as exc:
        raise ValueError(
            f"Unexpected NWIS response structure: {exc}") from exc

    # Extract dates, values, and qualifiers
    dates = []
    values = []
    qualifiers = []
    for entry in values_list:
        dates.append(datetime.strptime(
            entry['dateTime'][:10], '%Y-%m-%d').date())
        val = entry['value']
        # NWIS uses -999999 or empty string for missing
        if val == '' or val == '-999999' or val is None:
            values.append(np.nan)
        else:
            values.append(float(val))
        qualifiers.append(entry.get('qualifiers', ['']))

    return {
        'dates': np.array(dates),
        'values': np.array(values, dtype=np.float64),
        'units': unit,
        'site_id': site_id,
        'parameter': param_code,
        'statistic': statistic,
        'qualifiers': qualifiers,
    }
