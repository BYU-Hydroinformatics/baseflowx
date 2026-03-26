import numpy as np
from numba import njit


def clean_streamflow(series):
    """
    Cleans a streamflow time series by removing invalid values and keeping only years with at least 120 data points.

    Args:
        series (pandas.Series): The streamflow time series to be cleaned.

    Returns:
        tuple: A tuple containing the cleaned streamflow values and the corresponding dates.
    """
    date, Q = series.index, series.values.astype(float)
    has_value = np.isfinite(Q)
    date, Q = date[has_value], np.abs(Q[has_value])
    year_unique, counts = np.unique(date.year, return_counts=True)
    keep = np.isin(date.year, year_unique[counts >= 120])
    return Q[keep], date[keep]


def moving_average(x, w):
    res = np.convolve(x, np.ones(w)) / w
    return res[w - 1:-w + 1]


@njit
def multi_arange(starts, stops):
    """
    Generates a 1D numpy array containing all integers between the given start and stop values for each element in the input arrays.

    Args:
        starts (numpy.ndarray): A 1D numpy array of start values.
        stops (numpy.ndarray): A 1D numpy array of stop values, where each stop value corresponds to the start value at the same index.

    Returns:
        numpy.ndarray: A 1D numpy array containing all integers between the given start and stop values for each element in the input arrays.
    """
    pos = 0
    cnt = np.sum(stops - starts, dtype=np.int64)
    res = np.zeros((cnt,), dtype=np.int64)
    for i in range(starts.size):
        num = stops[i] - starts[i]
        res[pos:pos + num] = np.arange(starts[i], stops[i])
        pos += num
    return res


@njit
def backward(Q, b_LH, a):
    """
    Calculates the baseflow time series `b` from the discharge time series `Q` and the baseflow time series `b_LH` using a backward recursive approach.

    The function iterates through the discharge time series in reverse order, calculating the baseflow at each time step based on the baseflow at the next time step and the recession coefficient `a`. If the calculated baseflow exceeds the discharge at the current time step, the baseflow is set to the discharge.

    Args:
        Q (numpy.ndarray): The discharge time series.
        b_LH (numpy.ndarray): The baseflow time series.
        a (float): The recession coefficient.

    Returns:
        numpy.ndarray: The baseflow time series.
    """
    b = np.zeros(Q.shape[0])
    b[-1] = b_LH[-1]
    for i in range(Q.shape[0] - 1, 0, -1):
        b[i - 1] = b[i] / a
        if b[i] == 0:
            b[i - 1] = Q[i - 1]
        if b[i - 1] > Q[i - 1]:
            b[i - 1] = Q[i - 1]
    return b
