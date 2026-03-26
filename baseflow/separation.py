import numpy as np
from numba import njit, prange


def boughton(Q, a, C, initial_method='Q0', return_exceed=False):
    """
    Boughton doulbe-parameter filter (Boughton, 2004)
    Boughton W.C. (1993) - A hydrograph-based model for estimating water yield of ungauged catchments. Institute of Engineers Australia National Conference. Publ. 93/14, pp. 317-324.
    
    Args:
        Q (np.array): streamflow
        a (float): recession coefficient
        C (float): calibrated in baseflow.param_estimate
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be provided to directly set the initial baseflow value.
            Default is 'Q0'.
        return_exceed (bool, optional): if True, returns the number of times the
            baseflow exceeds the streamflow.

    Returns:
        b (np.array): baseflow
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    for i in range(Q.shape[0] - 1):
        b[i + 1] = a / (1 + C) * b[i] + C / (1 + C) * Q[i + 1]
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1

    return b

def chapman_maxwell(Q, a, initial_method='Q0' , return_exceed=False):
    """
    CM filter (Chapman & Maxwell, 1996)
    Chapman, T. G., Maxwell, A. I. (1996) - Baseflow separation - comparison of numerical methods with tracer experiments, in Hydrol. and Water Resour. Symp., Institution of Engineers Australia, Hobart. pp. 539-545.
    
    Args:
        Q (np.array): streamflow
        a (float): recession coefficient
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be oythoprovided to directly set the initial baseflow value.
            Default is 'Q0'.

    Returns:
        b (np.array): baseflow
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

   # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    for i in range(Q.shape[0] - 1):
        b[i + 1] = a / (2 - a) * b[i] + (1 - a) / (2 - a) * Q[i + 1]
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b

def chapman(Q, a = 0.925, initial_method='Q0', return_exceed=False):
    """Chapman filter (Chapman, 1991)
    Chapman, Tom G. "Comment on 'Evaluation of Automated Techniques for Base Flow and Recession Analyses' by R. J. Nathan and T. A. McMahon." Water Resources Research 27, no. 7 (1991): 1783–84. https://doi.org/10.1029/91WR01007.

    Args:
        Q (np.array): streamflow
        a (float): recession coefficient
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be provided to directly set the initial baseflow value.
            Default is 'Q0'.
        return_exceed (bool, optional): if True, returns the number of times the
            baseflow exceeds the streamflow.
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    for i in range(Q.shape[0] - 1):
        b[i + 1] = (3 * a - 1) / (3 - a) * b[i] + (1 - a) / (3 - a) * (Q[i + 1] + Q[i])
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b

def eckhardt(Q, a, BFImax, initial_method='Q0', return_exceed=False):
    """Eckhardt filter (Eckhardt, 2005)
    Eckhardt, K. “How to Construct Recursive Digital Filters for Baseflow Separation.” Hydrological Processes 19, no. 2 (2005): 507–15. https://doi.org/10.1002/hyp.5675.

    Args:
        Q (np.array): streamflow
        a (float): recession coefficient
        BFImax (float): maximum value of baseflow index (BFI)
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be provided to directly set the initial baseflow value.
            Default is 'Q0'.
        return_exceed (bool, optional): if True, returns the number of times the
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

   # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    for i in range(Q.shape[0] - 1):
        b[i + 1] = ((1 - BFImax) * a * b[i] + (1 - a) * BFImax * Q[i + 1]) / (1 - a * BFImax)
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b

def ewma(Q, e, initial_method='Q0', return_exceed=False):
    """exponential weighted moving average (EWMA) filter (Tularam & Ilahee, 2008)
    Tularam, Gurudeo Anand, and Mahbub Ilahee. “Exponential Smoothing Method of Base Flow Separation and Its Impact on Continuous Loss Estimates.” American Journal of Environmental Sciences 4, no. 2 (April 30, 2008): 136–44. https://doi.org/10.3844/ajessp.2008.136.144.
    
    Args:
        Q (np.array): streamflow
        e (float): smoothing parameter
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be provided to directly set the initial baseflow value.
            Default is 'Q0'.
        return_exceed (bool, optional): if True, returns the number of times the
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    for i in range(Q.shape[0] - 1):
        b[i + 1] = (1 - e) * b[i] + e * Q[i + 1]
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b

def fixed(Q, area=None):
    """Fixed interval graphical method from HYSEP program (Sloto & Crouse, 1996)
    Sloto, R. A., & Crouse, M. Y. (1996). HYSEP: A Computer Program for Streamflow Hydrograph Separation and Analysis (96-4040). Reston, VA: U.S. Geological Survey. https://doi.org/10.3133/wri964040.
    
    Args:
        Q (np.array): streamflow
        area (float): basin area in km^2
    """
    inN = hysep_interval(area)
    return fixed_interpolation(Q, inN)


@njit
def fixed_interpolation(Q, inN):
    b = np.zeros(Q.shape[0])
    n = Q.shape[0] // inN
    for i in prange(n):
        b[inN * i:inN * (i + 1)] = np.min(Q[inN * i:inN * (i + 1)])
    if n * inN != Q.shape[0]:
        b[n * inN:] = np.min(Q[n * inN:])
    return b

def furey(Q, a, A, initial_method='Q0', return_exceed=False):
    """Furey digital filter (Furey & Gupta, 2001, 2003)
    Furey, Peter R., and Vijay K. Gupta. “A Physically Based Filter for Separating Base Flow from Streamflow Time Series.” Water Resources Research 37, no. 11 (2001): 2709–22. https://doi.org/10.1029/2001WR000243.

    Args:
        Q (np.array): streamflow
        a (float): recession coefficient
        A (float): calibrated in baseflow.param_estimate
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be provided to directly set the initial baseflow value.
            Default is 'Q0'.
        return_exceed (bool, optional): if True, returns the number of times the
            baseflow exceeds the streamflow.
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    for i in range(Q.shape[0] - 1):
        b[i + 1] = (a - A * (1 - a)) * b[i] + A * (1 - a) * Q[i]
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b


def hyd_run(streamflow, k=0.9, passes=4):
    """
    Separates baseflow from a streamflow hydrograph using a digital filter method.

    Args:
        streamflow (numpy.ndarray): A numpy array of streamflow values in chronological order.
        k (float, optional): A filter coefficient between 0 and 1 (typically 0.9). Defaults to 0.9.
        passes (int, optional): Number of times the filter passes through the data (typically 4). Defaults to 4.

    Returns:
        numpy.ndarray: A numpy array of baseflow values.

    Example:
        >>> import numpy as np
        >>> streamflow = np.array([10, 15, 20, 18, 12])
        >>> baseflow = hyd_run(streamflow)
        >>> print(baseflow)
        [10.         10.90909091 13.27272727 15.18181818 14.36363636]
    """
    # Convert to numpy array and handle NaN values
    Q = np.array(streamflow)
    Q = Q[~np.isnan(Q)]

    # Initialize baseflow list
    baseflow = np.zeros_like(Q)
    baseflow[0] = Q[0]  # Set first baseflow value to first streamflow value

    for p in range(1, passes + 1):
        # Forward and backward pass
        if p % 2 == 1:
            start, end, step = 0, len(Q), 1
        else:
            start, end, step = len(Q) - 1, -1, -1

        for i in range(start + step, end, step):
            tmp = k * baseflow[i - step] + (1 - k) * (Q[i] + Q[i - step]) / 2
            baseflow[i] = min(tmp, Q[i])

    return baseflow

def lh_multi(Q, beta=0.925, num_pass=2, return_exceed=False):
    """
    Applies a low-pass filter to the input time series `Q` using the Lyne-Hollick (LH) recursive digital filter.
    
    The filter is applied in multiple passes, with the number of passes controlled by the `num_pass` parameter. The filter uses a smoothing parameter `beta` to control the degree of filtering.

    Spongberg, M. E. “Spectral Analysis of Base Flow Separation with Digital Filters.” Water Resources Research 36, no. 3 (2000): 745–52. https://doi.org/10.1029/1999WR900303.
   
    If `return_exceed` is True, the function will also return the number of times the filtered output `b` exceeds the original input `Q`.
    
    Lyne, V. and Hollick, M. (1979) Stochastic Time-Variable Rainfall-Runoff Modeling. Institute of Engineers Australia National Conference, 89-93. 
    
    Spongberg, M. E. “Spectral Analysis of Base Flow Separation with Digital Filters.” Water Resources Research 36, no. 3 (2000): 745–52. https://doi.org/10.1029/1999WR900303.

    Args:
        Q (numpy.ndarray): The input time series to be filtered.
        beta (float, optional): The smoothing parameter for the LH filter, between 0 and 1. Defaults to 0.925.
        num_pass (int, optional): The number of filter passes to apply. Defaults to 2.
        return_exceed (bool, optional): If True, the function will return the number of times the filtered output exceeds the original input. Defaults to False.
    
    Returns:
        numpy.ndarray: The filtered output time series.
        int (optional): The number of times the filtered output exceeds the original input, if `return_exceed` is True.
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    b[0] = Q[0]

    for n in range(num_pass):
        if n != 0:
            b = np.flip(b, axis=0)
            Q = b.copy()
        
        for i in range(Q.shape[0] - 1):
            b[i + 1] = beta * b[i] + (1 - beta) / 2 * (Q[i] + Q[i + 1])
            if b[i + 1] > Q[i + 1]:
                b[i + 1] = Q[i + 1]
                if return_exceed:
                    b[-1] += 1

    if num_pass % 2 == 0:
        b = np.flip(b, axis=0)

    return b

def lh(Q, beta=0.925, return_exceed=False):
    """LH digital filter (Lyne & Hollick, 1979)
    Lyne, V. and Hollick, M. (1979) Stochastic Time-Variable Rainfall-Runoff Modeling. Institute of Engineers Australia National Conference, 89-93.
    
    Args:
        Q (np.array): streamflow
        beta (float): filter parameter, 0.925 recommended by (Nathan & McMahon, 1990)
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # first pass
    b[0] = Q[0]
    for i in range(Q.shape[0] - 1):
        b[i + 1] = beta * b[i] + (1 - beta) / 2 * (Q[i] + Q[i + 1])
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1

    # second pass
    b1 = np.copy(b)
    for i in range(Q.shape[0] - 2, -1, -1):
        b[i] = beta * b[i + 1] + (1 - beta) / 2 * (b1[i + 1] + b1[i])
        if b[i] > b1[i]:
            b[i] = b1[i]
            if return_exceed:
                b[-1] += 1
    return b

def local(Q, b_LH, area=None, return_exceed=False):
    """Local minimum graphical method from HYSEP program (Sloto & Crouse, 1996)

    Args:
        Q (np.array): streamflow
        area (float): basin area in km^2
    """
    idx_turn = local_turn(Q, hysep_interval(area))
    if idx_turn.shape[0] < 3:
        raise IndexError('Less than 3 turning points found')
    b = linear_interpolation(Q, idx_turn, return_exceed=return_exceed)
    b[:idx_turn[0]] = b_LH[:idx_turn[0]]
    b[idx_turn[-1] + 1:] = b_LH[idx_turn[-1] + 1:]
    return b


def hysep_interval(area):
    # The duration of surface runoff is calculated from the empirical relation:
    # N=A^0.2, (1) where N is the number of days after which surface runoff ceases,
    # and A is the drainage area in square miles (Linsley and others, 1982, p. 210).
    # The interval 2N* used for hydrograph separations is the odd integer between
    # 3 and 11 nearest to 2N (Pettyjohn and Henning, 1979, p. 31).
    if area is None:
        N = 5
    else:
        N = np.power(0.3861022 * area, 0.2)
    inN = np.ceil(2 * N)
    if np.mod(inN, 2) == 0:
        inN = np.ceil(2 * N) - 1
    inN = np.int64(min(max(inN, 3), 11))
    return inN


@njit
def local_turn(Q, inN):
    idx_turn = np.zeros(Q.shape[0], dtype=np.int64)
    for i in prange(np.int64((inN - 1) / 2), np.int64(Q.shape[0] - (inN - 1) / 2)):
        if Q[i] == np.min(Q[np.int64(i - (inN - 1) / 2):np.int64(i + (inN + 1) / 2)]):
            idx_turn[i] = i
    return idx_turn[idx_turn != 0]

def slide(Q, area):
    """Slide interval graphical method from HYSEP program (Sloto & Crouse, 1996)
    Sloto, R. A., & Crouse, M. Y. (1996). HYSEP: A Computer Program for Streamflow Hydrograph Separation and Analysis (96-4040). Reston, VA: U.S. Geological Survey. https://doi.org/10.3133/wri964040.
    
    Args:
        Q (np.array): streamflow
        area (float): basin area in km^2
    """
    inN = hysep_interval(area)
    return slide_interpolation(Q, inN)

def slide_interpolation(Q, inN):
    b = np.zeros(Q.shape[0])
    for i in prange(np.int64((inN - 1) / 2), np.int64(Q.shape[0] - (inN - 1) / 2)):
        b[i] = np.min(Q[np.int64(i - (inN - 1) / 2):np.int64(i + (inN + 1) / 2)])
    b[:np.int64((inN - 1) / 2)] = np.min(Q[:np.int64((inN - 1) / 2)])
    b[np.int64(Q.shape[0] - (inN - 1) / 2):] = np.min(Q[np.int64(Q.shape[0] - (inN - 1) / 2):])
    return b

def ukih(Q, b_LH, return_exceed=False):
    """graphical method developed by UK Institute of Hydrology (UKIH, 1980)
    Aksoy, Hafzullah, Ilker Kurt, and Ebru Eris. “Filtered Smoothed Minima Baseflow Separation Method.” Journal of Hydrology 372, no. 1 (June 15, 2009): 94–101. https://doi.org/10.1016/j.jhydrol.2009.03.037.

    Args:
        Q (np.array): streamflow
        return_exceed (bool, optional): if True, returns the number of times the
            baseflow exceeds the streamflow.
    """
    N = 5
    block_end = Q.shape[0] // N * N
    idx_min = np.argmin(Q[:block_end].reshape(-1, N), axis=1)
    idx_min = idx_min + np.arange(0, block_end, N)
    idx_turn = ukih_turn(Q, idx_min)
    if idx_turn.shape[0] < 3:
        raise IndexError('Less than 3 turning points found')
    b = linear_interpolation(Q, idx_turn, return_exceed=return_exceed)
    b[:idx_turn[0]] = b_LH[:idx_turn[0]]
    b[idx_turn[-1] + 1:] = b_LH[idx_turn[-1] + 1:]
    return b


@njit
def ukih_turn(Q, idx_min):
    idx_turn = np.zeros(idx_min.shape[0], dtype=np.int64)
    for i in prange(idx_min.shape[0] - 2):
        if ((0.9 * Q[idx_min[i + 1]] < Q[idx_min[i]]) &
                (0.9 * Q[idx_min[i + 1]] < Q[idx_min[i + 2]])):
            idx_turn[i] = idx_min[i + 1]
    return idx_turn[idx_turn != 0]


@njit
def linear_interpolation(Q, idx_turn, return_exceed=False):
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    n = 0
    for i in range(idx_turn[0], idx_turn[-1] + 1):
        if i == idx_turn[n + 1]:
            n += 1
            b[i] = Q[i]
        else:
            b[i] = Q[idx_turn[n]] + (Q[idx_turn[n + 1]] - Q[idx_turn[n]]) / \
                (idx_turn[n + 1] - idx_turn[n]) * (i - idx_turn[n])
        if b[i] > Q[i]:
            b[i] = Q[i]
            if return_exceed:
                b[-1] += 1
    return b

def what(streamflow, BFImax, alpha):
    """
    Separates baseflow and quickflow from a streamflow time series using the WHAT method.

    Args:
        streamflow (numpy.ndarray): A numpy array of streamflow values.
        BFImax (float): The maximum baseflow index (BFI) value.
        alpha (float): A filter parameter.

    Returns:
        tuple: A tuple containing two numpy arrays: baseflow and quickflow.

    Example:
        >>> import numpy as np
        >>> streamflow = np.array([10, 15, 20, 18, 12])
        >>> baseflow, quickflow = what(streamflow, 0.8, 0.98)
        >>> print(baseflow)
        [ 9.8         12.74        15.68        16.544       15.3712    ]
        >>> print(quickflow)
        [ 0.2          2.26         4.32         1.456       -3.3712    ]
    """
    baseflow = np.zeros_like(streamflow)

    for t in range(1, len(streamflow)):
        baseflow[t] = ((1 - BFImax) * alpha * baseflow[t-1] + (1 - alpha) * BFImax * streamflow[t]) / (1 - alpha * BFImax)

    quickflow = streamflow - baseflow

    return baseflow

def willems(Q, a, w, initial_method='Q0', return_exceed=False):
    """digital filter (Willems, 2009)

    Args:
        Q (np.array): streamflow
        a (float): recession coefficient
        w (float): case-speciﬁc average proportion of the quick ﬂow
                   in the streamflow, calibrated in baseflow.param_estimate
        initial_method (str or float, optional): method to calculate the initial baseflow value.
            Accepted string values are:
            - 'Q0': Use Q[0] as the initial baseflow value.
            - 'min': Use np.min(Q) as the initial baseflow value.
            - 'LH': Calculate the initial baseflow value using the LH method.
            Alternatively, a float value can be provided to directly set the initial baseflow value.
            Default is 'Q0'.
        return_exceed (bool, optional): if True, returns the number of times the
            baseflow exceeds the streamflow.
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # Set initial value for b based on the specified method
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            b[0] = Q[0]
        elif initial_method == 'min':
            b[0] = np.min(Q)
        elif initial_method == 'LH':
            b[0] = lh(Q)[0]  # Calculate the initial value using the LH method
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    else:
        b[0] = initial_method

    v = (1 - w) * (1 - a) / (2 * w)
    for i in range(Q.shape[0] - 1):
        b[i + 1] = (a - v) / (1 + v) * b[i] + v / (1 + v) * (Q[i] + Q[i + 1])
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b

def strict_baseflow(Q, ice=None, quantile=0.9):
    """
    Identify the strict baseflow component of a flow time series.
    
    This function applies a series of heuristic rules to identify the strict baseflow
    component of a flow time series. The rules are based on the behavior of the
    derivative of the flow time series, as well as the magnitude of the flow values.
    
    The function returns a boolean mask indicating the time steps that correspond to
    the strict baseflow component.
    
    Parameters:
        Q (numpy.ndarray): The flow time series.
        ice (numpy.ndarray, optional): A boolean mask indicating time steps with ice
            conditions, which can invalidate the groundwater-baseflow relationship.
        quantile (float, optional): The quantile value used to identify major events.
            Default is 0.9 (90th percentile).
    
    Returns:
        numpy.ndarray: A boolean mask indicating the time steps that correspond to
            the strict baseflow component.
    """
    dQ = (Q[2:] - Q[:-2]) / 2

    # 1. flow data associated with positive and zero values of dy / dt
    wet1 = np.concatenate([[True], dQ >= 0, [True]])

    # 2. previous 2 points before points with dy/dt≥0, as well as the next 3 points
    idx_first = np.where(wet1[1:].astype(int) - wet1[:-1].astype(int) == 1)[0] + 1
    idx_last = np.where(wet1[1:].astype(int) - wet1[:-1].astype(int) == -1)[0]
    idx_before = np.repeat([idx_first], 2) - np.tile(range(1, 3), idx_first.shape)
    idx_next = np.repeat([idx_last], 3) + np.tile(range(1, 4), idx_last.shape)
    idx_remove = np.concatenate([idx_before, idx_next])
    wet2 = np.full(Q.shape, False)
    wet2[idx_remove.clip(min=0, max=Q.shape[0] - 1)] = True

    # 3. five data points after major events (quantile)
    growing = np.concatenate([[True], (Q[1:] - Q[:-1]) >= 0, [True]])
    idx_major = np.where((Q >= np.quantile(Q, quantile)) & growing[:-1] & ~growing[1:])[0]
    idx_after = np.repeat([idx_major], 5) + np.tile(range(1, 6), idx_major.shape)
    wet3 = np.full(Q.shape, False)
    wet3[idx_after.clip(min=0, max=Q.shape[0] - 1)] = True

    # 4. flow data followed by a data point with a larger value of -dy / dt
    wet4 = np.concatenate([[True], dQ[1:] - dQ[:-1] < 0, [True, True]])

    # dry points, namely strict baseflow
    dry = ~(wet1 + wet2 + wet3 + wet4)

    # avoid ice conditions which invalidate the groundwater-baseflow relationship
    if ice is not None:
        dry[ice] = False

    return dry




def bn77(Q, L_min, snow_freeze_period, observational_precision, quantile=0.9):
    """
    Identifies the drought flow points in the discharge time series.
    Cheng, Lei, Lu Zhang, and Wilfried Brutsaert. “Automated Selection of Pure Base Flows from Regular Daily Streamflow Data: Objective Algorithm.” Journal of Hydrologic Engineering 21, no. 11 (November 1, 2016): 06016008. https://doi.org/10.1061/(ASCE)HE.1943-5584.0001427.
    
    Args:
        Q (numpy.ndarray): The discharge time series.
        L_min (int): Minimum number of points to be eliminated at the beginning and end of recession episode.
        snow_freeze_period (tuple): Start and end indices of the snow and/or freeze period.
        observational_precision (float): Observational precision threshold.
        quantile (float): Quanti le for identifying major events, default is 0.9.

    Returns:
        numpy.ndarray: The indices of the drought flow points.
    """
    # Step 1: Time series
    S = _estimate_recession_slope(Q)
    
    # Step 2: Recession episodes
    recession_episodes = _identify_recession_episodes(S, L_min)
    
    # Step 3: Drought flow points
    drought_flow_points = _eliminate_points(recession_episodes, L_min, snow_freeze_period, observational_precision, Q, quantile)
    
    return drought_flow_points

@njit
# step 1
def _estimate_recession_slope(Q): # FIXME: fix the slopes of the first and the last timestamp
    """
    Estimates the recession slope S(t)
    """
    N = len(Q)
    S = np.zeros(N)
    for i in range(N-1):
        if i == 0:
            pass
        else:
            S[i] = (Q[i-1] - Q[i+1]) / 2
    return S


@njit
# step 2
def _identify_recession_episodes(S, L_min):
    """
    Identifies the preliminary recession episodes.

    Args:
        S (numpy.ndarray): The recession slope time series.
        L_min (int): Minimum length of a recession episode.

    Returns:
        list: A list of arrays, each containing the indices of a preliminary recession episode.
    """
    i = 0
    N = len(S)
    recession_episodes = []

    while i < N - 1:
        if S[i] <= 0 and S[i+1] > 0:
            episode_start = i + 1 ##
            l = 0
            i=i+1 ##
            # how many points's S > 0
            while i < N - 1 and S[i] > 0:
                i += 1
                l += 1
            if l >= L_min:
                recession_episodes.append(np.arange(episode_start, i))
            else:
                i = i + l
        else:
            i += 1

    return recession_episodes


@njit
# step 3
def _eliminate_points(recession_episodes, L_min, snow_freeze_period, observational_precision, Q, S, quantile):
    """
    Eliminates the points at the beginning, end, and during the snow/freeze period, as well as anomalous and low-precision points.

    Args:
        recession_episodes (list): A list of arrays, each containing the indices of a preliminary recession episode.
        L_min (int): Minimum number of points to be eliminated at the beginning and end of recession episode.
        snow_freeze_period (tuple): Start and end indices of the snow and/or freeze period.
        observational_precision (float): Observational precision threshold.
        Q (numpy.ndarray): The discharge time series.
        S (numpy.ndarray): The recession slope time series.
        quantile (float): Quantile for identifying major events.

    Returns:
        numpy.ndarray: The indices of the drought flow points.
    """
    drought_flow_points = []
    major_event_threshold = np.quantile(Q, quantile)

    for episode in recession_episodes:
        # Check the value of the first point in the episode
        if Q[episode[0]] > major_event_threshold:
            # C4: Remove the first three points if the first point's value is greater than the 90% quantile
            if len(episode) > 3:
                episode = episode[3:]
            else:
                continue  # Skip episodes that are too short
        else:
            # C3: Remove the first two points otherwise
            if len(episode) > 2:
                episode = episode[2:]
            else:
                continue  # Skip episodes that are too short

        # C5: Remove the last point of each episode
        if len(episode) > 1:
            episode = episode[:-1]
        else:
            continue  # Skip episodes that are too short

        # C6: Remove points where Si/Si-1 >= 2
        if len(episode) > 1:
            episode = episode[1:][S[episode[1:]] / S[episode[:-1]] < 2]

        # C7: Remove points where Si < Si+1
        if len(episode) > 1:
            episode = episode[:-1][S[episode[:-1]] >= S[episode[1:]]]
        
        # C8: Eliminating the data points during the snow and/or freeze periods (i.e. C8)
        episode = episode[(episode < snow_freeze_period[0]) | (episode > snow_freeze_period[1])]
        
        # C9: Eliminating the data points of which Q(t) are smaller than observational precision (i.e. C9)
        episode = episode[Q[episode] >= observational_precision]
        
        drought_flow_points.append(episode)

    return np.concatenate(drought_flow_points) if drought_flow_points else np.array([])
