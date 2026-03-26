import numpy as np
from numba import njit, prange
from pybaseflow.utils import moving_average, multi_arange, backward


def recession_coefficient(Q, strict):
    """
    Calculates the recession coefficient `K` from the given discharge `Q` and a boolean mask `strict` indicating which values to use.
    
    The recession coefficient `K` is calculated as follows:
    1. Extract the middle values of `Q` (`cQ`) and the centered finite difference of `Q` (`dQ`) using the `strict` mask.
    2. Sort `dQ / cQ` in descending order and take the value at the 5th percentile.
    3. Calculate `K` as the negative ratio of `cQ` to `dQ` at the selected index.
    4. Return the exponential of `-1 / K` as the final recession coefficient.
    
    Args:
        Q (numpy.ndarray): Array of discharge values.
        strict (numpy.ndarray): Boolean mask indicating which values of `Q` to use.
    
    Returns:
        float: The calculated recession coefficient.
    """
    cQ, dQ = Q[1:-1], (Q[2:] - Q[:-2]) / 2
    cQ, dQ = cQ[strict[1:-1]], dQ[strict[1:-1]]

    idx = np.argsort(-dQ / cQ)[np.floor(dQ.shape[0] * 0.05).astype(int)]
    K = - cQ[idx] / dQ[idx]
    return np.exp(-1 / K)


def param_calibrate(param_range, method, Q, b_LH, a):
    """
    Calibrates the parameters for a baseflow estimation method.
    
    Args:
        param_range (numpy.ndarray): The range of parameter values to test.
        method (callable): The baseflow estimation method to use.
        Q (numpy.ndarray): The discharge values.
        b_LH (numpy.ndarray): The low-flow baseflow values.
        a (float): The parameter for the baseflow estimation method.
    
    Returns:
        float: The optimal parameter value from the given range.
    """
    idx_rec = recession_period(Q)
    idx_oth = np.full(Q.shape[0], True)
    idx_oth[idx_rec] = False
    return param_calibrate_jit(param_range, method, Q, b_LH, a, idx_rec, idx_oth)


@njit(parallel=True)
def param_calibrate_jit(param_range, method, Q, b_LH, a, idx_rec, idx_oth):
    """
    Calibrates the parameters for a baseflow estimation method using the Numba-accelerated `param_calibrate_jit` function.
    
    The function takes in the range of parameter values to test, the baseflow estimation method, the discharge values, the low-flow baseflow values, and the parameter for the baseflow estimation method. It then calculates the recession period indices and other indices, and uses the `param_calibrate_jit` function to find the optimal parameter value from the given range.
    
    Args:
        param_range (numpy.ndarray): The range of parameter values to test.
        method (callable): The baseflow estimation method to use.
        Q (numpy.ndarray): The discharge values.
        b_LH (numpy.ndarray): The low-flow baseflow values.
        a (float): The parameter for the baseflow estimation method.
    
    Returns:
        float: The optimal parameter value from the given range.
    """
    logQ = np.log1p(Q)
    loss = np.zeros(param_range.shape)
    for i in prange(param_range.shape[0]):
        p = param_range[i]
        b_exceed = method(Q, b_LH, a, p, return_exceed=True)
        f_exd, logb = b_exceed[-1] / Q.shape[0], np.log1p(b_exceed[:-1])

        # NSE for recession part
        Q_obs, Q_sim = logQ[idx_rec], logb[idx_rec]
        SS_res = np.sum(np.square(Q_obs - Q_sim))
        SS_tot = np.sum(np.square(Q_obs - np.mean(Q_obs)))
        NSE_rec = (1 - SS_res / (SS_tot + 1e-10)) - 1e-10

        # NSE for other part
        Q_obs, Q_sim = logQ[idx_oth], logb[idx_oth]
        SS_res = np.sum(np.square(Q_obs - Q_sim))
        SS_tot = np.sum(np.square(Q_obs - np.mean(Q_obs)))
        NSE_oth = (1 - SS_res / (SS_tot + 1e-10)) - 1e-10

        loss[i] = 1 - (1 - (1 - NSE_rec) / (1 - NSE_oth)) * (1 - f_exd)
    return param_range[np.argmin(loss)]


def recession_period(Q):
    """
    Identifies the recession periods in the discharge time series.
    
    The function takes the discharge time series `Q` as input and returns the indices of the beginning and end of the recession periods. The recession periods are identified by finding the local maxima in the 3-point moving average of the discharge time series. The function keeps only the recession periods that are at least 10 time steps long, and trims the beginning of each recession period by 60% of the recession period duration.
    
    Args:
        Q (numpy.ndarray): The discharge time series.
    
    Returns:
        numpy.ndarray: The indices of the beginning and end of the recession periods.
    """
    idx_dec = np.zeros(Q.shape[0] - 1, dtype=np.int64)
    Q_ave = moving_average(Q, 3)
    idx_dec[1:-1] = (Q_ave[:-1] - Q_ave[1:]) > 0
    idx_beg = np.where(idx_dec[:-1] - idx_dec[1:] == -1)[0] + 1
    idx_end = np.where(idx_dec[:-1] - idx_dec[1:] == 1)[0] + 1
    idx_keep = (idx_end - idx_beg) >= 10
    idx_beg = idx_beg[idx_keep]
    idx_end = idx_end[idx_keep]
    duration = idx_end - idx_beg
    idx_beg = idx_beg + np.ceil(duration * 0.6).astype(np.int64)
    return multi_arange(idx_beg, idx_end)


def maxmium_BFI(Q, b_LH, a, date=None):
    """
    Calculates the maximum baseflow index (BFI) for a given discharge time series.
    
    The function takes the discharge time series `Q`, the baseflow time series `b_LH`, and the recession coefficient `a` as input. It calculates the annual baseflow and discharge, and then computes the maximum BFI. If the maximum BFI is greater than 0.9, the function returns the ratio of the total baseflow to the total discharge instead.
    
    Args:
        Q (numpy.ndarray): The discharge time series.
        b_LH (numpy.ndarray): The baseflow time series.
        a (float): The recession coefficient.
        date (datetime.datetime, optional): The date associated with the discharge time series. If provided, the function will compute the annual BFI for each year.
    
    Returns:
        float: The maximum baseflow index.
    """
    b = backward(Q, b_LH, a)

    if date is None:
        idx_end = b.shape[0] // 365 * 365
        annual_b = np.mean(b[:idx_end].reshape(-1, 365), axis=1)
        annual_Q = np.mean(Q[:idx_end].reshape(-1, 365), axis=1)
        annual_BFI = annual_b / annual_Q
    else:
        idx_year = date.year - date.year.min()
        counts = np.bincount(idx_year)
        idx_valid = counts > 0
        annual_b = np.bincount(idx_year, weights=b)[idx_valid] / counts[idx_valid]
        annual_Q = np.bincount(idx_year, weights=Q)[idx_valid] / counts[idx_valid]
        annual_BFI = annual_b / annual_Q

    BFI_max = np.max(annual_BFI)
    BFI_max = BFI_max if BFI_max < 0.9 else np.sum(annual_b) / np.sum(annual_Q)
    return BFI_max


def bflow_recession_analysis(Q, b, ndmin=10, ndmax=300):
    """BFlow recession analysis for SWAT alpha factor estimation.

    Identifies baseflow recession segments, assembles a master recession
    curve, and fits the exponential recession constant (alpha factor).
    This maps directly to SWAT's ALPHA_BF parameter.

    Arnold, J.G. and Allen, P.M. (1999) Automated methods for estimating
    baseflow and ground water recharge from streamflow records. JAWRA
    35(2), 411-424.

    Args:
        Q (numpy.ndarray): Total streamflow.
        b (numpy.ndarray): Separated baseflow (e.g., from lh_multi).
        ndmin (int): Minimum recession segment length (days). Default 10.
        ndmax (int): Maximum recession segment length (days). Default 300.

    Returns:
        dict: Dictionary with keys:
            - 'alpha_factor': exponential recession constant
            - 'baseflow_days': days for recession to decline one log cycle
            - 'n_segments': number of recession segments used
            - 'BFI': baseflow fraction (sum(b)/sum(Q))
    """
    n = len(b)

    # Identify recession segments: periods where baseflow is declining
    declining = np.zeros(n, dtype=bool)
    declining[1:] = b[1:] < b[:-1]

    # Find contiguous declining segments
    segments = []
    seg_start = None
    for i in range(n):
        if declining[i]:
            if seg_start is None:
                seg_start = i - 1  # include the point before decline starts
        else:
            if seg_start is not None:
                seg_len = i - seg_start
                if ndmin <= seg_len <= ndmax:
                    segments.append((seg_start, i))
                seg_start = None
    # Handle segment ending at array end
    if seg_start is not None:
        seg_len = n - seg_start
        if ndmin <= seg_len <= ndmax:
            segments.append((seg_start, n))

    if len(segments) == 0:
        return {
            'alpha_factor': np.nan,
            'baseflow_days': np.nan,
            'n_segments': 0,
            'BFI': np.sum(b) / np.sum(Q) if np.sum(Q) > 0 else np.nan,
        }

    # Fit exponential recession to each segment and compute alpha
    # For each segment: b(t) = b(0) * exp(-alpha * t)
    # Taking log: ln(b(t)) = ln(b(0)) - alpha * t
    # Fit slope via least squares across all segments (master recession curve)
    all_t = []
    all_logb = []
    for start, end in segments:
        seg_b = b[start:end]
        valid = seg_b > 0
        if np.sum(valid) < 2:
            continue
        t_rel = np.arange(end - start)[valid].astype(float)
        logb = np.log(seg_b[valid])
        # Normalize: shift t so each segment starts at t=0
        all_t.append(t_rel)
        # Normalize: shift logb so each segment starts at logb=0
        all_logb.append(logb - logb[0])

    if len(all_t) == 0:
        return {
            'alpha_factor': np.nan,
            'baseflow_days': np.nan,
            'n_segments': len(segments),
            'BFI': np.sum(b) / np.sum(Q) if np.sum(Q) > 0 else np.nan,
        }

    t_all = np.concatenate(all_t)
    logb_all = np.concatenate(all_logb)

    # Least squares fit: logb = -alpha * t  (no intercept since normalized)
    # alpha = -sum(t * logb) / sum(t^2)
    alpha = -np.sum(t_all * logb_all) / (np.sum(t_all ** 2) + 1e-30)
    alpha = max(alpha, 1e-10)  # ensure positive

    # Baseflow days: days for one log cycle decline
    # One log cycle = ln(10), so alpha * BFD = ln(10)
    baseflow_days = np.log(10) / alpha

    BFI = np.sum(b) / np.sum(Q) if np.sum(Q) > 0 else np.nan

    return {
        'alpha_factor': alpha,
        'baseflow_days': baseflow_days,
        'n_segments': len(segments),
        'BFI': BFI,
    }


def bflow(Q, beta=0.925):
    """BFlow automated baseflow separation and recession analysis.

    Runs the Lyne-Hollick filter with 3 passes (the BFlow protocol)
    and performs recession analysis to compute the SWAT alpha factor.

    Arnold, J.G. and Allen, P.M. (1999) Automated methods for estimating
    baseflow and ground water recharge from streamflow records. JAWRA
    35(2), 411-424.

    Args:
        Q (numpy.ndarray): Daily streamflow.
        beta (float): Filter parameter. Default 0.925.

    Returns:
        dict: Dictionary with keys:
            - 'baseflow': numpy array of baseflow values (3-pass result)
            - 'alpha_factor': exponential recession constant (SWAT ALPHA_BF)
            - 'baseflow_days': days for recession to decline one log cycle
            - 'n_segments': number of recession segments used
            - 'BFI': baseflow fraction from the 3-pass filter
            - 'BFI_pass1': baseflow fraction from pass 1
            - 'BFI_pass2': baseflow fraction from pass 2
    """
    from pybaseflow.separation import lh_multi

    total_Q = np.sum(Q)

    # Run 3 passes
    b1 = lh_multi(Q, beta=beta, num_pass=1)
    b2 = lh_multi(Q, beta=beta, num_pass=2)
    b3 = lh_multi(Q, beta=beta, num_pass=3)

    # Recession analysis on the 3-pass result
    recession = bflow_recession_analysis(Q, b3)

    return {
        'baseflow': b3,
        'alpha_factor': recession['alpha_factor'],
        'baseflow_days': recession['baseflow_days'],
        'n_segments': recession['n_segments'],
        'BFI': np.sum(b3) / total_Q if total_Q > 0 else np.nan,
        'BFI_pass1': np.sum(b1) / total_Q if total_Q > 0 else np.nan,
        'BFI_pass2': np.sum(b2) / total_Q if total_Q > 0 else np.nan,
    }


