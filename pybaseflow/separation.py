import numpy as np
from numba import njit, prange


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _init_baseflow(Q, initial_method='Q0'):
    """Compute the initial baseflow value b[0].

    Args:
        Q (numpy.ndarray): Streamflow array.
        initial_method (str or float): 'Q0' (default), 'min', 'LH', or a float.

    Returns:
        float: The initial baseflow value.
    """
    if isinstance(initial_method, str):
        if initial_method == 'Q0':
            return Q[0]
        elif initial_method == 'min':
            return np.min(Q)
        elif initial_method == 'LH':
            return lh(Q)[0]
        else:
            raise ValueError(f"Invalid initial_method: {initial_method}")
    return float(initial_method)


def _recursive_digital_filter(Q, alpha, beta, gamma=0.0,
                              initial_method='Q0', return_exceed=False):
    """Generalized recursive digital filter for baseflow separation.

    Computes baseflow using the general form:
        b[t] = alpha * b[t-1] + beta * (Q[t] + gamma * Q[t-1])
    with the constraint b[t] <= Q[t].

    Two structural families of filters map onto this form:
      - gamma = 0: Chapman-Maxwell, Boughton, Eckhardt, EWMA (linear reservoir)
      - gamma = 1: Lyne-Hollick, Chapman (1991), Willems (signal processing)

    Args:
        Q (numpy.ndarray): Streamflow array.
        alpha (float): Coefficient on previous baseflow b[t-1].
        beta (float): Coefficient on streamflow terms.
        gamma (float): Weight on previous streamflow Q[t-1]. 0 or 1 for
            standard filters; other values for IHACRES-type filters.
        initial_method (str or float): Method to set b[0]. See _init_baseflow().
        return_exceed (bool): If True, the last element of the returned array
            holds the count of timesteps where computed baseflow exceeded Q.

    Returns:
        numpy.ndarray: Baseflow array (length Q.shape[0], or Q.shape[0]+1
            if return_exceed is True).
    """
    n = Q.shape[0]
    if return_exceed:
        b = np.zeros(n + 1)
    else:
        b = np.zeros(n)

    b[0] = _init_baseflow(Q, initial_method)

    for i in range(n - 1):
        b[i + 1] = alpha * b[i] + beta * (Q[i + 1] + gamma * Q[i])
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b


# ---------------------------------------------------------------------------
# Recursive digital filters — gamma = 0 family (linear reservoir based)
# ---------------------------------------------------------------------------

def boughton(Q, a, C, initial_method='Q0', return_exceed=False):
    """Boughton two-parameter filter (Boughton, 1993).

    b[t] = a/(1+C) * b[t-1] + C/(1+C) * Q[t]

    Boughton W.C. (1993) - A hydrograph-based model for estimating water yield
    of ungauged catchments. Inst. Eng. Aust. Natl. Conf. Publ. 93/14, 317-324.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient.
        C (float): Calibrated parameter (see param_calibrate).
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    return _recursive_digital_filter(
        Q, alpha=a / (1 + C), beta=C / (1 + C), gamma=0.0,
        initial_method=initial_method, return_exceed=return_exceed)


def chapman_maxwell(Q, a, initial_method='Q0', return_exceed=False):
    """Chapman-Maxwell filter (Chapman & Maxwell, 1996).

    b[t] = a/(2-a) * b[t-1] + (1-a)/(2-a) * Q[t]

    Equivalent to Eckhardt with BFImax = 0.5.

    Chapman, T.G., Maxwell, A.I. (1996) - Baseflow separation - comparison of
    numerical methods with tracer experiments. Hydrol. Water Resour. Symp.,
    Inst. Eng. Australia, Hobart, 539-545.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    return _recursive_digital_filter(
        Q, alpha=a / (2 - a), beta=(1 - a) / (2 - a), gamma=0.0,
        initial_method=initial_method, return_exceed=return_exceed)


def eckhardt(Q, a, BFImax, initial_method='Q0', return_exceed=False):
    """Eckhardt two-parameter filter (Eckhardt, 2005).

    b[t] = ((1-BFImax)*a*b[t-1] + (1-a)*BFImax*Q[t]) / (1 - a*BFImax)

    Eckhardt, K. (2005) How to construct recursive digital filters for
    baseflow separation. Hydrol. Process. 19, 507-515.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient.
        BFImax (float): Maximum baseflow index.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    denom = 1 - a * BFImax
    return _recursive_digital_filter(
        Q, alpha=(1 - BFImax) * a / denom,
        beta=(1 - a) * BFImax / denom, gamma=0.0,
        initial_method=initial_method, return_exceed=return_exceed)


def what(Q, BFImax, a):
    """WHAT method (Lim et al., 2005). Alias for eckhardt().

    Mathematically identical to the Eckhardt filter.

    Lim, K.J., et al. (2005) Automated Web GIS based hydrograph analysis tool,
    WHAT. JAWRA 41(6), 1407-1416.

    Args:
        Q (numpy.ndarray): Streamflow.
        BFImax (float): Maximum baseflow index.
        a (float): Recession coefficient.

    Returns:
        numpy.ndarray: Baseflow.
    """
    return eckhardt(Q, a, BFImax)


def ewma(Q, e, initial_method='Q0', return_exceed=False):
    """Exponential weighted moving average filter (Tularam & Ilahee, 2008).

    b[t] = (1-e) * b[t-1] + e * Q[t]

    Tularam, G.A. and Ilahee, M. (2008) Exponential smoothing method of base
    flow separation and its impact on continuous loss estimates. Am. J. Environ.
    Sci. 4(2), 136-144.

    Args:
        Q (numpy.ndarray): Streamflow.
        e (float): Smoothing parameter.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    return _recursive_digital_filter(
        Q, alpha=(1 - e), beta=e, gamma=0.0,
        initial_method=initial_method, return_exceed=return_exceed)


def furey(Q, a, A, initial_method='Q0', return_exceed=False):
    """Furey-Gupta physically-based filter (Furey & Gupta, 2001).

    b[t] = (a - A*(1-a)) * b[t-1] + A*(1-a) * Q[t-1]

    Note: this filter uses Q[t-1] (not Q[t]) in the second term, so gamma
    acts on the *previous* timestep's streamflow via the general form with
    beta applied to Q[t] + gamma*Q[t-1]. Here we set beta on Q[t-1] by
    mapping: alpha = a - A*(1-a), beta = A*(1-a), gamma = 0, and shifting
    the Q index.

    Furey, P.R. and Gupta, V.K. (2001) A physically based filter for
    separating base flow from streamflow time series. Water Resour. Res.
    37(11), 2709-2722.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient.
        A (float): Calibrated parameter.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    # Furey uses Q[i] (previous timestep) rather than Q[i+1], so it doesn't
    # fit the standard general form directly. Implement the loop explicitly.
    n = Q.shape[0]
    if return_exceed:
        b = np.zeros(n + 1)
    else:
        b = np.zeros(n)

    b[0] = _init_baseflow(Q, initial_method)

    for i in range(n - 1):
        b[i + 1] = (a - A * (1 - a)) * b[i] + A * (1 - a) * Q[i]
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1
    return b


# ---------------------------------------------------------------------------
# Recursive digital filters — gamma = 1 family (signal processing based)
# ---------------------------------------------------------------------------

def chapman(Q, a=0.925, initial_method='Q0', return_exceed=False):
    """Chapman filter (Chapman, 1991).

    b[t] = (3a-1)/(3-a) * b[t-1] + (1-a)/(3-a) * (Q[t] + Q[t-1])

    Chapman, T.G. (1991) Comment on 'Evaluation of automated techniques for
    base flow and recession analyses' by R.J. Nathan and T.A. McMahon.
    Water Resour. Res. 27(7), 1783-1784.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient. Default 0.925.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    return _recursive_digital_filter(
        Q, alpha=(3 * a - 1) / (3 - a), beta=(1 - a) / (3 - a), gamma=1.0,
        initial_method=initial_method, return_exceed=return_exceed)


def willems(Q, a, w, initial_method='Q0', return_exceed=False):
    """Willems digital filter (Willems, 2009).

    v = (1-w)*(1-a) / (2*w)
    b[t] = (a-v)/(1+v) * b[t-1] + v/(1+v) * (Q[t] + Q[t-1])

    Willems, P. (2009) A time series tool to support the multi-criteria
    performance evaluation of rainfall-runoff models. Environ. Model. Softw.
    24(3), 311-321.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient.
        w (float): Average proportion of quick flow in streamflow.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    v = (1 - w) * (1 - a) / (2 * w)
    return _recursive_digital_filter(
        Q, alpha=(a - v) / (1 + v), beta=v / (1 + v), gamma=1.0,
        initial_method=initial_method, return_exceed=return_exceed)


# ---------------------------------------------------------------------------
# Recursive digital filters — variable gamma (hybrid)
# ---------------------------------------------------------------------------

def ihacres(Q, a, C, alpha_s, initial_method='Q0', return_exceed=False):
    """Jakeman-Hornberger / IHACRES baseflow filter (Jakeman & Hornberger, 1993).

    b[t] = a/(1+C) * b[t-1] + C/(1+C) * (Q[t] + alpha_s * Q[t-1])

    Three-parameter extension of the Boughton filter. When alpha_s = 0,
    reduces exactly to boughton(). The alpha_s term adds dependence on
    the previous timestep's total flow, accounting for delayed baseflow
    response.

    In the generalized filter taxonomy: alpha = a/(1+C), beta = C/(1+C),
    gamma = alpha_s. This filter occupies the continuum between the
    gamma=0 and gamma=1 families.

    Jakeman, A.J. and Hornberger, G.M. (1993) How much complexity is
    warranted in a rainfall-runoff model? Water Resour. Res. 29(8),
    2637-2649.

    Args:
        Q (numpy.ndarray): Streamflow.
        a (float): Recession coefficient for slow flow (0.9-0.999).
        C (float): Partitioning parameter (0.1-1.0).
        alpha_s (float): Previous-timestep dependence (-0.99 to 0).
            Set to 0 to recover the Boughton filter.
        initial_method (str or float): Method to set b[0]. Default 'Q0'.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    return _recursive_digital_filter(
        Q, alpha=a / (1 + C), beta=C / (1 + C), gamma=alpha_s,
        initial_method=initial_method, return_exceed=return_exceed)


# ---------------------------------------------------------------------------
# Lyne-Hollick filter (single-pass and multi-pass)
# ---------------------------------------------------------------------------

def lh(Q, beta=0.925, return_exceed=False):
    """Lyne-Hollick digital filter, 2-pass (Lyne & Hollick, 1979).

    Forward pass:  b[t] = beta*b[t-1] + (1-beta)/2 * (Q[t] + Q[t-1])
    Backward pass: applied to the forward result for smoothing.

    The default beta=0.925 follows the recommendation of Nathan & McMahon
    (1990). Using lh_multi(Q, beta=0.925, num_pass=3) reproduces the
    Nathan-McMahon recommended 3-pass protocol exactly.

    Lyne, V. and Hollick, M. (1979) Stochastic time-variable rainfall-runoff
    modelling. Inst. Eng. Aust. Natl. Conf., 89-93.

    Nathan, R.J. and McMahon, T.A. (1990) Evaluation of automated techniques
    for base flow and recession analyses. Water Resour. Res. 26(7), 1465-1473.

    Args:
        Q (numpy.ndarray): Streamflow.
        beta (float): Filter parameter. Default 0.925.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    if return_exceed:
        b = np.zeros(Q.shape[0] + 1)
    else:
        b = np.zeros(Q.shape[0])

    # first pass (forward)
    b[0] = Q[0]
    for i in range(Q.shape[0] - 1):
        b[i + 1] = beta * b[i] + (1 - beta) / 2 * (Q[i] + Q[i + 1])
        if b[i + 1] > Q[i + 1]:
            b[i + 1] = Q[i + 1]
            if return_exceed:
                b[-1] += 1

    # second pass (backward)
    b1 = np.copy(b)
    for i in range(Q.shape[0] - 2, -1, -1):
        b[i] = beta * b[i + 1] + (1 - beta) / 2 * (b1[i + 1] + b1[i])
        if b[i] > b1[i]:
            b[i] = b1[i]
            if return_exceed:
                b[-1] += 1
    return b


def lh_multi(Q, beta=0.925, num_pass=2, return_exceed=False):
    """Lyne-Hollick digital filter, n-pass (Spongberg, 2000).

    Applies the LH filter with alternating forward/backward passes.

    Spongberg, M.E. (2000) Spectral analysis of base flow separation with
    digital filters. Water Resour. Res. 36(3), 745-752.

    Args:
        Q (numpy.ndarray): Streamflow.
        beta (float): Filter parameter. Default 0.925.
        num_pass (int): Number of passes. Default 2.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
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


# ---------------------------------------------------------------------------
# Graphical / interval-based methods
# ---------------------------------------------------------------------------

def _hysep_interval(area):
    """Compute the HYSEP separation interval from drainage area.

    N = A^0.2 where A is in square miles. The interval 2N* is the odd
    integer between 3 and 11 nearest to 2N.

    Args:
        area (float or None): Basin area in km^2. If None, defaults to N=5.

    Returns:
        int: The HYSEP interval (odd integer, 3-11).
    """
    if area is None:
        N = 5
    else:
        N = np.power(0.3861022 * area, 0.2)
    inN = np.ceil(2 * N)
    if np.mod(inN, 2) == 0:
        inN = np.ceil(2 * N) - 1
    inN = np.int64(min(max(inN, 3), 11))
    return inN


# Keep the old name as an alias for backward compatibility
hysep_interval = _hysep_interval


def fixed(Q, area=None):
    """Fixed interval graphical method from HYSEP (Sloto & Crouse, 1996).

    Sloto, R.A. & Crouse, M.Y. (1996) HYSEP: A computer program for
    streamflow hydrograph separation and analysis. USGS WRI 96-4040.

    Args:
        Q (numpy.ndarray): Streamflow.
        area (float): Basin area in km^2.

    Returns:
        numpy.ndarray: Baseflow.
    """
    inN = _hysep_interval(area)
    return _fixed_interpolation(Q, inN)


@njit
def _fixed_interpolation(Q, inN):
    b = np.zeros(Q.shape[0])
    n = Q.shape[0] // inN
    for i in prange(n):
        b[inN * i:inN * (i + 1)] = np.min(Q[inN * i:inN * (i + 1)])
    if n * inN != Q.shape[0]:
        b[n * inN:] = np.min(Q[n * inN:])
    return b


def slide(Q, area):
    """Sliding interval graphical method from HYSEP (Sloto & Crouse, 1996).

    Sloto, R.A. & Crouse, M.Y. (1996) HYSEP: A computer program for
    streamflow hydrograph separation and analysis. USGS WRI 96-4040.

    Args:
        Q (numpy.ndarray): Streamflow.
        area (float): Basin area in km^2.

    Returns:
        numpy.ndarray: Baseflow.
    """
    inN = _hysep_interval(area)
    return _slide_interpolation(Q, inN)


def _slide_interpolation(Q, inN):
    b = np.zeros(Q.shape[0])
    for i in prange(np.int64((inN - 1) / 2), np.int64(Q.shape[0] - (inN - 1) / 2)):
        b[i] = np.min(Q[np.int64(i - (inN - 1) / 2):np.int64(i + (inN + 1) / 2)])
    b[:np.int64((inN - 1) / 2)] = np.min(Q[:np.int64((inN - 1) / 2)])
    b[np.int64(Q.shape[0] - (inN - 1) / 2):] = np.min(
        Q[np.int64(Q.shape[0] - (inN - 1) / 2):])
    return b


def local(Q, b_LH, area=None, return_exceed=False):
    """Local minimum graphical method from HYSEP (Sloto & Crouse, 1996).

    Sloto, R.A. & Crouse, M.Y. (1996) HYSEP: A computer program for
    streamflow hydrograph separation and analysis. USGS WRI 96-4040.

    Args:
        Q (numpy.ndarray): Streamflow.
        b_LH (numpy.ndarray): Lyne-Hollick baseflow (used for edge filling).
        area (float): Basin area in km^2.
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    idx_turn = _local_turn(Q, _hysep_interval(area))
    if idx_turn.shape[0] < 3:
        raise IndexError('Less than 3 turning points found')
    b = _linear_interpolation(Q, idx_turn, return_exceed=return_exceed)
    b[:idx_turn[0]] = b_LH[:idx_turn[0]]
    b[idx_turn[-1] + 1:] = b_LH[idx_turn[-1] + 1:]
    return b


@njit
def _local_turn(Q, inN):
    idx_turn = np.zeros(Q.shape[0], dtype=np.int64)
    for i in prange(np.int64((inN - 1) / 2), np.int64(Q.shape[0] - (inN - 1) / 2)):
        if Q[i] == np.min(Q[np.int64(i - (inN - 1) / 2):np.int64(i + (inN + 1) / 2)]):
            idx_turn[i] = i
    return idx_turn[idx_turn != 0]


def ukih(Q, b_LH, return_exceed=False):
    """UK Institute of Hydrology smoothed minima method (UKIH, 1980).

    Aksoy, H., Kurt, I. and Eris, E. (2009) Filtered smoothed minima baseflow
    separation method. J. Hydrol. 372(1), 94-101.

    Args:
        Q (numpy.ndarray): Streamflow.
        b_LH (numpy.ndarray): Lyne-Hollick baseflow (used for edge filling).
        return_exceed (bool): If True, append exceed count.

    Returns:
        numpy.ndarray: Baseflow.
    """
    N = 5
    block_end = Q.shape[0] // N * N
    idx_min = np.argmin(Q[:block_end].reshape(-1, N), axis=1)
    idx_min = idx_min + np.arange(0, block_end, N)
    idx_turn = _ukih_turn(Q, idx_min)
    if idx_turn.shape[0] < 3:
        raise IndexError('Less than 3 turning points found')
    b = _linear_interpolation(Q, idx_turn, return_exceed=return_exceed)
    b[:idx_turn[0]] = b_LH[:idx_turn[0]]
    b[idx_turn[-1] + 1:] = b_LH[idx_turn[-1] + 1:]
    return b


@njit
def _ukih_turn(Q, idx_min):
    idx_turn = np.zeros(idx_min.shape[0], dtype=np.int64)
    for i in prange(idx_min.shape[0] - 2):
        if ((0.9 * Q[idx_min[i + 1]] < Q[idx_min[i]]) &
                (0.9 * Q[idx_min[i + 1]] < Q[idx_min[i + 2]])):
            idx_turn[i] = idx_min[i + 1]
    return idx_turn[idx_turn != 0]


@njit
def _linear_interpolation(Q, idx_turn, return_exceed=False):
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


# ---------------------------------------------------------------------------
# PART method (Rutledge, 1998)
# ---------------------------------------------------------------------------

def part(Q, area, log_cycle_threshold=0.1):
    """PART baseflow separation (Rutledge, 1998).

    Recession-based graphical method that identifies days where streamflow
    equals baseflow using an antecedent recession requirement, then
    interpolates in log-space between qualifying days.

    The procedure is run three times with different recession durations
    (floor(N), ceil(N), ceil(N)+1) and combined via linear interpolation
    evaluated at the exact N = A^0.2.

    Rutledge, A.T. (1998) Computer programs for describing the recession of
    ground-water discharge and for estimating mean ground-water recharge and
    discharge from streamflow records — Update. USGS WRI 98-4148.

    Args:
        Q (numpy.ndarray): Daily streamflow.
        area (float): Drainage area in km^2.
        log_cycle_threshold (float): Maximum daily log10 decline for a
            qualifying day. Default 0.1.

    Returns:
        numpy.ndarray: Daily baseflow array.
    """
    # Compute real-valued N (days of antecedent recession)
    area_sqmi = 0.3861022 * area
    N_real = area_sqmi ** 0.2

    N_floor = max(int(np.floor(N_real)), 1)
    N_ceil = max(int(np.ceil(N_real)), 2)
    if N_ceil == N_floor:
        N_ceil = N_floor + 1
    N_values = [N_floor, N_ceil, N_ceil + 1]

    results = []
    for N in N_values:
        b = _part_single(Q, N, log_cycle_threshold)
        results.append(b)

    # Linear interpolation between floor and ceil results
    if N_ceil > N_floor:
        frac = (N_real - N_floor) / (N_ceil - N_floor)
    else:
        frac = 0.0
    b_out = results[0] * (1 - frac) + results[1] * frac

    return b_out


def _part_single(Q, N, log_cycle_threshold):
    """Run the PART procedure for a single value of N."""
    n = len(Q)
    eps = 1e-99  # avoid log(0)
    Q_safe = np.where(Q > 0, Q, eps)
    logQ = np.log10(Q_safe)

    # Step 1: identify qualifying days (antecedent recession for N days)
    is_gw = np.zeros(n, dtype=bool)
    for t in range(N, n):
        recession = True
        for j in range(t - N, t):
            if Q_safe[j] < Q_safe[j + 1]:
                recession = False
                break
        if recession:
            is_gw[t] = True

    # Step 2: apply log-cycle check — disqualify if followed by steep decline
    for t in range(n - 1):
        if is_gw[t]:
            if logQ[t] - logQ[t + 1] > log_cycle_threshold:
                is_gw[t] = False

    # Step 3: iterative log-linear interpolation with correction
    b = np.zeros(n)
    b = _part_interpolate(Q_safe, logQ, is_gw, b)

    # Cap small values to zero
    b[b < 1e-6] = 0.0
    return b


def _part_interpolate(Q, logQ, is_gw, b):
    """Log-linear interpolation with iterative correction."""
    n = len(Q)
    is_anchor = is_gw.copy()

    for _iteration in range(n):  # convergence guaranteed
        # Set baseflow at anchor points
        b[:] = 0.0
        anchor_indices = np.where(is_anchor)[0]

        if len(anchor_indices) == 0:
            b[:] = Q[:]
            return b

        # Constant extrapolation at edges
        if anchor_indices[0] > 0:
            b[:anchor_indices[0]] = Q[anchor_indices[0]]
        if anchor_indices[-1] < n - 1:
            b[anchor_indices[-1] + 1:] = Q[anchor_indices[-1]]

        # Set anchor values
        for idx in anchor_indices:
            b[idx] = Q[idx]

        # Log-linear interpolation between anchors
        for k in range(len(anchor_indices) - 1):
            i0 = anchor_indices[k]
            i1 = anchor_indices[k + 1]
            if i1 - i0 <= 1:
                continue
            log_b0 = logQ[i0]
            log_b1 = logQ[i1]
            for t in range(i0 + 1, i1):
                frac = (t - i0) / (i1 - i0)
                b[t] = 10.0 ** (log_b0 + frac * (log_b1 - log_b0))

        # Check for violations (baseflow > streamflow)
        violations = np.where((b > Q) & ~is_anchor)[0]
        if len(violations) == 0:
            break

        # Find worst violation in each inter-anchor interval
        new_anchors = []
        for k in range(len(anchor_indices) - 1):
            i0 = anchor_indices[k]
            i1 = anchor_indices[k + 1]
            mask = (violations >= i0) & (violations <= i1)
            interval_violations = violations[mask]
            if len(interval_violations) > 0:
                log_excess = np.log10(np.maximum(b[interval_violations], 1e-99)) - \
                             logQ[interval_violations]
                worst = interval_violations[np.argmax(log_excess)]
                new_anchors.append(worst)

        if not new_anchors:
            break

        for idx in new_anchors:
            is_anchor[idx] = True

    # Final cap: baseflow cannot exceed streamflow
    b = np.minimum(b, Q)
    return b


# ---------------------------------------------------------------------------
# Strict baseflow identification
# ---------------------------------------------------------------------------

def strict_baseflow(Q, ice=None, quantile=0.9):
    """Identify strict baseflow periods in a streamflow time series.

    Applies heuristic rules based on the derivative of the hydrograph to
    identify timesteps that are unambiguously baseflow-dominated. Used by
    recession_coefficient() for recession constant estimation.

    Args:
        Q (numpy.ndarray): Streamflow.
        ice (numpy.ndarray, optional): Boolean mask for ice-affected periods.
        quantile (float): Quantile threshold for major events. Default 0.9.

    Returns:
        numpy.ndarray: Boolean mask (True = strict baseflow).
    """
    dQ = (Q[2:] - Q[:-2]) / 2

    # 1. flow data associated with positive and zero values of dy / dt
    wet1 = np.concatenate([[True], dQ >= 0, [True]])

    # 2. previous 2 points before points with dy/dt>=0, plus the next 3 points
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

    # dry points = strict baseflow
    dry = ~(wet1 + wet2 + wet3 + wet4)

    if ice is not None:
        dry[ice] = False

    return dry


# ---------------------------------------------------------------------------
# Brutsaert-Nieber recession analysis (bn77)
# ---------------------------------------------------------------------------

def bn77(Q, L_min, snow_freeze_period, observational_precision, quantile=0.9):
    """Brutsaert-Nieber drought flow identification (Cheng et al., 2016).

    Identifies drought (pure baseflow) points in a discharge time series
    using recession slope analysis and multiple elimination criteria.

    Cheng, L., Zhang, L. and Brutsaert, W. (2016) Automated selection of pure
    base flows from regular daily streamflow data: objective algorithm.
    J. Hydrol. Eng. 21(11), 06016008.

    Args:
        Q (numpy.ndarray): Discharge time series.
        L_min (int): Minimum recession episode length.
        snow_freeze_period (tuple): (start_index, end_index) of snow/freeze period.
        observational_precision (float): Minimum Q threshold.
        quantile (float): Quantile for major event identification. Default 0.9.

    Returns:
        numpy.ndarray: Indices of drought flow points.
    """
    S = _estimate_recession_slope(Q)
    recession_episodes = _identify_recession_episodes(S, L_min)
    drought_flow_points = _eliminate_points(
        recession_episodes, L_min, snow_freeze_period,
        observational_precision, Q, S, quantile)
    return drought_flow_points


@njit
def _estimate_recession_slope(Q):
    """Estimate recession slope S(t) = (Q(t-1) - Q(t+1)) / 2."""
    N = len(Q)
    S = np.zeros(N)
    for i in range(1, N - 1):
        S[i] = (Q[i - 1] - Q[i + 1]) / 2
    return S


def _identify_recession_episodes(S, L_min):
    """Identify preliminary recession episodes where S > 0 for >= L_min steps."""
    i = 0
    N = len(S)
    recession_episodes = []

    while i < N - 1:
        if S[i] <= 0 and S[i + 1] > 0:
            episode_start = i + 1
            l = 0
            i = i + 1
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


def _eliminate_points(recession_episodes, L_min, snow_freeze_period,
                      observational_precision, Q, S, quantile):
    """Apply elimination criteria C3-C9 to recession episodes."""
    drought_flow_points = []
    major_event_threshold = np.quantile(Q, quantile)

    for episode in recession_episodes:
        # C4: Remove first 3 points if starting above quantile threshold
        if Q[episode[0]] > major_event_threshold:
            if len(episode) > 3:
                episode = episode[3:]
            else:
                continue
        else:
            # C3: Remove first 2 points otherwise
            if len(episode) > 2:
                episode = episode[2:]
            else:
                continue

        # C5: Remove the last point
        if len(episode) > 1:
            episode = episode[:-1]
        else:
            continue

        # C6: Remove points where S[i]/S[i-1] >= 2
        if len(episode) > 1:
            episode = episode[1:][S[episode[1:]] / S[episode[:-1]] < 2]

        # C7: Remove points where S[i] < S[i+1]
        if len(episode) > 1:
            episode = episode[:-1][S[episode[:-1]] >= S[episode[1:]]]

        # C8: Remove points during snow/freeze period
        episode = episode[(episode < snow_freeze_period[0]) |
                          (episode > snow_freeze_period[1])]

        # C9: Remove points where Q < observational precision
        episode = episode[Q[episode] >= observational_precision]

        if len(episode) > 0:
            drought_flow_points.append(episode)

    return np.concatenate(drought_flow_points) if drought_flow_points else np.array([], dtype=np.int64)
