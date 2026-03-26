"""Tracer-based baseflow separation methods."""

import numpy as np


def estimate_endmembers(SC, bf_percentile=99, ro_percentile=1):
    """Estimate baseflow and runoff end-member conductivities.

    Args:
        SC (numpy.ndarray): Specific conductance time series (uS/cm).
        bf_percentile (float): Percentile for baseflow end-member. Default 99.
        ro_percentile (float): Percentile for runoff end-member. Default 1.

    Returns:
        tuple: (SC_BF, SC_RO) estimated end-member conductivities.
    """
    SC_valid = SC[np.isfinite(SC)]
    SC_BF = np.percentile(SC_valid, bf_percentile)
    SC_RO = np.percentile(SC_valid, ro_percentile)
    return SC_BF, SC_RO


def cmb(Q, SC, SC_BF=None, SC_RO=None, sc_bf_percentile=99,
        sc_ro_percentile=1):
    """Conductivity Mass Balance baseflow separation (Stewart et al., 2007).

    Two-component mixing model using specific conductance as a tracer:
        Q_b(t) = Q(t) * (SC(t) - SC_RO) / (SC_BF - SC_RO)

    If end-member conductivities are not provided, they are estimated from
    the SC record using percentiles.

    Stewart, M.T., Cimino, J. and Ross, M. (2007) Calibration of base flow
    separation methods with streamflow conductivity. Groundwater 45(1), 17-27.

    Args:
        Q (numpy.ndarray): Streamflow.
        SC (numpy.ndarray): Concurrent specific conductance (uS/cm).
        SC_BF (float, optional): Baseflow end-member conductivity. If None,
            estimated as the sc_bf_percentile of SC.
        SC_RO (float, optional): Surface runoff end-member conductivity.
            If None, estimated as the sc_ro_percentile of SC.
        sc_bf_percentile (float): Percentile for SC_BF estimation. Default 99.
        sc_ro_percentile (float): Percentile for SC_RO estimation. Default 1.

    Returns:
        numpy.ndarray: Daily baseflow array.
    """
    if SC_BF is None or SC_RO is None:
        est_BF, est_RO = estimate_endmembers(
            SC, bf_percentile=sc_bf_percentile,
            ro_percentile=sc_ro_percentile)
        if SC_BF is None:
            SC_BF = est_BF
        if SC_RO is None:
            SC_RO = est_RO

    if SC_BF == SC_RO:
        raise ValueError("SC_BF and SC_RO are equal; cannot separate.")

    # Core CMB equation
    b = Q * (SC - SC_RO) / (SC_BF - SC_RO)

    # Physical constraints
    b = np.clip(b, 0.0, Q)

    # Propagate NaN from SC
    b[~np.isfinite(SC)] = np.nan

    return b


def calibrate_eckhardt_from_cmb(Q, SC, a=None, SC_BF=None, SC_RO=None,
                                sc_bf_percentile=99, sc_ro_percentile=1):
    """Calibrate Eckhardt BFImax using CMB as a reference.

    Runs CMB to get a reference baseflow, computes BFI, and returns a
    BFImax suitable for use with eckhardt(). Optionally estimates the
    recession coefficient if not provided.

    Args:
        Q (numpy.ndarray): Streamflow.
        SC (numpy.ndarray): Concurrent specific conductance (uS/cm).
        a (float, optional): Recession coefficient. If None, estimated
            from the streamflow using strict_baseflow + recession_coefficient.
        SC_BF (float, optional): Baseflow end-member conductivity.
        SC_RO (float, optional): Surface runoff end-member conductivity.
        sc_bf_percentile (float): Percentile for SC_BF estimation. Default 99.
        sc_ro_percentile (float): Percentile for SC_RO estimation. Default 1.

    Returns:
        dict: Dictionary with keys:
            - 'BFImax': calibrated maximum baseflow index
            - 'BFI_cmb': BFI computed from CMB
            - 'a': recession coefficient used
            - 'SC_BF': baseflow end-member conductivity used
            - 'SC_RO': runoff end-member conductivity used
    """
    from pybaseflow.separation import strict_baseflow
    from pybaseflow.estimate import recession_coefficient

    # Get CMB baseflow
    b_cmb = cmb(Q, SC, SC_BF=SC_BF, SC_RO=SC_RO,
                sc_bf_percentile=sc_bf_percentile,
                sc_ro_percentile=sc_ro_percentile)

    # Use only valid (non-NaN) pairs
    valid = np.isfinite(b_cmb) & np.isfinite(Q) & (Q > 0)
    BFI_cmb = np.sum(b_cmb[valid]) / np.sum(Q[valid])
    BFI_cmb = np.clip(BFI_cmb, 0.01, 0.99)

    # Estimate recession coefficient if not provided
    if a is None:
        strict = strict_baseflow(Q)
        a = recession_coefficient(Q, strict)

    # Resolve end-members used
    if SC_BF is None or SC_RO is None:
        est_BF, est_RO = estimate_endmembers(
            SC, bf_percentile=sc_bf_percentile,
            ro_percentile=sc_ro_percentile)
        SC_BF = SC_BF if SC_BF is not None else est_BF
        SC_RO = SC_RO if SC_RO is not None else est_RO

    return {
        'BFImax': BFI_cmb,
        'BFI_cmb': BFI_cmb,
        'a': a,
        'SC_BF': SC_BF,
        'SC_RO': SC_RO,
    }
