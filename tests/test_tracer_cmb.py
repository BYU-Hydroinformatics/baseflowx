"""Implementation correctness checks for baseflowx.tracer (CMB).

This is the one method left with zero validation coverage after Task 1's
Xie(2020) cross-check and the three literature benchmarks -- Xie(2020)
doesn't implement CMB at all, and no per-gage published CMB result was
found with reproducible end-members and period of record. This module
closes that gap the same way PART's core logic was checked in
test_internal_consistency.py: not against an external reference, but
against a synthetic example where the correct answer is known by
construction, plus the physical/structural invariants the method is
documented to satisfy.

This validates that the *code* correctly implements the two-endmember
mixing equation (Stewart, Cimino & Ross, 2007). It says nothing about
whether CMB's answer is closer to the true baseflow for any real
catchment -- that is Task 3's question, not this one.
"""

import numpy as np
import pytest

from baseflowx import tracer


# ---------------------------------------------------------------------------
# Exact recovery on a synthetic mixing example with a known answer
# ---------------------------------------------------------------------------

def test_recovers_exact_baseflow_from_constructed_mixture():
    """Construct Q_baseflow and Q_runoff directly, so the true split is
    known by construction. Derive the SC that a perfect two-endmember
    mixture would produce, then confirm cmb() recovers the true baseflow
    split exactly from (Q, SC) alone -- the only inputs it is actually
    given.
    """
    rng = np.random.default_rng(0)
    n = 500
    SC_BF, SC_RO = 800.0, 60.0  # uS/cm; groundwater saltier than runoff

    Q_bf = 10.0 + 5.0 * np.sin(np.linspace(0, 6 * np.pi, n)) ** 2  # slow-varying
    Q_ro = rng.gamma(shape=1.5, scale=8.0, size=n)  # flashy, mostly near zero
    Q = Q_bf + Q_ro
    true_bfi = Q_bf / Q

    # Flow-weighted mixing: SC*Q = SC_BF*Q_bf + SC_RO*Q_ro
    SC = (SC_BF * Q_bf + SC_RO * Q_ro) / Q

    b = tracer.cmb(Q, SC, SC_BF=SC_BF, SC_RO=SC_RO)

    np.testing.assert_allclose(b, Q_bf, rtol=1e-10)
    np.testing.assert_allclose(b / Q, true_bfi, rtol=1e-10)


def test_pure_baseflow_and_pure_runoff_endpoints():
    """At the two endpoints of the mixing line -- SC exactly at the
    baseflow end-member, or exactly at the runoff end-member -- the split
    must be unambiguous: all baseflow, or none."""
    Q = np.array([10.0, 10.0])
    SC_BF, SC_RO = 800.0, 60.0
    SC = np.array([SC_BF, SC_RO])  # day 0: pure baseflow; day 1: pure runoff

    b = tracer.cmb(Q, SC, SC_BF=SC_BF, SC_RO=SC_RO)
    assert b[0] == pytest.approx(10.0)
    assert b[1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# estimate_endmembers(): hand-checkable percentile behavior
# ---------------------------------------------------------------------------

def test_estimate_endmembers_matches_hand_computed_percentiles():
    SC = np.arange(1, 101, dtype=float)  # 1..100, so percentiles are exact
    SC_BF, SC_RO = tracer.estimate_endmembers(SC, bf_percentile=99, ro_percentile=1)
    assert SC_BF == pytest.approx(np.percentile(SC, 99))
    assert SC_RO == pytest.approx(np.percentile(SC, 1))
    assert SC_BF > SC_RO


def test_estimate_endmembers_ignores_nan():
    SC = np.concatenate([np.arange(1, 101, dtype=float), [np.nan, np.nan]])
    SC_BF, SC_RO = tracer.estimate_endmembers(SC, bf_percentile=99, ro_percentile=1)
    assert np.isfinite(SC_BF) and np.isfinite(SC_RO)


# ---------------------------------------------------------------------------
# Physical and structural invariants
# ---------------------------------------------------------------------------

def test_cmb_never_exceeds_streamflow_or_goes_negative():
    rng = np.random.default_rng(1)
    n = 1000
    Q = rng.gamma(shape=2.0, scale=20.0, size=n) + 0.1
    # SC with occasional out-of-endmember-range noise, to exercise the clip
    SC = rng.uniform(-500, 2000, size=n)
    SC_BF, SC_RO = 800.0, 60.0

    b = tracer.cmb(Q, SC, SC_BF=SC_BF, SC_RO=SC_RO)
    assert np.all(b <= Q + 1e-9)
    assert np.all(b >= 0)


def test_cmb_propagates_nan_from_missing_sc():
    Q = np.array([10.0, 10.0, 10.0])
    SC = np.array([800.0, np.nan, 60.0])
    b = tracer.cmb(Q, SC, SC_BF=800.0, SC_RO=60.0)
    assert np.isnan(b[1])
    assert np.isfinite(b[0]) and np.isfinite(b[2])


def test_cmb_rejects_degenerate_endmembers():
    Q = np.array([10.0, 10.0])
    SC = np.array([500.0, 500.0])
    with pytest.raises(ValueError):
        tracer.cmb(Q, SC, SC_BF=500.0, SC_RO=500.0)


def test_cmb_defaults_to_estimated_endmembers_when_not_given():
    """When SC_BF/SC_RO aren't passed explicitly, cmb() should fall back to
    estimate_endmembers() -- confirm the two paths agree."""
    rng = np.random.default_rng(2)
    Q = rng.gamma(shape=2.0, scale=20.0, size=200) + 0.1
    SC = rng.uniform(50, 900, size=200)

    b_explicit_default = tracer.cmb(
        Q, SC,
        SC_BF=np.percentile(SC, 99), SC_RO=np.percentile(SC, 1))
    b_auto = tracer.cmb(Q, SC)  # should estimate the same 99th/1st percentiles
    np.testing.assert_allclose(b_explicit_default, b_auto, rtol=1e-10)


# ---------------------------------------------------------------------------
# calibrate_eckhardt_from_cmb(): end-to-end sanity
# ---------------------------------------------------------------------------

def test_calibrate_eckhardt_from_cmb_recovers_known_bfi():
    """On the same constructed mixture used above, BFI_cmb should match the
    true baseflow index, and BFImax should land in a valid range."""
    rng = np.random.default_rng(3)
    n = 500
    SC_BF, SC_RO = 800.0, 60.0

    Q_bf = 10.0 + 5.0 * np.sin(np.linspace(0, 6 * np.pi, n)) ** 2
    Q_ro = rng.gamma(shape=1.5, scale=8.0, size=n)
    Q = Q_bf + Q_ro
    SC = (SC_BF * Q_bf + SC_RO * Q_ro) / Q
    true_bfi = np.sum(Q_bf) / np.sum(Q)

    # Recession coefficient estimation needs a realistic hydrograph shape
    # (sustained recessions) and isn't what this test is exercising -- pass
    # a fixed value so the test isolates the CMB-derived BFImax/BFI_cmb path.
    result = tracer.calibrate_eckhardt_from_cmb(Q, SC, a=0.9, SC_BF=SC_BF, SC_RO=SC_RO)

    assert result["BFI_cmb"] == pytest.approx(true_bfi, rel=1e-9)
    assert 0.01 <= result["BFImax"] <= 0.99
    assert result["SC_BF"] == SC_BF and result["SC_RO"] == SC_RO
