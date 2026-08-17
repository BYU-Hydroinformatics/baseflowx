"""Validate baseflowx against an independent reference implementation.

The `baseflow` PyPI package (Xie et al., 2020 -- "Evaluation of typical
methods for baseflow separation in the contiguous United States", J.
Hydrol. 583, 124628) implements 12 of baseflowx's 17 methods. baseflowx
was originally refactored and extended from that codebase, so this is
not an independent derivation of the underlying equations -- but it is
a genuinely independent *implementation*: separately packaged, separately
maintained, running through Numba-JIT-compiled code rather than
baseflowx's plain NumPy loops. Agreement to floating-point precision
confirms the refactor did not silently change the recursion, the edge
handling, or the initial-condition logic for any of the twelve shared
methods, across six gages spanning distinct hydroclimatic settings.

Requires the optional `baseflow` package: `pip install baseflowx[validation]`.
Skips cleanly if it is not installed.
"""

import numpy as np
import pytest

xie = pytest.importorskip("baseflow", reason="pip install baseflowx[validation]")

from baseflowx import estimate, separation

RTOL = 1e-10
ATOL = 1e-10


def _metrics(a, b):
    diff = np.asarray(a) - np.asarray(b)
    return {
        "max_abs_dev": float(np.max(np.abs(diff))),
        "rmse": float(np.sqrt(np.mean(diff**2))),
    }


# ---------------------------------------------------------------------------
# Deterministic graphical methods: HYSEP fixed / sliding / local, UKIH.
# No free parameters beyond drainage area -- exact agreement is expected.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", ["fixed", "slide", "local", "ukih"])
def test_graphical_methods_match_xie2020(validation_gages, name):
    for site_id, Q, area, _ in validation_gages:
        b_lh_bfx = separation.lh(Q)
        b_lh_xie = xie.LH(Q)
        # Both packages implement LH identically; this just confirms it
        # before using it as the edge-fill input to local()/ukih().
        np.testing.assert_allclose(b_lh_bfx, b_lh_xie, rtol=RTOL, atol=ATOL)

        if name == "fixed":
            b_bfx = separation.fixed(Q, area)
            b_xie = xie.Fixed(Q, area)
        elif name == "slide":
            b_bfx = separation.slide(Q, area)
            b_xie = xie.Slide(Q, area)
        elif name == "local":
            b_bfx = separation.local(Q, b_lh_bfx, area)
            b_xie = xie.Local(Q, b_lh_xie, area)
        elif name == "ukih":
            b_bfx = separation.ukih(Q, b_lh_bfx)
            b_xie = xie.UKIH(Q, b_lh_xie)

        m = _metrics(b_bfx, b_xie)
        assert m["max_abs_dev"] < ATOL, (
            f"{name} @ {site_id}: max abs deviation {m['max_abs_dev']:.3e} cfs "
            f"exceeds floating-point tolerance"
        )


# ---------------------------------------------------------------------------
# Recursive digital filters (gamma=0 and gamma=1 families). Parameters are
# fixed (not calibrated) so the test isolates the recursion itself; a
# separate test below checks that calibration/estimation routines
# (recession_coefficient, param_calibrate) also agree between packages.
# ---------------------------------------------------------------------------

_FIXED_PARAMS = {
    "chapman_maxwell": {"a": 0.93},
    "eckhardt_08": {"a": 0.93, "BFImax": 0.80},
    "eckhardt_05": {"a": 0.93, "BFImax": 0.50},
    "boughton": {"a": 0.93, "C": 0.10},
    "ewma": {"e": 0.06},
    "furey": {"a": 0.93, "A": 0.10},
    "chapman": {"a": 0.93},
    "willems": {"a": 0.93, "w": 0.10},
}


@pytest.mark.parametrize("name", list(_FIXED_PARAMS))
def test_digital_filters_match_xie2020(validation_gages, name):
    """Reservoir-filter families (gamma=0 and gamma=1) require b[0] to be
    seeded consistently. Xie et al.'s package hardcodes b[0] = b_LH[0], the
    two-pass Lyne-Hollick result -- *not* Q[0] (the backward pass moves
    b_LH[0] away from Q[0]). baseflowx exposes this as
    ``initial_method='LH'`` (baseflowx defaults to the simpler 'Q0' instead,
    which is a legitimate but different design choice, not a bug -- see
    ``_init_baseflow`` in separation.py). Matching Xie's convention here
    isolates the recursion itself as the thing under test.
    """
    for site_id, Q, area, _ in validation_gages:
        b_lh = xie.LH(Q)
        p = _FIXED_PARAMS[name]

        if name == "chapman_maxwell":
            b_bfx = separation.chapman_maxwell(Q, p["a"], initial_method="LH")
            b_xie = xie.CM(Q, b_lh, p["a"])
        elif name.startswith("eckhardt"):
            b_bfx = separation.eckhardt(Q, p["a"], p["BFImax"], initial_method="LH")
            b_xie = xie.Eckhardt(Q, b_lh, p["a"], p["BFImax"])
        elif name == "boughton":
            b_bfx = separation.boughton(Q, p["a"], p["C"], initial_method="LH")
            b_xie = xie.Boughton(Q, b_lh, p["a"], p["C"])
        elif name == "ewma":
            b_bfx = separation.ewma(Q, p["e"], initial_method="LH")
            b_xie = xie.EWMA(Q, b_lh, a=0.93, e=p["e"])
        elif name == "furey":
            b_bfx = separation.furey(Q, p["a"], p["A"], initial_method="LH")
            b_xie = xie.Furey(Q, b_lh, p["a"], p["A"])
        elif name == "chapman":
            b_bfx = separation.chapman(Q, p["a"], initial_method="LH")
            b_xie = xie.Chapman(Q, b_lh, p["a"])
        elif name == "willems":
            b_bfx = separation.willems(Q, p["a"], p["w"], initial_method="LH")
            b_xie = xie.Willems(Q, b_lh, p["a"], p["w"])

        m = _metrics(b_bfx, b_xie)
        assert m["max_abs_dev"] < ATOL, (
            f"{name} @ {site_id}: max abs deviation {m['max_abs_dev']:.3e} cfs "
            f"exceeds floating-point tolerance"
        )


def test_what_alias_matches_eckhardt(fish_river_Q):
    """what() is documented as an alias for eckhardt(); confirm the alias holds."""
    Q = fish_river_Q
    a, BFImax = 0.93, 0.5
    np.testing.assert_array_equal(
        separation.what(Q, BFImax, a), separation.eckhardt(Q, a, BFImax)
    )


# ---------------------------------------------------------------------------
# Parameter estimation: recession_coefficient() and param_calibrate() are
# also shared with the reference package. If these disagree, every filter
# above would receive different parameters in practice even though the
# recursions themselves match.
# ---------------------------------------------------------------------------

def test_recession_coefficient_matches_xie2020(validation_gages):
    for site_id, Q, area, _ in validation_gages:
        strict = separation.strict_baseflow(Q)
        a_bfx = estimate.recession_coefficient(Q, strict)
        a_xie = xie.recession_coefficient(Q, strict)
        assert a_bfx == pytest.approx(a_xie, rel=RTOL), (
            f"recession_coefficient @ {site_id}: {a_bfx} vs {a_xie}"
        )


def test_param_calibrate_matches_xie2020(fish_river_Q):
    """Regression test for the param_calibrate() signature fix.

    baseflowx.estimate.param_calibrate previously called
    ``method(Q, b_LH, a, p, return_exceed=True)`` -- the calling convention
    of the *original* Xie et al. (2020) filter functions, which take b_LH
    as an explicit initial-condition argument. baseflowx's refactored
    filters (separation.boughton, etc.) instead self-seed b[0] via
    ``initial_method='Q0'`` and no longer accept b_LH positionally, so the
    old call raised ValueError for every baseflowx filter -- param_calibrate
    was unusable. Fixed to call ``method(Q, a, p, return_exceed=True)``.
    This test both confirms the fix (no exception) and cross-checks the
    result against the reference package's own calibration loop, since
    the two now use different calling conventions for `method` but should
    still search the same loss surface and find the same optimum.
    """
    Q = fish_river_Q
    strict = separation.strict_baseflow(Q)
    a = estimate.recession_coefficient(Q, strict)
    b_lh = separation.lh(Q)
    C_range = np.arange(0.01, 1.00, 0.01)

    C_bfx = estimate.param_calibrate(C_range, separation.boughton, Q, a)
    C_xie = xie.param_calibrate(C_range, xie.Boughton, Q, b_lh, a)
    assert C_bfx == pytest.approx(C_xie, abs=1e-9)
