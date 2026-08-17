"""Correctness checks for methods with no external reference implementation.

Twelve of baseflowx's 17 methods are cross-checked against the independently
packaged Xie et al. (2020) `baseflow` implementation in
`test_reference_implementations.py`. The remaining five -- PART, IHACRES,
bn77, strict_baseflow, and the WHAT/Eckhardt and Chapman-Maxwell/Eckhardt
aliases -- have no such external Python reference available. This module
validates them a different way: against documented mathematical identities
(two filters that are provably the same equation in different parameter
dress), against hand-derived expected output on small synthetic series where
the correct answer can be worked out by hand, and against physical
invariants that must hold regardless of implementation (baseflow can never
exceed streamflow; a pure monotonic recession must be entirely classified
as baseflow).
"""

import numpy as np
import pytest

from baseflowx import estimate, separation


# ---------------------------------------------------------------------------
# Documented algebraic identities between filters
# ---------------------------------------------------------------------------

def test_chapman_maxwell_equals_eckhardt_at_BFImax_half(fish_river_Q):
    """docs/methods/digital-filters.md and the chapman_maxwell() docstring
    both state this equivalence; confirm it numerically."""
    Q = fish_river_Q
    a = 0.93
    b_cm = separation.chapman_maxwell(Q, a)
    b_eck = separation.eckhardt(Q, a, BFImax=0.5)
    np.testing.assert_allclose(b_cm, b_eck, rtol=1e-12, atol=1e-12)


def test_ihacres_reduces_to_boughton_at_alpha_s_zero(fish_river_Q):
    """ihacres() docstring: 'When alpha_s = 0, reduces exactly to boughton().'"""
    Q = fish_river_Q
    a, C = 0.93, 0.15
    b_boughton = separation.boughton(Q, a, C)
    b_ihacres = separation.ihacres(Q, a, C, alpha_s=0.0)
    np.testing.assert_allclose(b_boughton, b_ihacres, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Physical constraint: baseflow can never exceed streamflow
# ---------------------------------------------------------------------------

_DIGITAL_FILTERS = {
    "chapman_maxwell": lambda Q: separation.chapman_maxwell(Q, 0.93),
    "eckhardt": lambda Q: separation.eckhardt(Q, 0.93, 0.8),
    "boughton": lambda Q: separation.boughton(Q, 0.93, 0.2),
    "ewma": lambda Q: separation.ewma(Q, 0.05),
    "furey": lambda Q: separation.furey(Q, 0.93, 0.1),
    "chapman": lambda Q: separation.chapman(Q, 0.93),
    "willems": lambda Q: separation.willems(Q, 0.93, 0.2),
    "ihacres": lambda Q: separation.ihacres(Q, 0.93, 0.2, -0.3),
    "lh": lambda Q: separation.lh(Q),
    "lh_multi_3pass": lambda Q: separation.lh_multi(Q, num_pass=3),
}


@pytest.mark.parametrize("name", list(_DIGITAL_FILTERS))
def test_filter_never_exceeds_streamflow(validation_gages, name):
    for site_id, Q, area, _ in validation_gages:
        b = _DIGITAL_FILTERS[name](Q)
        assert np.all(b <= Q + 1e-9), (
            f"{name} @ {site_id}: baseflow exceeds streamflow by "
            f"{np.max(b - Q):.3e} at the worst timestep"
        )
        assert np.all(b >= 0), f"{name} @ {site_id}: negative baseflow produced"


@pytest.mark.parametrize("name", ["fixed", "slide", "local", "ukih", "part"])
def test_graphical_method_never_exceeds_streamflow(validation_gages, name):
    for site_id, Q, area, _ in validation_gages:
        b_lh = separation.lh(Q)
        if name == "fixed":
            b = separation.fixed(Q, area)
        elif name == "slide":
            b = separation.slide(Q, area)
        elif name == "local":
            b = separation.local(Q, b_lh, area)
        elif name == "ukih":
            b = separation.ukih(Q, b_lh)
        elif name == "part":
            b = separation.part(Q, area)
        assert np.all(b <= Q + 1e-6), (
            f"{name} @ {site_id}: baseflow exceeds streamflow by "
            f"{np.max(b - Q):.3e} at the worst timestep"
        )


# ---------------------------------------------------------------------------
# PART: hand-derived small example (Rutledge, 1998, as implemented per
# docs/methods/graphical-methods.md)
# ---------------------------------------------------------------------------

def test_part_pure_monotonic_recession_is_entirely_baseflow():
    """If a record declines monotonically from start to end, every day (after
    the first N) satisfies PART's antecedent-recession qualifying condition,
    and the 0.1 log-cycle safeguard is never triggered for a gentle enough
    decay. Every qualifying day gets b_t = Q_t exactly, and log-linear
    interpolation between two equal-valued anchors changes nothing --
    so baseflow should equal streamflow (to numerical precision) everywhere
    except a short warm-up at the very start, which PART fills by constant
    extrapolation from the first anchor.
    """
    n = 200
    Q = 1000.0 * 0.995 ** np.arange(n)  # ~0.5%/day decline: well under 0.1 log-cycle/day
    area = 259.0  # km^2 -> N = (0.3861*259)^0.2 ~ 2.5 days, small warm-up

    b = separation.part(Q, area)

    # After the initial anchor, PART should track Q almost exactly.
    idx_first_anchor = np.argmax(b == Q)
    np.testing.assert_allclose(b[idx_first_anchor:], Q[idx_first_anchor:], rtol=1e-6)
    # Never exceeds Q, matching the physical constraint.
    assert np.all(b <= Q + 1e-6)


def test_part_log_cycle_threshold_disqualifies_steep_decline():
    """A single-day drop steeper than 0.1 log cycles (~21%) must disqualify
    that day as a PART anchor, per Rutledge (1998) / Barnes (1939) as
    described in docs/methods/graphical-methods.md. Construct a record with
    one such drop and confirm the qualifying-day logic (not just the
    interpolation) reacts to it: baseflow at that day must come from
    interpolation, not from b_t = Q_t.
    """
    n = 60
    Q = 1000.0 * 0.995 ** np.arange(n)
    steep_day = 30
    Q[steep_day] = Q[steep_day - 1] * 0.5  # 50% drop >> 21% (0.1 log-cycle) threshold
    area = 259.0

    b = separation.part(Q, area)
    # The day *before* the steep drop is disqualified by the log-cycle check
    # (its decline into day `steep_day` exceeds the threshold), so its
    # baseflow comes from interpolation rather than equaling Q exactly.
    assert not np.isclose(b[steep_day - 1], Q[steep_day - 1], rtol=1e-9)


# ---------------------------------------------------------------------------
# strict_baseflow / bn77: recession-only invariants
# ---------------------------------------------------------------------------

def test_strict_baseflow_points_are_declining(validation_gages):
    """Every point flagged as strict baseflow must lie on a declining or
    flat stretch of the hydrograph -- strict_baseflow's own first rule
    excludes any point with a non-negative centered derivative, so no
    flagged point should coincide with a local rise wider than the
    smoothing window."""
    for site_id, Q, area, _ in validation_gages:
        strict = separation.strict_baseflow(Q)
        dQ = np.concatenate([[0], (Q[2:] - Q[:-2]) / 2, [0]])
        flagged = np.where(strict)[0]
        assert np.all(dQ[flagged] < 0), (
            f"strict_baseflow @ {site_id}: flagged a non-declining point"
        )


def test_bn77_drought_points_are_subset_of_record(validation_gages):
    for site_id, Q, area, _ in validation_gages:
        idx = separation.bn77(
            Q, L_min=10, snow_freeze_period=(335, 90),
            observational_precision=0.01, quantile=0.9)
        assert idx.min() >= 0 and idx.max() < len(Q)
        assert len(idx) == len(set(idx.tolist())), "bn77 returned duplicate indices"


# ---------------------------------------------------------------------------
# bflow(): each successive Lyne-Hollick pass should remove more quickflow
# ---------------------------------------------------------------------------

def test_bflow_bfi_decreases_across_passes(fish_river_Q):
    result = estimate.bflow(fish_river_Q)
    assert result["BFI_pass1"] >= result["BFI_pass2"] >= result["BFI"], (
        "docs/methods/recession-analysis.md documents BFI_pass1 >= BFI_pass2 "
        f">= BFI (three-pass); got {result['BFI_pass1']}, "
        f"{result['BFI_pass2']}, {result['BFI']}"
    )
