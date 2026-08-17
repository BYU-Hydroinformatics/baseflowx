#!/usr/bin/env python3
"""Regenerate docs/validation.md from a live run against the reference package.

Run from the repo root:

    pip install -e .[validation]
    python scripts/validation_report.py

Requires network access on first run of the six validation gages'
discharge fetch is already cached under tests/data/ and baseflowx/data/,
so this script itself does not touch the network.
"""

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import baseflow as xie  # noqa: E402
from baseflowx import estimate, separation  # noqa: E402
from conftest import VALIDATION_GAGES, _load_discharge  # noqa: E402


def metrics(a, b):
    diff = np.asarray(a) - np.asarray(b)
    return float(np.max(np.abs(diff))), float(np.sqrt(np.mean(diff**2)))


def load_gages():
    out = []
    for site_id, path, area, ecoregion in VALIDATION_GAGES:
        out.append((site_id, _load_discharge(path), area, ecoregion))
    return out


GRAPHICAL = ["fixed", "slide", "local", "ukih"]
FIXED_PARAMS = {
    "chapman_maxwell": {"a": 0.93},
    "eckhardt (BFImax=0.80)": {"a": 0.93, "BFImax": 0.80},
    "eckhardt (BFImax=0.50)": {"a": 0.93, "BFImax": 0.50},
    "boughton": {"a": 0.93, "C": 0.10},
    "ewma": {"e": 0.06},
    "furey": {"a": 0.93, "A": 0.10},
    "chapman": {"a": 0.93},
    "willems": {"a": 0.93, "w": 0.10},
}


def run_graphical(name, Q, area, b_lh_bfx, b_lh_xie):
    if name == "fixed":
        return separation.fixed(Q, area), xie.Fixed(Q, area)
    if name == "slide":
        return separation.slide(Q, area), xie.Slide(Q, area)
    if name == "local":
        return separation.local(Q, b_lh_bfx, area), xie.Local(Q, b_lh_xie, area)
    if name == "ukih":
        return separation.ukih(Q, b_lh_bfx), xie.UKIH(Q, b_lh_xie)


def run_filter(name, Q, p, b_lh):
    a = p.get("a")
    if name == "chapman_maxwell":
        return (separation.chapman_maxwell(Q, a, initial_method="LH"),
                xie.CM(Q, b_lh, a))
    if name.startswith("eckhardt"):
        return (separation.eckhardt(Q, a, p["BFImax"], initial_method="LH"),
                xie.Eckhardt(Q, b_lh, a, p["BFImax"]))
    if name == "boughton":
        return (separation.boughton(Q, a, p["C"], initial_method="LH"),
                xie.Boughton(Q, b_lh, a, p["C"]))
    if name == "ewma":
        return (separation.ewma(Q, p["e"], initial_method="LH"),
                xie.EWMA(Q, b_lh, a=0.93, e=p["e"]))
    if name == "furey":
        return (separation.furey(Q, a, p["A"], initial_method="LH"),
                xie.Furey(Q, b_lh, a, p["A"]))
    if name == "chapman":
        return (separation.chapman(Q, a, initial_method="LH"),
                xie.Chapman(Q, b_lh, a))
    if name == "willems":
        return (separation.willems(Q, a, p["w"], initial_method="LH"),
                xie.Willems(Q, b_lh, a, p["w"]))


def main():
    gages = load_gages()
    rows = []

    for name in GRAPHICAL:
        worst_dev, worst_rmse = 0.0, 0.0
        for site_id, Q, area, _ in gages:
            b_lh_bfx = separation.lh(Q)
            b_lh_xie = xie.LH(Q)
            b_bfx, b_xie = run_graphical(name, Q, area, b_lh_bfx, b_lh_xie)
            dev, rmse = metrics(b_bfx, b_xie)
            worst_dev, worst_rmse = max(worst_dev, dev), max(worst_rmse, rmse)
        rows.append((name, worst_dev, worst_rmse))

    for name, p in FIXED_PARAMS.items():
        worst_dev, worst_rmse = 0.0, 0.0
        for site_id, Q, area, _ in gages:
            b_lh = xie.LH(Q)
            b_bfx, b_xie = run_filter(name, Q, p, b_lh)
            dev, rmse = metrics(b_bfx, b_xie)
            worst_dev, worst_rmse = max(worst_dev, dev), max(worst_rmse, rmse)
        rows.append((name, worst_dev, worst_rmse))

    print(f"{'method':<26}{'max |dev| (cfs)':>18}{'max RMSE (cfs)':>18}")
    print("-" * 62)
    for name, dev, rmse in rows:
        print(f"{name:<26}{dev:>18.2e}{rmse:>18.2e}")

    n_gages = len(gages)
    print(f"\n{n_gages} gages: " + ", ".join(f"{s}" for s, *_ in gages))


if __name__ == "__main__":
    main()
