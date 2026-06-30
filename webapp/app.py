"""Baseflow Explorer — Flask + Leaflet web app."""
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from dataretrieval import nwis
from flask import Flask, jsonify, render_template, request

import baseflowx as bf

app = Flask(__name__)

DATA_DIR = Path(__file__).parent / "data"
SITES_PATH = DATA_DIR / "nwis_dv_sites.parquet"
SQMI_TO_KM2 = 2.58999

LABELED_SITES = {"01013500"}

METHODS = {
    "lh":              {"label": "Lyne-Hollick"},
    "eckhardt":        {"label": "Eckhardt"},
    "chapman":         {"label": "Chapman"},
    "chapman_maxwell": {"label": "Chapman-Maxwell"},
    "bflow":           {"label": "BFlow"},
    "fixed":           {"label": "Fixed interval"},
    "slide":           {"label": "Sliding interval"},
    "local":           {"label": "Local minimum"},
    "ukih":            {"label": "UKIH"},
    "part":            {"label": "PART"},
}

_sites_df = pd.read_parquet(SITES_PATH)
_dv_cache: dict[tuple, pd.DataFrame] = {}


def _fetch_dv(site_id: str, start: str, end: str) -> pd.DataFrame:
    key = (site_id, start, end)
    if key in _dv_cache:
        return _dv_cache[key]
    df, _ = nwis.get_dv(sites=site_id, start=start, end=end, parameterCd="00060")
    _dv_cache[key] = df
    return df


@app.route("/")
def index():
    return render_template(
        "index.html",
        version=bf.__version__,
        methods=METHODS,
        labeled_sites=list(LABELED_SITES),
    )


@app.route("/sites.json")
def sites_json():
    rows = _sites_df[["site_no", "station_nm", "dec_lat_va", "dec_long_va", "drain_area_sqmi"]]
    rows = rows.replace({np.nan: None})
    return jsonify(rows.to_dict("records"))


@app.route("/analyze")
def analyze():
    try:
        site_id = request.args["site_id"]
        start = request.args.get("start", "2020-10-01")
        end = request.args.get("end", date.today().isoformat())
        methods_arg = request.args.get("methods", "lh,eckhardt,chapman_maxwell,ukih,part")
        methods_req = [m for m in methods_arg.split(",") if m]
        beta = float(request.args.get("beta", 0.925))
        a_coef = float(request.args.get("a", 0.925))
        bfi_max = float(request.args.get("bfi_max", 0.80))
        area_km2 = float(request.args.get("area", 500.0))
    except (KeyError, ValueError) as e:
        return jsonify({"error": f"Bad request: {e}"}), 400

    try:
        df = _fetch_dv(site_id, start, end)
    except Exception as e:
        return jsonify({"error": f"NWIS fetch failed: {e}"}), 502

    if df is None or df.empty:
        return jsonify({"error": "No daily-values data for this gage and date range."}), 404

    q_col = next(
        (c for c in df.columns if c.startswith("00060") and "Mean" in c and "cd" not in c),
        None,
    )
    if q_col is None:
        return jsonify({"error": "No discharge column in NWIS response."}), 500

    q_raw = df[q_col].astype(float).to_numpy()
    dates_idx = df.index
    mask = np.isfinite(q_raw) & (q_raw > 0)
    Q = q_raw[mask]
    dates_used = pd.to_datetime(dates_idx[mask])

    if len(Q) < 30:
        return jsonify({"error": f"Only {len(Q)} valid data points — too few to analyze."}), 400

    b_lh = bf.lh(Q, beta=beta)
    runners = {
        "lh":              lambda: bf.lh(Q, beta=beta),
        "eckhardt":        lambda: bf.eckhardt(Q, a=a_coef, BFImax=bfi_max),
        "chapman":         lambda: bf.chapman(Q, a=a_coef),
        "chapman_maxwell": lambda: bf.chapman_maxwell(Q, a=a_coef),
        "bflow":           lambda: bf.bflow(Q)["baseflow"],
        "fixed":           lambda: bf.fixed(Q, area=area_km2),
        "slide":           lambda: bf.slide(Q, area=area_km2),
        "local":           lambda: bf.local(Q, b_lh, area=area_km2),
        "ukih":            lambda: bf.ukih(Q, b_lh),
        "part":            lambda: bf.part(Q, area=area_km2),
    }

    baseflow: dict[str, list[float]] = {}
    bfi: dict[str, float | None] = {}
    failed: dict[str, str] = {}
    total_q = float(np.sum(Q))
    for m in methods_req:
        if m not in runners:
            continue
        try:
            b = np.asarray(runners[m]())[: len(Q)]
            baseflow[m] = b.tolist()
            bfi[m] = float(np.sum(b) / total_q)
        except Exception as e:
            failed[m] = str(e)

    site_row = _sites_df[_sites_df.site_no == site_id]
    station_nm = site_row.iloc[0].station_nm if not site_row.empty else ""

    return jsonify({
        "site_id": site_id,
        "station_nm": station_nm,
        "dates": [d.isoformat() for d in dates_used.to_pydatetime()],
        "Q": Q.tolist(),
        "baseflow": baseflow,
        "bfi": bfi,
        "failed": failed,
        "n_days": len(Q),
        "mean_q": float(np.mean(Q)),
    })


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080, debug=True)
