"""Run the conforming model over the whole stance, solving the sink each row.

For every STEP-th row of foot_motion.csv: find the sink (mm) at which the summed
slice forces equal the measured F (brentq), then keep each slice's force and
peak column strain. Writes full_stance.csv (long format) and a strain map.

Run:  python full_stance.py          (every 2nd row)
      python full_stance.py 1        (every row)
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from fill_midsole import H, L, T, THETA, read_volume
from foot_conform import ES, conform, foot_bodies, pose, slice_beam
from lattice_law import x2_law

HERE = OUTPUTS
OUT = HERE / "full_stance.csv"
F_MIN = 30.0         # N, skip rows with almost no load
SINK_MAX = 20.0      # mm, search limit


def main():
    step = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    mo = pd.read_csv(HERE / "foot_motion.csv")
    feet = foot_bodies()
    outline, z_top = read_volume()
    law = x2_law(THETA, H, L, T, ES)
    rows, t0 = [], time.time()
    for k in range(0, len(mo), step):
        r = mo.loc[k]
        if r.F < F_MIN:
            continue
        sink, tops = sink_solve(feet, outline, z_top, law, r)
        tot = sum(t["force"] for t in tops)
        print(f"row {k:3d} ({r.pct:5.1f} %): F {r.F:6.0f} N, sink {sink:5.2f} mm, "
              f"model {tot:6.0f} N, peak strain {max(t['eps'].max() for t in tops) * 100:4.1f} %")
        for t in tops:
            rows.append({"row": k, "pct": r.pct, "F": r.F, "sink": sink, "x": t["x"],
                         "force": t["force"], "eps_max": t["eps"].max()})
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, float_format="%.4f")
    print(f"done in {time.time() - t0:.0f} s -> {OUT.name}")
    summary(df)


def total_force(feet, outline, z_top, law, r, sink):
    tops = conform(pose(feet, r, sink), outline, z_top)
    for t in tops:
        slice_beam(t, law(t["x"]) if callable(law) else law)   # law per zone or one law
    return sum(t["force"] for t in tops), tops


def sink_solve(feet, outline, z_top, law, r):
    """Sink where the slice forces add up to the measured F."""
    res = lambda s: total_force(feet, outline, z_top, law, r, s)[0] - r.F
    lo = -5.0                                  # foot raised: surely no contact
    if res(SINK_MAX) < 0:
        print(f"  warning: even {SINK_MAX} mm sink does not carry {r.F:.0f} N")
        s = SINK_MAX
    else:
        s = brentq(res, lo, SINK_MAX, xtol=0.01)
    return s, total_force(feet, outline, z_top, law, r, s)[1]


def summary(df):
    """Peak per slice over the stance, and a strain map (slice x vs stance %)."""
    pk = df.groupby("x").agg(force_max=("force", "max"), eps_max=("eps_max", "max"))
    print("\npeak over stance, per slice (every 3rd):")
    print(pk.iloc[::3].round(3).to_string())
    m = df.pivot(index="x", columns="pct", values="eps_max") * 100
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.pcolormesh(m.columns, m.index, m.values, cmap="magma", shading="nearest")
    fig.colorbar(im, label="peak column strain (%)")
    ax.set_xlabel("stance (%)"); ax.set_ylabel("slice x (mm)  toe -> heel")
    ax.set_title("Lattice strain over the stance")
    fig.tight_layout()
    fig.savefig(HERE / "full_stance.png", dpi=120)


if __name__ == "__main__":
    main()
