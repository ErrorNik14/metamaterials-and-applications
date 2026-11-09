"""Pick lattice parameters per zone (forefoot / midfoot / heel) so the peak
column strain over the whole stance sits near a target.

Loop:  run the stance with the current zone lattices (sink solved each row)
       -> peak strain per zone -> for each zone pick the candidate whose
       E2 is closest to E2_now * (peak / target), most auxetic among near ties
       -> repeat until nothing changes.
Candidates: grid over theta, l, h/l, t with M&E E2 / nu21, same limits as before.
Writes zone_pick.csv (the chosen lattice per zone).

Run:  python zone_search.py
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

from pathlib import Path

import numpy as np
import pandas as pd

from foot_conform import ES, foot_bodies
from fill_midsole import read_volume
from full_stance import F_MIN, sink_solve
from lattice_law import x2_law

HERE = OUTPUTS
ZONES = {"fore": (0.0, 63.1), "mid": (63.1, 186.6), "heel": (186.6, 300.0)}   # hinge x
TARGET = 0.175          # peak column strain aimed for
E2_FLOOR = 0.3          # MPa, softest lattice allowed (lightly loaded heel)
STEP = 3                # every 3rd stance row while searching
START = (-30.0, 4.0, 8.0, 1.4)   # theta, l, h, t: the current single lattice


def main():
    cand = candidates()
    print(f"{len(cand)} candidate lattices, E2 {cand.E2.min():.2f}-{cand.E2.max():.2f} MPa")
    pick = {z: cand.loc[closest(cand, START)] for z in ZONES}
    mo = pd.read_csv(HERE / "foot_motion.csv")
    feet = foot_bodies()
    outline, z_top = read_volume()
    for it in range(5):
        peak = stance_peaks(pick, mo, feet, outline, z_top)
        print(f"iter {it}: " + "  ".join(
            f"{z} E2 {pick[z].E2:.2f} -> peak {peak[z] * 100:4.1f} %" for z in ZONES))
        new = {}
        for z in ZONES:
            want = max(pick[z].E2 * peak[z] / TARGET, E2_FLOOR)
            near = cand[(cand.E2 > 0.9 * want) & (cand.E2 < 1.1 * want)]
            new[z] = near.loc[near.nu21.idxmin()] if len(near) else \
                cand.loc[(cand.E2 - want).abs().idxmin()]
        if all(new[z].name == pick[z].name for z in ZONES):
            break
        pick = new
    out = pd.DataFrame({z: p for z, p in pick.items()}).T
    out["peak_strain"] = [peak[z] for z in out.index]
    out.to_csv(HERE / "zone_pick.csv", float_format="%.4f")
    print("\nzone pick:\n", out.round(3).to_string())


def me(th, h, l, t):
    """Masters & Evans E2 and nu21 (rigid joints), as masters_evans.m."""
    a, g = h / l, t / l
    s, c = np.sin(np.radians(th)), np.cos(np.radians(th))
    Cp, Ca, Cv = 1.0, g * g, 2 * a * g * g
    D2 = c * c * Cp + s * s * Ca + Cv
    return ES * g**3 * (a + s) / (c * D2), s * (a + s) * (Cp - Ca) / D2


def candidates():
    rows = []
    for th in (-10, -15, -20, -25, -30):
        for l in (3.5, 4.0, 5.0, 6.0, 7.0):
            for ha in (1.0, 1.3, 1.6, 2.0):
                for t in np.round(np.arange(1.0, 2.61, 0.2), 2):
                    h = ha * l
                    s = np.sin(np.radians(th))
                    if t / l > 0.40 or h + 2 * l * s - t < 2 or 2 * (h + l * s) > 16:
                        continue
                    E2, nu = me(th, h, l, t)
                    rows.append((th, l, h, t, E2, nu))
    return pd.DataFrame(rows, columns=["theta", "l", "h", "t", "E2", "nu21"])


def closest(cand, p):
    th, l, h, t = p
    return ((cand.theta - th)**2 + (cand.l - l)**2 + (cand.h - h)**2 + (cand.t - t)**2).idxmin()


def zone_of(x):
    return next(z for z, (a, b) in ZONES.items() if a <= x < b)


def stance_peaks(pick, mo, feet, outline, z_top):
    laws = {z: x2_law(p.theta, p.h, p.l, p.t, ES) for z, p in pick.items()}
    law = lambda x: laws[zone_of(x)]
    peak = {z: 0.0 for z in ZONES}
    for k in range(0, len(mo), STEP):
        r = mo.loc[k]
        if r.F < F_MIN:
            continue
        _, tops = sink_solve(feet, outline, z_top, law, r)
        for t in tops:
            z = zone_of(t["x"])
            peak[z] = max(peak[z], t["eps"].max())
    return peak


if __name__ == "__main__":
    main()
