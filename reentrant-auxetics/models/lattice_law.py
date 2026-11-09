"""Ohtaki large-deflection curve, ported from ohtaki.m (same equations, same stops).

Run:  python lattice_law.py      -> x2 compression to 30 %, prints and caches it
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

CACHE = OUTPUTS / "ohtaki_cache"


def oh_curve(theta_deg, h, l, t, Es, dirn, eps_end, f_step=0.05):
    """Strain, stress (MPa) and secant nu along x1 (dirn 1) or x2 (dirn 2).
    eps_end > 0 tension, < 0 compression. See ohtaki.m for the derivation."""
    a, g = h / l, t / l
    th = np.radians(theta_deg)
    s, c = np.sin(th), np.cos(th)
    sgn = np.sign(eps_end)

    def run(F, k0):
        rhs = lambda _, y: [y[1], F[0] * np.sin(y[0]) - F[1] * np.cos(y[0]),
                            np.cos(y[0]), np.sin(y[0])]
        return solve_ivp(rhs, [0, 1], [th, k0, 0, 0], rtol=1e-9, atol=1e-12).y[:, -1]

    eps, sig, nu = [0.0], [0.0], [np.nan]
    k0, f = 0.0, 0.0
    step = f_step
    while True:
        F = [0.0, 0.0]
        F[dirn - 1] = f + sgn * step
        k_try = fsolve(lambda k: run(F, k[0])[1], [k0], xtol=1e-10)[0]
        _, _, dx, dy = run(F, k_try)
        e1 = (dx - c) / c
        e2 = (2 * dy - 2 * s) / (2 * a + 2 * s)
        fn = F[dirn - 1]
        if dirn == 1:
            sg, e_load, e_lat = fn * Es * g**3 / (3 * (a + s)), e1, e2
        else:
            sg, e_load, e_lat = fn * Es * g**3 / (3 * c), e2, e1
        # too big a strain jump: halve the force step and retry (the .m file
        # just stopped; x2 is soft, so a fixed step overshoots at the start)
        if abs(e_load - eps[-1]) > 0.01 or np.sign(e_load) != sgn:
            if step > f_step / 256:
                step /= 2
                continue
            print("  ohtaki stopped: solver jumped / strain limit")
            break
        f, k0 = fn, k_try
        eps.append(e_load); sig.append(sg); nu.append(-e_lat / e_load)
        if abs(e_load) >= abs(eps_end):
            break
    return np.array(eps), np.array(sig), np.array(nu)


def x2_law(theta_deg, h, l, t, Es, eps_max=0.3):
    """Compression stress magnitude vs compression strain along x2, cached."""
    key = f"{theta_deg}_{h:.3f}_{l:.3f}_{t:.3f}_{Es:.3f}_{eps_max}"
    CACHE.mkdir(exist_ok=True)
    fn = CACHE / f"x2_{key}.npy"
    if fn.exists():
        e, s = np.load(fn)
    else:
        e, s, _ = oh_curve(theta_deg, h, l, t, Es, 2, -eps_max)
        e, s = -e, -s
        if len(e) > 2:          # don't cache a failed run
            np.save(fn, np.vstack([e, s]))
    return e, s


if __name__ == "__main__":
    e, s = x2_law(-20, 1.3 * 7, 7, 1.1, 19.3)
    print(f"x2 compression: {len(e)} points, up to eps {e[-1]:.3f}, sigma {s[-1]:.4f} MPa")
    print(f"  initial E2 ~ {s[1] / e[1]:.4f} MPa")
