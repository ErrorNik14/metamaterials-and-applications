"""Make the midsole top follow the posed pseudo-foot (start of the slice-loading step).

For one row of foot_motion.csv:
  1. pose the 3D pseudo-foot (same rules as foot_anim.py),
  2. for every honeycomb slice, sample the foot underside z_foot(y) by casting
     rays straight down, and compare it with the midsole top z_top(y),
  3. where the foot dips below the top, the top is pushed down to it:
        new top = min(z_top, z_foot),  depth = z_top - new top  (>= 0)
     Where the foot is above (e.g. heel lifted) the midsole stays as it is.
  4. each slice top is then solved as a beam on the lattice (Winkler foundation,
     Ohtaki x2 compression) pushed by the rigid foot -> deflection, slice force.
Shows the result in pyvista: posed foot, midsole, top colored by deflection.

The motion data does not say how far the foot sinks into the sole; that comes
later from the honeycomb stiffness (sink until the slice forces add up to F).
For now SINK (mm) is set by hand.

Run:  python foot_conform.py 30        (row)
      python foot_conform.py 30 3      (row, sink in mm)
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyvista as pv
import trimesh

from fill_midsole import B, GAP, H, L, PLATE_BOT, PLATE_TOP, T, THETA, read_volume
from lattice_law import x2_law

HERE = OUTPUTS
ANKLE = np.array([231.2, 82.9])
MTP = np.array([63.1, 52.0])
NY = 40                    # sample points across each slice
ES = 19.3                  # solid modulus (MPa), same as the model tests
K_PEN = 1e4                # contact penalty spring (N/mm), stiff = foot is rigid


def main():
    row = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    sink = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    mo = pd.read_csv(HERE / "foot_motion.csv")
    feet = foot_bodies()
    outline, z_top = read_volume()
    posed = pose(feet, mo.loc[row], sink)
    tops = conform(posed, outline, z_top)
    d = np.concatenate([t["depth"] for t in tops])
    g = np.concatenate([t["gap"] for t in tops])
    print(f"smallest foot-to-top gap after sinking: {np.nanmin(g):.1f} mm (negative = overlap)")
    print(f"row {row} ({mo.pct[row]:.0f} % stance), sink {sink} mm: {len(tops)} slices, "
          f"max depth {d.max():.1f} mm, in contact at {np.mean(d > 0) * 100:.0f} % of points")
    for t in tops[::3]:
        if t["depth"].max() > 0:
            print(f"  slice x {t['x']:5.1f}: max depth {t['depth'].max():5.1f} mm")
    law = x2_law(THETA, H, L, T, ES)
    for t in tops:
        slice_beam(t, law)
    tot = sum(t["force"] for t in tops)
    print(f"beam solve: total slice force {tot:.0f} N  vs  F = {mo.F[row]:.0f} N in the data")
    for t in tops:
        if t["force"] > 1:
            print(f"  slice x {t['x']:5.1f}: force {t['force']:6.1f} N, "
                  f"max deflection {t['w'].max():4.2f} mm, max strain {t['eps'].max() * 100:4.1f} %")
    show(posed, tops, row, mo)


def foot_bodies():
    """Pseudo-foot as three trimesh bodies in midsole coordinates."""
    t = (INPUTS / "pseudo_foot.step").read_text(encoding="latin-1")
    pts = {int(i): [float(v) for v in x.split(",")]
           for i, x in re.findall(r"#(\d+)=CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", t)}
    V = np.array([pts[int(i)] for i in re.findall(r"=VERTEX_POINT\('[^']*',#(\d+)\)", t)])
    heel = trimesh.creation.icosphere(4, 37.5)
    heel.apply_translation([231.173, 69.357, 82.897])
    return {"toe": trimesh.convex.convex_hull(V[V[:, 0] < 62]),
            "mid": trimesh.convex.convex_hull(V[(V[:, 0] > 62) & (V[:, 0] < 195)]),
            "heel": heel}


def rot_y(c, a_deg):
    """4x4 rotation in the x-z view about a hinge line along y through c = (x, z),
    counter-clockwise as drawn in foot_anim.py (heel side goes up for a > 0)."""
    a = np.radians(a_deg)
    R = np.eye(4)
    R[0, 0], R[0, 2], R[2, 0], R[2, 2] = np.cos(a), -np.sin(a), np.sin(a), np.cos(a)
    T, Ti = np.eye(4), np.eye(4)
    T[[0, 2], 3], Ti[[0, 2], 3] = [c[0], c[1]], [-c[0], -c[1]]
    return T @ R @ Ti


def pose(feet, r, sink=0.0):
    rear = -r.pitch                       # heel+midfoot rigid, about the MTP
    toe = rear - r.mtp
    out = {}
    for name, m in feet.items():
        m = m.copy()
        m.apply_transform(rot_y(MTP, toe if name == "toe" else rear))
        m.apply_translation([0, 0, -sink])
        out[name] = m
    return out


def underside(posed):
    """Function (x, y arrays) -> lowest z of the foot there (inf where no foot).
    Every body is convex, so a vertical line enters it at the highest of the
    downward-facing face planes and leaves at the lowest upward-facing one."""
    planes = []
    for m in posed.values():
        n = m.face_normals
        d = np.einsum("ij,ij->i", n, m.vertices[m.faces[:, 0]])
        planes.append((n, d))

    def z_foot(x, y):
        x, y = np.atleast_1d(x)[:, None], np.atleast_1d(y)[:, None]
        best = np.full(x.shape[0], np.inf)
        for n, d in planes:
            dn, up = n[:, 2] < -1e-9, n[:, 2] > 1e-9
            # plane n.p = d solved for z at (x, y)
            z_lo = ((d[dn] - n[dn, 0] * x - n[dn, 1] * y) / n[dn, 2]).max(axis=1)
            z_hi = ((d[up] - n[up, 0] * x - n[up, 1] * y) / n[up, 2]).min(axis=1)
            hit = z_lo <= z_hi
            best = np.where(hit, np.minimum(best, z_lo), best)
        return best

    return z_foot


def conform(posed, outline, z_top):
    """Per slice: y samples, original top, conformed top, depth."""
    from shapely.geometry import LineString
    z_foot = underside(posed)
    x_min, _, x_max, _ = outline.bounds
    tops = []
    for x in np.arange(x_min + B / 2, x_max - B / 2 + 1e-9, B + GAP):
        cut = outline.intersection(LineString([(x, -1e3), (x, 1e3)]))
        _, ya, _, yb = cut.bounds
        ys = np.linspace(ya, yb, NY)
        zt = z_top(np.full_like(ys, x), ys)
        zf = z_foot(np.full_like(ys, x), ys)
        new = np.minimum(zt, zf)
        tops.append({"x": x, "y": ys, "z_top": zt, "z_new": new, "depth": zt - new,
                     "gap": np.where(np.isinf(zf), np.nan, zf - zt)})
    return tops


def beam_K(n, dl, EI):
    """Euler-Bernoulli beam stiffness, n nodes, dofs (w, theta) per node."""
    k = EI / dl**3 * np.array([[12, 6 * dl, -12, 6 * dl],
                               [6 * dl, 4 * dl**2, -6 * dl, 2 * dl**2],
                               [-12, -6 * dl, 12, -6 * dl],
                               [6 * dl, 2 * dl**2, -6 * dl, 4 * dl**2]])
    K = np.zeros((2 * n, 2 * n))
    for e in range(n - 1):
        K[2 * e:2 * e + 4, 2 * e:2 * e + 4] += k
    return K


def slice_beam(t, law):
    """Top of one slice as a beam across y on a lattice (Winkler) foundation,
    pushed down by the rigid foot.
      beam:        top plate strip (B + GAP wide) + skin wall (B wide), free ends
      foundation:  each node = a lattice column, stress from Ohtaki x2 compression
                   at strain w / column height
      contact:     penalty springs where the foot would sit below the top
    Adds w (deflection, + down), eps (column strain), force (N) to t."""
    e_law, s_law = law
    E0 = s_law[1] / e_law[1]
    y, n = t["y"], len(t["y"])
    dl = y[1] - y[0]
    EI = ES * ((B + GAP) * PLATE_TOP**3 + B * T**3) / 12
    K0 = beam_K(n, dl, EI)
    trib = np.full(n, dl); trib[[0, -1]] = dl / 2
    hcol = np.maximum(t["z_top"] - PLATE_BOT - PLATE_TOP, 1.0)
    g = np.where(np.isnan(t["gap"]), -1e9, -t["gap"])   # needed push-down, > 0 = overlap

    def col_force(w):
        eps = np.abs(w) / hcol
        sig = np.where(eps <= e_law[-1], np.interp(eps, e_law, s_law),
                       s_law[-1] + (eps - e_law[-1]) * (s_law[-1] - s_law[-2]) /
                       (e_law[-1] - e_law[-2]))
        sig = np.where(w < 0, E0 * eps, sig)        # lifted: linear tension
        return np.sign(w) * sig * B * trib

    w = np.zeros(n)
    active = g > 0
    for _ in range(200):
        ks = np.where(np.abs(w) > 1e-9, col_force(w) / np.where(w == 0, 1, w),
                      E0 * B * trib / hcol)           # secant spring per node
        K = K0.copy()
        K[::2, ::2][np.diag_indices(n)] += ks + K_PEN * active
        f = np.zeros(2 * n)
        f[::2] = K_PEN * g * active
        w_new = np.linalg.solve(K + 1e-9 * np.eye(2 * n), f)[::2]
        new_active = g > w_new
        done = np.allclose(w_new, w, atol=1e-5) and np.array_equal(new_active, active)
        w, active = w_new, new_active | (active & (g > w_new - 1e-6))
        if done:
            break
    t["w"], t["eps"] = w, np.maximum(w, 0) / hcol
    t["force"] = float(np.sum(K_PEN * np.maximum(g - w, 0)))
    t["z_new"] = t["z_top"] - w


def show(posed, tops, row, mo):
    p = pv.Plotter()
    colors = {"toe": "#534AB7", "mid": "#0F6E56", "heel": "#993C1D"}
    for name, m in posed.items():
        p.add_mesh(pv.wrap(m), color=colors[name], opacity=0.35)
    sole = pv.read(HERE / "midsole_filled.stl")
    p.add_mesh(sole, color="#BA7517", opacity=0.15)
    pts = np.vstack([np.c_[np.full_like(t["y"], t["x"]), t["y"], t["z_new"]] for t in tops])
    cloud = pv.PolyData(pts)
    cloud["deflection (mm)"] = np.concatenate([t["w"] for t in tops])
    p.add_mesh(cloud, scalars="deflection (mm)", point_size=6, render_points_as_spheres=True,
               cmap="viridis")
    p.add_text(f"row {row}: {mo.pct[row]:.0f} % stance", font_size=10)
    p.show_axes()
    p.view_xz()
    p.show()


if __name__ == "__main__":
    main()
