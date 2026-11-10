"""Batch plane-stress grid search for ANSYS Mechanical APDL (run A for every cell).

For every candidate cell (the 24 of the analytical grid search + the two zone cells) this
writes one self-contained APDL input file:
  4 x 4 re-entrant block, ROUNDED joints and ends, as in the hand-run model, right edge closed,
  top mirrored from the bottom, loading strips top and bottom (thickness t);
  PLANE182 plane stress, 5 mm thick, Es 19.3 MPa, nu 0.45, mesh size t/4;
  bottom strip UY = 0 (+ one node UX = 0), top strip UY = -1 % of the wall height;
  E2 and nu21 measured on the middle cells from the four joint probes, as in run A.
Each run appends one line to work/grid_plane_results.csv (run_all.bat copies it up when done).

Run:  python make_grid_plane.py      -> inp/*.inp, grid_plane_results.csv (header), preview.png
Then: double-click run_all.bat       (runs every .inp with ANSYS in batch mode)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import LineString, box
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
GRID = HERE.parents[1] / "outputs" / "grid_eval.csv"
ZONES = HERE.parents[1] / "outputs" / "zone_pick.csv"
NX, NY, DEPTH, ES, NU = 4, 4, 5.0, 19.3, 0.45


def cells():
    """The 24 analytical-grid candidates plus the forefoot and midfoot/heel cells."""
    g = pd.read_csv(GRID)[["theta", "l", "h", "t", "E2", "nu21", "dip", "energy_J"]]
    z = pd.read_csv(ZONES, index_col=0)[["theta", "l", "h", "t", "E2", "nu21"]]
    z = z.loc[["fore", "mid"]].reset_index(names="zone")
    c = pd.concat([g.assign(zone=""), z], ignore_index=True)
    c = c.drop_duplicates(subset=["theta", "l", "h", "t"], keep="last").reset_index(drop=True)
    c.insert(0, "id", np.arange(1, len(c) + 1))
    return c


def block(theta, l, h, t):
    """Sharp-jointed block outline + the four probe joints and gauge lengths."""
    th = np.radians(abs(theta))                     # re-entrant angle, positive
    lc, ls = l * np.cos(th), l * np.sin(th)
    px, py = 2 * lc, 2 * (h - ls)
    P = {"A": (0, 0), "B": (0, h), "C": (lc, h - ls), "D": (px, h),
         "E": (lc, 2 * h - ls), "F": (0, py), "G": (px, py)}
    walls = [("A", "B"), ("B", "C"), ("C", "D"), ("C", "E"), ("E", "F"), ("E", "G")]

    def seg(a, b, i, j):
        return LineString([(P[a][0] + i * px, P[a][1] + j * py), (P[b][0] + i * px, P[b][1] + j * py)])

    lines = [seg(a, b, i, j) for i in range(NX) for j in range(NY) for a, b in walls]
    lines += [seg("A", "B", NX, j) for j in range(NY)]                       # close right edge
    top0 = NY * py
    lines += [LineString([(i * px, top0), (i * px, top0 + h)]) for i in range(NX + 1)]  # mirror top
    H, W = top0 + h, NX * px
    w = unary_union(lines).buffer(t / 2, quad_segs=8)   # rounded joints and ends, as the hand-run model
    region = unary_union([w, box(-t / 2, -t, W + t / 2, 0), box(-t / 2, H, W + t / 2, H + t)])
    region = orient(region.simplify(0.002), 1.0)
    if region.geom_type != "Polygon":
        region = max(region.geoms, key=lambda g: g.area)
    probes = {"P1": (lc + px, 2 * h - ls + py), "P2": (lc + px, 2 * h - ls + 2 * py),
              "P3": (lc + px, h - ls + 2 * py), "P4": (lc + 2 * px, h - ls + 2 * py)}
    return region, probes, px, py, W, H


def apdl(c, region, probes, px, py, W, H):
    """One self-contained APDL input file for run A of this cell."""
    t = c.t
    L = ["/BATCH", f"/TITLE, cell {c.id}: theta {c.theta:g} l {c.l:g} h {c.h:g} t {c.t:g}",
         "/PREP7", "ET,1,PLANE182", "KEYOPT,1,1,0", "KEYOPT,1,3,3", f"R,1,{DEPTH}",
         f"MP,EX,1,{ES}", f"MP,PRXY,1,{NU}"]
    k = ln = 0
    rings = [region.exterior, *region.interiors]
    ranges = []
    for ring in rings:
        pts = list(ring.coords)[:-1]
        k0, l0 = k + 1, ln + 1
        for x, y in pts:
            k += 1
            L.append(f"K,{k},{x:.5f},{y:.5f},0")
        for i in range(len(pts)):
            ln += 1
            L.append(f"L,{k0 + i},{k0 + (i + 1) % len(pts)}")
        ranges.append((l0, ln))
    for i, (a, b) in enumerate(ranges):
        L += [f"LSEL,S,LINE,,{a},{b}", "AL,ALL"]
    L += ["ALLSEL", f"ASEL,S,AREA,,2,{len(rings)}", "CM,HOLES,AREA", "ALLSEL", "ASBA,1,HOLES",
          "ALLSEL", f"ESIZE,{t / 4:.4f}", "MSHAPE,0,2D", "MSHKEY,0", "AMESH,ALL"]
    tol = 1e-3
    disp = 0.01 * H
    L += ["/SOLU", "ANTYPE,STATIC",
          f"NSEL,S,LOC,Y,{-t - tol:.5f},{-t + tol:.5f}", "D,ALL,UY,0",
          f"N0=NODE({-t / 2:.5f},{-t:.5f},0)", "NSEL,ALL", "D,N0,UX,0",
          f"NSEL,S,LOC,Y,{H + t - tol:.5f},{H + t + tol:.5f}", f"D,ALL,UY,{-disp:.6f}",
          "ALLSEL", "SOLVE", "FINISH", "/POST1", "SET,LAST"]
    for n, (x, y) in probes.items():
        L += [f"N{n}=NODE({x:.5f},{y:.5f},0)", f"*GET,UY{n},NODE,N{n},U,Y", f"*GET,UX{n},NODE,N{n},U,X"]
    L += [f"NSEL,S,LOC,Y,{-t - tol:.5f},{-t + tol:.5f}", "FSUM", "*GET,FY,FSUM,0,ITEM,FY", "ALLSEL",
          f"EPS2=(UYP2-UYP1)/{py:.6f}", f"EPS1=(UXP4-UXP3)/{px:.6f}",
          f"SIG=ABS(FY)/({W:.5f}*{DEPTH})", "E2=SIG/ABS(EPS2)", "NU21=-EPS1/EPS2",
          f"EPSG={disp / H:.6f}", "E2G=SIG/EPSG",
          "*GET,NNODE,NODE,0,COUNT",
          "*CFOPEN,grid_plane_results,csv,,APPEND",
          f"*VWRITE,{c.id},{c.theta},{c.l},{c.h},{c.t},E2,NU21,E2G,FY,EPS2,EPS1,NNODE",
          "(F4.0,',',F6.1,',',F5.2,',',F6.2,',',F5.2,',',G14.6,',',G14.6,',',G14.6,',',"
          "G14.6,',',G14.6,',',G14.6,',',F9.0)",
          "*CFCLOS", "FINISH"]
    return "\n".join(L) + "\n"


def preview(items, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(items)
    cols = 7
    fig, axs = plt.subplots(int(np.ceil(n / cols)), cols, figsize=(14, 2.9 * np.ceil(n / cols)))
    for ax in axs.flat:
        ax.axis("off")
    for ax, (c, region, probes) in zip(axs.flat, items):
        ax.fill(*region.exterior.xy, color="#d9c9b0", lw=0)
        for r in region.interiors:
            ax.fill(*r.xy, color="white", lw=0)
        ax.plot(*region.exterior.xy, color="k", lw=0.3)
        for x, y in probes.values():
            ax.plot(x, y, "o", ms=2, color="#c0392b")
        ax.set_title(f"{c.id}: θ{c.theta:g} l{c.l:g} h{c.h:g} t{c.t:g} {c.zone}", fontsize=7)
        ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def main():
    out = HERE / "inp"
    out.mkdir(exist_ok=True)
    for f in out.glob("*.inp"):
        f.unlink()
    (HERE / "work").mkdir(exist_ok=True)              # ANSYS runs here, so the CSV lives here too
    (HERE / "work" / "grid_plane_results.csv").write_text(
        "id,theta,l,h,t,E2_MPa,nu21,E2_global_MPa,FY_N,eps2,eps1,nodes\n")
    c = cells()
    c.to_csv(HERE / "grid_plane_cells.csv", index=False)
    items = []
    for _, r in c.iterrows():
        region, probes, px, py, W, H = block(r.theta, r.l, r.h, r.t)
        (out / f"cell_{int(r.id):02d}.inp").write_text(apdl(r, region, probes, px, py, W, H))
        items.append((r, region, probes))
        print(f"cell {int(r.id):2d}: theta {r.theta:4.0f} l {r.l} h {r.h:4.1f} t {r.t:.1f} {r.zone:4s} "
              f"holes {len(region.interiors):3d}  outline pts {sum(len(x.coords) for x in [region.exterior, *region.interiors])}")
    preview(items, HERE / "preview.png")
    print(f"{len(c)} input files in {out}")


if __name__ == "__main__":
    main()
