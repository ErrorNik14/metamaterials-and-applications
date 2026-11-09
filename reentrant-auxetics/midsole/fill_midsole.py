"""Fill the midsole volume with vertical re-entrant honeycomb slices.

Slices lie across the foot in the y-z plane (y = width, z = vertical) and
are spaced from heel to toe (x) with gaps. Optional bottom and top plates
join them into one body. Output: midsole_filled.stl

Run:  python fill_midsole.py
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

import re
from pathlib import Path

import numpy as np
import shapely
import trimesh
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

# ---- inputs (lattice: stiffer start pick, M&E E2 ~1.1 MPa, nu21 ~ -0.5) ----
STEP_FILE = INPUTS / "midsole_base.step"
OUT_FILE = OUTPUTS / "midsole_filled.stl"
LENGTH = 270.0          # midsole length after scaling (mm)
THETA = -30.0           # lattice angle (deg), negative = re-entrant
L = 4.0                 # inclined wall (mm)
H = 8.0                 # vertical wall (mm), h/l = 2
T = 1.4                 # wall thickness (mm), t/l = 0.35
B = 5.0                 # slice thickness along the foot (mm)
GAP = 3.0               # gap between slices (mm)
PLATE_BOT = 2.0         # bottom plate thickness joining the slices (mm), 0 = none
PLATE_TOP = 2.0         # top plate thickness, follows the footbed (mm), 0 = none
STRIP = 1.0             # top plate is built from strips this wide along x (mm)


def step_vertices(path):
    """Vertex coordinates and straight edges (vertex index pairs) of the STEP file.
    Works here because every face in the file is flat."""
    txt = path.read_text(encoding="latin-1")
    pts = {int(i): [float(v) for v in xyz.split(",")]
           for i, xyz in re.findall(r"#(\d+)=CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", txt)}
    vtx = re.findall(r"#(\d+)=VERTEX_POINT\('[^']*',#(\d+)\)", txt)
    index = {int(v): k for k, (v, _) in enumerate(vtx)}
    P = np.array([pts[int(p)] for _, p in vtx])
    edges = [(index[int(a)], index[int(b)])
             for a, b in re.findall(r"=EDGE_CURVE\('[^']*',#(\d+),#(\d+),", txt)]
    return P, np.array(edges)


def read_volume():
    """Bottom outline (polygon in x-y) and a top-height function z_top(x, y)."""
    P, edges = step_vertices(STEP_FILE)
    P = P - P.min(axis=0)
    P = P * LENGTH / P[:, 0].max()
    on_bot = P[:, 2] < 1e-6
    top = P[~on_bot]
    # outline = the edges lying in the flat bottom face, joined into a loop
    bot_edges = [LineString(P[[a, b], :2]) for a, b in edges if on_bot[a] and on_bot[b]]
    outline = max(shapely.polygonize(bot_edges).geoms, key=lambda g: g.area)
    lin = LinearNDInterpolator(top[:, :2], top[:, 2])
    near = NearestNDInterpolator(top[:, :2], top[:, 2])

    def z_top(x, y):
        z = lin(x, y)
        return np.where(np.isnan(z), near(x, y), z)

    return outline, z_top


def zone_lattice(x):
    """(theta, l, h, t) for the slice at x: from zone_pick.csv (zone_search.py)
    if it exists, else the single lattice above."""
    zp = OUTPUTS / "zone_pick.csv"
    if zp.exists():
        import pandas as pd
        bounds = {"fore": (0.0, 63.1), "mid": (63.1, 186.6), "heel": (186.6, 1e9)}
        z = pd.read_csv(zp, index_col=0)
        for name, (a, b) in bounds.items():
            if a <= x < b:
                r = z.loc[name]
                return r.theta, r.l, r.h, r.t
    return THETA, L, H, T


def lattice_walls(region, p=None):
    """Re-entrant walls (same cell as build_lat), clipped to the region, thickened."""
    THETA, L, H, T = p if p else (globals()["THETA"], globals()["L"], globals()["H"], globals()["T"])
    c, s = np.cos(np.radians(THETA)), np.sin(np.radians(THETA))
    px, py = 2 * L * c, 2 * (H + L * s)
    cell = [((0, 0), (0, H)),                                  # lower vertical
            ((0, H), (L * c, H + L * s)),                      # lower-left inclined
            ((L * c, H + L * s), (px, H)),                     # lower-right inclined
            ((L * c, H + L * s), (L * c, 2 * H + L * s)),      # middle vertical
            ((L * c, 2 * H + L * s), (0, py)),                 # upper-left inclined
            ((L * c, 2 * H + L * s), (px, py))]                # upper-right inclined
    x0, z0, x1, z1 = region.bounds
    lines = []
    for ox in np.arange(x0 - px, x1 + px, px):
        for oz in np.arange(z0, z1 + py, py):
            for (a, b) in cell:
                lines.append(LineString([(a[0] + ox, a[1] + oz), (b[0] + ox, b[1] + oz)]))
    walls = unary_union(lines).intersection(region)    # keep the parts inside
    return walls.buffer(T / 2)


def slice_solid(outline, z_top, x):
    """2D slice across the foot at length position x, in the (y, z) plane:
    lattice walls plus a skin around the edge."""
    cut = outline.intersection(LineString([(x, -1e3), (x, 1e3)]))
    _, ya, _, yb = cut.bounds
    pieces = []                      # the line could cross the outline in pieces
    for seg in getattr(cut, "geoms", [cut]):
        _, a, _, b = seg.bounds
        ys = np.linspace(a, b, 150)
        zs = z_top(np.full_like(ys, x), ys)
        pieces.append(Polygon([(a, 0), *zip(ys, zs), (b, 0)]).buffer(0))
    region = unary_union(pieces)
    p = zone_lattice(x)
    skin = region.difference(region.buffer(-p[3]))
    solid = unary_union([lattice_walls(region, p), skin]).intersection(region)
    return solid, ya, yb


def top_band(outline, z_top, x):
    """2D band of thickness PLATE_TOP under the footbed, across the foot at x."""
    cut = outline.intersection(LineString([(x, -1e3), (x, 1e3)]))
    bands = []
    for seg in getattr(cut, "geoms", [cut]):
        if seg.length < 1e-6:
            continue
        _, a, _, b = seg.bounds
        ys = np.linspace(a, b, 100)
        zs = z_top(np.full_like(ys, x), ys)
        bands.append(Polygon([*zip(ys, zs - PLATE_TOP), *zip(ys[::-1], zs[::-1])]).buffer(0))
    return unary_union(bands)


def plates(outline, z_top):
    """Bottom plate (flat) and top plate (follows the footbed, made of thin strips)."""
    parts = []
    if PLATE_BOT > 0:
        parts.append(trimesh.creation.extrude_polygon(outline, PLATE_BOT))
    if PLATE_TOP > 0:
        x_min, _, x_max, _ = outline.bounds
        for x in np.arange(x_min + STRIP / 2, x_max, STRIP):
            band = top_band(outline, z_top, x)
            if band.area > 1e-6:
                # strips overlap a little so the union merges them cleanly
                parts.append(extrude_yz(band, x - STRIP / 2 - 0.05, x + STRIP / 2 + 0.05))
    return parts


def extrude_yz(solid, x0, x1):
    """Extrude a 2D shape in (y, z) along x from x0 to x1."""
    # snap tiny gaps, else the mesh can leak; coarsen the snap until it is closed
    # last resort: close hairline slivers with a tiny buffer out-and-in first
    tries = [(g, solid) for g in (0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05)]
    tries += [(g, solid.buffer(0.02).buffer(-0.02)) for g in (0.001, 0.01, 0.02, 0.05)]
    tries += [(g, solid.simplify(0.01)) for g in (0.001, 0.002, 0.005)]
    for grid, shape in tries:
        clean = shapely.set_precision(shape, grid).buffer(0)
        polys = getattr(clean, "geoms", [clean])
        meshes = [trimesh.creation.extrude_polygon(p, x1 - x0) for p in polys if p.area > 1e-6]
        m = trimesh.util.concatenate(meshes)
        if m.is_watertight:
            break
    # (u, v, w) -> (x0 + w, u, v)
    m.apply_transform([[0, 0, 1, x0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    if m.volume < 0:
        m.invert()      # keep the faces pointing outward
    return m


def main():
    outline, z_top = read_volume()
    x_min, y_min, x_max, y_max = outline.bounds
    xs = np.arange(x_min + B / 2, x_max - B / 2 + 1e-9, B + GAP)
    print(f"midsole {x_max - x_min:.0f} x {y_max - y_min:.0f} mm, {len(xs)} slices")

    parts = []
    for x in xs:
        solid, ya, yb = slice_solid(outline, z_top, x)
        parts.append(extrude_yz(solid, x - B / 2, x + B / 2))
        print(f"  slice x = {x:5.1f} mm: y {ya:5.1f} to {yb:5.1f}, solid area {solid.area:6.0f} mm2")

    parts += plates(outline, z_top)
    print(f"joining {len(parts)} pieces (slices + plates) into one body ...")
    mesh = trimesh.boolean.union(parts, engine="manifold")
    mesh.export(OUT_FILE)
    print(f"wrote {OUT_FILE.name}: {len(mesh.faces)} triangles, "
          f"watertight {mesh.is_watertight}, volume {mesh.volume / 1e3:.1f} cm3")
    mesh.show()      # 3D preview window (close it to finish)


if __name__ == "__main__":
    main()
