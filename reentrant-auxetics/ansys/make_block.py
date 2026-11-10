"""Lattice test block for ANSYS, as an exact 2D outline in DXF (mm).

Same block as unit_cell_setup.md: 4 x 4 cells, right edge closed, top verticals
mirrored from the bottom (180 deg symmetric), loading strips top and bottom.
Walls are the centrelines thickened to t (round joints), all merged into one
region. Open the DXF in SpaceClaim / Discovery and Fill it -> one surface.

Run:  python make_block.py            (forefoot cell -> block_fore.dxf / .step / .stl / .png)
      python make_block.py fore side  (strips on left/right for run B -> block_fore_side.*)
      python make_block.py mid        (midfoot / heel cell -> block_mid.dxf)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString, box
from shapely.ops import unary_union

CELLS = {"fore": dict(theta=30.0, l=4.0, h=8.0, t=1.6),
         "mid": dict(theta=30.0, l=5.0, h=10.0, t=1.0)}
NX, NY = 4, 4


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "fore"
    side = "side" in sys.argv[2:]          # strips left/right (run B, x1 loading)
    region, info = block(**CELLS[name], side=side)
    out = Path(__file__).parent / (f"block_{name}" + ("_side" if side else ""))
    write_dxf(region, out.with_suffix(".dxf"))
    write_step(region, out.with_suffix(".step"))
    import trimesh                          # 5 mm solid, for viewing / printing
    trimesh.creation.extrude_polygon(region, 5.0).export(out.with_suffix(".stl"))
    preview(region, out.with_suffix(".png"))
    print(f"{name}: {info}")
    print(f"  one region: {region.geom_type}, holes {len(region.interiors)}, "
          f"area {region.area:.1f} mm2 -> {out.name}.dxf / .png")


def block(theta, l, h, t, side=False):
    """Merged walls + strips. theta = re-entrant angle (positive, deg)."""
    lc, ls = l * np.cos(np.radians(theta)), l * np.sin(np.radians(theta))
    px, py = 2 * lc, 2 * (h - ls)
    P = {"A": (0, 0), "B": (0, h), "C": (lc, h - ls), "D": (px, h),
         "E": (lc, 2 * h - ls), "F": (0, py), "G": (px, py)}
    walls = [("A", "B"), ("B", "C"), ("C", "D"), ("C", "E"), ("E", "F"), ("E", "G")]

    def seg(a, b, i, j):
        return LineString([(P[a][0] + i * px, P[a][1] + j * py),
                           (P[b][0] + i * px, P[b][1] + j * py)])

    lines = [seg(a, b, i, j) for i in range(NX) for j in range(NY) for a, b in walls]
    lines += [seg("A", "B", NX, j) for j in range(NY)]                 # close right edge
    top0 = NY * py                                                     # valleys of top row
    lines += [LineString([(i * px, top0), (i * px, top0 + h)]) for i in range(NX + 1)]
    H = top0 + h                                                       # top of the walls
    W = NX * px
    walls = unary_union(lines).buffer(t / 2, quad_segs=8)
    if side:                                                           # run B: x1 loading
        strips = [box(-t / 2 - t, 0, -t / 2, H), box(W + t / 2, 0, W + t / 2 + t, H)]
    else:                                                              # runs A, C: x2
        strips = [box(-t / 2, -t, W + t / 2, 0), box(-t / 2, H, W + t / 2, H + t)]
    region = unary_union([walls, *strips])
    info = (f"theta -{theta:g}, l {l:g}, h {h:g}, t {t:g} | px {px:.3f}, py {py:.3f} | "
            f"block {W:.2f} wide, walls y 0..{H:.1f}, strips {t:g} thick")
    return region, info


def write_dxf(region, path):
    """Minimal ASCII DXF (R12): one closed POLYLINE per ring, units mm."""
    out = ["0", "SECTION", "2", "HEADER", "9", "$INSUNITS", "70", "4", "0", "ENDSEC",
           "0", "SECTION", "2", "ENTITIES"]
    for ring in [region.exterior, *region.interiors]:
        out += ["0", "POLYLINE", "8", "0", "66", "1", "70", "1"]
        for x, y in list(ring.coords)[:-1]:
            out += ["0", "VERTEX", "8", "0", "10", f"{x:.5f}", "20", f"{y:.5f}", "30", "0.0"]
        out += ["0", "SEQEND"]
    out += ["0", "ENDSEC", "0", "EOF"]
    path.write_text("\n".join(out) + "\n")


def write_step(region, path):
    """AP214 STEP: one planar face (z = 0) bounded by the outline and its holes,
    as a surface body (SHELL_BASED_SURFACE_MODEL), units mm. Arcs are kept as
    the short straight segments of the outline."""
    from shapely.geometry.polygon import orient
    region = orient(region, 1.0)              # outer CCW, holes CW (normal +z)
    L, n = [], [0]

    def add(s):
        n[0] += 1
        L.append(f"#{n[0]}={s};")
        return f"#{n[0]}"

    ctx = add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
    add(f"APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,{ctx})")
    pc = add(f"PRODUCT_CONTEXT('',{ctx},'mechanical')")
    prod = add(f"PRODUCT('{path.stem}','{path.stem}','',({pc}))")
    pdc = add(f"PRODUCT_DEFINITION_CONTEXT('part definition',{ctx},'design')")
    pdf = add(f"PRODUCT_DEFINITION_FORMATION('','',{prod})")
    pd = add(f"PRODUCT_DEFINITION('design','',{pdf},{pdc})")
    pds = add(f"PRODUCT_DEFINITION_SHAPE('','',{pd})")
    lu = add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
    au = add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
    su = add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
    unc = add(f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),{lu},'distance_accuracy_value','')")
    gctx = add(f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({unc}))"
               f"GLOBAL_UNIT_ASSIGNED_CONTEXT(({lu},{au},{su}))REPRESENTATION_CONTEXT('',''))")

    o = add("CARTESIAN_POINT('',(0.,0.,0.))")
    zd = add("DIRECTION('',(0.,0.,1.))")
    xd = add("DIRECTION('',(1.,0.,0.))")
    ax = add(f"AXIS2_PLACEMENT_3D('',{o},{zd},{xd})")
    plane = add(f"PLANE('',{ax})")

    bounds = []
    for k, ring in enumerate([region.exterior, *region.interiors]):
        pts = list(ring.coords)[:-1]
        vs = [add(f"VERTEX_POINT('',{add(f'CARTESIAN_POINT({chr(39)}{chr(39)},({x:.6f},{y:.6f},0.))')})")
              for x, y in pts]
        oes = []
        for i, (x0, y0) in enumerate(pts):
            x1, y1 = pts[(i + 1) % len(pts)]
            d = np.hypot(x1 - x0, y1 - y0)
            p0 = add(f"CARTESIAN_POINT('',({x0:.6f},{y0:.6f},0.))")
            dr = add(f"DIRECTION('',({(x1 - x0) / d:.9f},{(y1 - y0) / d:.9f},0.))")
            ln = add(f"LINE('',{p0},{add(f'VECTOR({chr(39)}{chr(39)},{dr},{d:.6f})')})")
            ec = add(f"EDGE_CURVE('',{vs[i]},{vs[(i + 1) % len(pts)]},{ln},.T.)")
            oes.append(add(f"ORIENTED_EDGE('',*,*,{ec},.T.)"))
        loop = add(f"EDGE_LOOP('',({','.join(oes)}))")
        kind = "FACE_OUTER_BOUND" if k == 0 else "FACE_BOUND"
        bounds.append(add(f"{kind}('',{loop},.T.)"))
    face = add(f"ADVANCED_FACE('',({','.join(bounds)}),{plane},.T.)")
    shell = add(f"OPEN_SHELL('',({face}))")
    sbsm = add(f"SHELL_BASED_SURFACE_MODEL('',({shell}))")
    rep = add(f"MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{path.stem}',({sbsm},{ax}),{gctx})")
    add(f"SHAPE_DEFINITION_REPRESENTATION({pds},{rep})")

    head = ["ISO-10303-21;", "HEADER;", "FILE_DESCRIPTION(('re-entrant lattice test block'),'2;1');",
            f"FILE_NAME('{path.name}','2026-11-08T00:00:00',(''),(''),'make_block.py','','');",
            "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));", "ENDSEC;", "DATA;"]
    path.write_text("\n".join(head + L + ["ENDSEC;", "END-ISO-10303-21;"]) + "\n")


def preview(region, path):
    fig, ax = plt.subplots(figsize=(4, 7))
    ax.fill(*region.exterior.xy, color="#d9c9b0")
    ax.plot(*region.exterior.xy, color="k", lw=0.5)
    for r in region.interiors:
        ax.fill(*r.xy, color="white")
        ax.plot(*r.xy, color="k", lw=0.5)
    ax.set_aspect("equal")
    ax.grid(lw=0.3)
    fig.savefig(path, dpi=130, bbox_inches="tight")


if __name__ == "__main__":
    main()
