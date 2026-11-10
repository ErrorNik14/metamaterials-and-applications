"""Figures for report_v2 (Part 1: re-entrant models, Part 2: midsole optimisation).

Run:  python make_figs_v2.py              (all; grid evaluation is cached)
      python make_figs_v2.py cell barrel  (only some, by name)
Output: report_v2/fig/
"""
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT / "midsole"), str(ROOT / "models"), str(ROOT / "report")]
from paths import ANSYS_RESULTS, OUTPUTS  # noqa: E402

FIG = HERE / "fig"
INK, MUTED, PAPER, ACCENT, CLAY, WARM, GRID = ("#1b1b1a", "#6b6a66", "#fcfcfb", "#c0392b",
                                              "#d9c9b0", "#e9e5de", "#d9d8d3")
HEX_C, RE_C = "#2a78d6", "#c0392b"


def style():
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 8.5,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.5, "axes.spines.top": False,
        "axes.spines.right": False, "axes.facecolor": PAPER, "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.4, "grid.linestyle": (0, (1, 3)), "legend.frameon": False})


def save(fig, name):
    fig.savefig(FIG / name, dpi=220, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print("  wrote", name)


def trim(path, pad=30):
    from PIL import Image, ImageChops
    im = Image.open(path).convert("RGB")
    bb = ImageChops.difference(im, Image.new("RGB", im.size, im.getpixel((2, 2)))).getbbox()
    if bb:
        im.crop((max(bb[0] - pad, 0), max(bb[1] - pad, 0), min(bb[2] + pad, im.width),
                 min(bb[3] + pad, im.height))).save(path)


# ---- geometry helpers ------------------------------------------------------------
def walls(theta, l, h, t, nx, ny, box_clip=True):
    """Thickened lattice walls; theta < 0 re-entrant, > 0 hexagonal (degrees)."""
    from shapely.geometry import LineString, box
    from shapely.ops import unary_union
    c, s = np.cos(np.radians(theta)), np.sin(np.radians(theta))
    px, py = 2 * l * c, 2 * (h + l * s)
    cell = [((0, 0), (0, h)), ((0, h), (l * c, h + l * s)), ((l * c, h + l * s), (px, h)),
            ((l * c, h + l * s), (l * c, 2 * h + l * s)), ((l * c, 2 * h + l * s), (0, py)),
            ((l * c, 2 * h + l * s), (px, py))]
    lines = [LineString([(a[0] + i * px, a[1] + j * py), (b[0] + i * px, b[1] + j * py)])
             for i in range(nx) for j in range(ny) for a, b in cell]
    lines += [LineString([(nx * px, j * py), (nx * px, j * py + h)]) for j in range(ny)]
    w = unary_union(lines).buffer(t / 2, quad_segs=6)
    if box_clip:
        w = w.intersection(box(-t / 2, 0, nx * px + t / 2, ny * py))
    return w, px, py


def draw_poly(ax, g, fc, ec=INK, lw=0.5):
    for p in getattr(g, "geoms", [g]):
        if not hasattr(p, "exterior"):
            continue
        ax.fill(*p.exterior.xy, color=fc, lw=0)
        ax.plot(*p.exterior.xy, color=ec, lw=lw)
        for r in p.interiors:
            ax.fill(*r.xy, color=PAPER, lw=0)
            ax.plot(*r.xy, color=ec, lw=lw)


# ---- Part 1 ------------------------------------------------------------------------
def cell():
    """Annotated unit cell, hexagonal vs re-entrant side by side."""
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for ax, th, title, col in ((axs[0], 30, "hexagonal  (θ > 0)", HEX_C),
                               (axs[1], -30, "re-entrant  (θ < 0)", RE_C)):
        l, h, t = 4.0, 8.0, 1.0
        w, px, py = walls(th, l, h, t, 2, 2)
        draw_poly(ax, w, "#ece9e2", lw=0.4)
        c, s = np.cos(np.radians(th)), np.sin(np.radians(th))
        # highlight one cell's walls
        pts = [(px, h), (px + l * c, h + l * s), (px + l * c, 2 * h + l * s)]
        ax.plot([px, px], [0, h], color=col, lw=2)
        ax.plot([pts[0][0], pts[1][0]], [pts[0][1], pts[1][1]], color=col, lw=2)
        ax.annotate("h", (px - 0.6, h / 2), ha="right", va="center", fontsize=9, color=col, style="italic")
        mx, my = (pts[0][0] + pts[1][0]) / 2, (pts[0][1] + pts[1][1]) / 2
        ax.annotate("l", (mx, my + (0.9 if th > 0 else -1.4)), ha="center", fontsize=9, color=col,
                    style="italic")
        # angle arc
        a0 = 0 if th > 0 else th
        arc = np.radians(np.linspace(a0, a0 + abs(th), 20))
        ax.plot(px + 1.8 * np.cos(arc), h + 1.8 * np.sin(arc), color=col, lw=0.7)
        ax.plot([px, px + 2.4], [h, h], color=col, lw=0.5, ls=(0, (2, 1.5)))
        ax.text(px + 2.6, h + (0.5 if th > 0 else -1.2), "θ", fontsize=9, color=col, style="italic")
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        ax.set_aspect("equal"); ax.axis("off")
    # axes arrows
    ax = axs[1]
    ax.annotate("", (-3, -2.5), (2, -2.5), arrowprops=dict(arrowstyle="<-", lw=0.6, color=MUTED))
    ax.annotate("", (-3, -2.5), (-3, 2.5), arrowprops=dict(arrowstyle="<-", lw=0.6, color=MUTED))
    ax.text(2.3, -2.9, "x₁", fontsize=8, color=MUTED); ax.text(-3.4, 2.8, "x₂", fontsize=8, color=MUTED)
    fig.tight_layout()
    save(fig, "unit_cell.png")


def timeline():
    """Short history of re-entrant auxetics, one line."""
    ev = [(1944, "Love — negative ν noted\nin some crystals"),
          (1982, "Gibson et al. — honeycomb\nmechanics, re-entrant cells"),
          (1987, "Lakes — re-entrant\nauxetic foam (Science)"),
          (1991, "Evans — the word\n'auxetic'"),
          (1996, "Masters & Evans —\nstretching + hinging"),
          (2004, "Ohtaki et al. —\nlarge deflection"),
          (2015, "Malek & Gibson —\nrigid joints, shear"),
          (2018, "Reviews: AM-printed\nauxetics, footwear")]
    fig, ax = plt.subplots(figsize=(7.4, 1.9))
    ax.axis("off")
    # events evenly spaced (not to scale) so labels never collide; a break marks the long gap
    xs = np.arange(len(ev), dtype=float)
    ax.plot([-0.4, len(ev) - 0.6], [0, 0], color=INK, lw=0.8)
    ax.plot([0.42, 0.52], [-0.08, 0.08], color=PAPER, lw=4)                # gap 1944 -> 1982
    ax.plot([0.40, 0.50], [-0.08, 0.08], color=INK, lw=0.6)
    ax.plot([0.46, 0.56], [-0.08, 0.08], color=INK, lw=0.6)
    for i, (x, (y, txt)) in enumerate(zip(xs, ev)):
        up = i % 2 == 0
        ax.plot([x, x], [0, 0.35 if up else -0.35], color=MUTED, lw=0.5)
        ax.plot(x, 0, "o", ms=4, color=ACCENT if y in (1987, 2004) else INK)
        ax.text(x, 0.42 if up else -0.42, f"{y}\n{txt}" if up else f"{txt}\n{y}", ha="center",
                va="bottom" if up else "top", fontsize=6.3, color=INK, linespacing=1.25)
    ax.set_xlim(-0.7, len(ev) - 0.3); ax.set_ylim(-1.3, 1.3)
    save(fig, "timeline.png")


def models():
    """Copy the model comparison figures (report mode, no own titles)."""
    import os, subprocess
    subprocess.run([sys.executable, str(ROOT / "models" / "model_compare.py")], cwd=OUTPUTS,
                   env={**os.environ, "REPORT": "1", "MPLBACKEND": "Agg"}, capture_output=True)
    for n in ("model_compare_moduli.png", "model_compare_curves.png"):
        shutil.copy(OUTPUTS / n, FIG / n)
    shutil.copy(ROOT / "ansys" / "geometry" / "block_fore.png", FIG / "ansys_block.png")
    shutil.copy(ROOT / "ansys" / "geometry" / "block_fore_side.png", FIG / "ansys_block_side.png")
    print("  copied model + ANSYS block figures")


def ansys_curve():
    """Full ANSYS run C including the collapse, vs Ohtaki."""
    from lattice_law import x2_law
    d = pd.read_csv(ANSYS_RESULTS / "x2_compression_full.csv")
    e, s = x2_law(-30.0, 8.0, 4.0, 1.6, 19.3)
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.plot(e * 100, s, color=INK, lw=1.6, label="Ohtaki")
    ax.plot(d.global_strain * 100, d.stress_MPa, "o-", color=ACCENT, ms=3, lw=0.8,
            label="ANSYS, whole block")
    k = d.stress_MPa.idxmax()
    ax.annotate("peak 0.36 MPa:\ncells buckle", (d.global_strain[k] * 100, d.stress_MPa[k]),
                xytext=(19, 0.62), fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))
    ax.set_xlim(0, 30); ax.set_ylim(0, 1.0)
    ax.set_xlabel("compression strain (%)"); ax.set_ylabel("stress (MPa)")
    ax.legend(fontsize=7, loc="upper left")
    save(fig, "ansys_collapse.png")


# ---- Part 2 ------------------------------------------------------------------------
def barrel():
    """3D: hexagonal block barrels out, re-entrant block draws in (same strain)."""
    import pyvista as pv
    import trimesh
    from shapely.affinity import translate
    pv.OFF_SCREEN = True
    eps, depth = 0.22, 14.0

    def block(theta, nu):
        w, px, py = walls(theta, 4.0, 8.0, 1.2, 4, 3)
        W, H = 4 * px, 3 * py
        g = translate(w, -W / 2, 0).segmentize(0.4)   # dense, so walls can curve
        from shapely.ops import transform

        def deform(x, y, z=None):
            x, y = np.asarray(x), np.asarray(y)
            prof = np.sin(np.pi * np.clip(y / H, 0, 1))          # ends held by plates
            return x * (1 + nu * eps * prof), y * (1 - eps)
        return g, transform(deform, g), W, H

    def mesh(poly, dz=0.0):
        m = trimesh.util.concatenate([trimesh.creation.extrude_polygon(q, depth)
                                      for q in getattr(poly, "geoms", [poly]) if q.area > 1e-3])
        m.apply_transform([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])  # y -> z
        return pv.wrap(m)

    from shapely.geometry import box
    p = pv.Plotter(off_screen=True, window_size=(3200, 1500), shape=(1, 2), border=False)
    for col, (theta, nu, name, c) in enumerate(((30, 1.6, "hexagonal", HEX_C),
                                                (-30, -1.0, "re-entrant", RE_C))):
        p.subplot(0, col)
        p.set_background(PAPER)
        g0, g1, W, H = block(theta, nu)
        for poly, dx, colr, op in ((g0, -0.75 * W, "#ece9e2", 1.0), (g1, 0.75 * W, CLAY, 1.0)):
            from shapely.affinity import translate as tr
            m = mesh(tr(poly, dx, 0))
            p.add_mesh(m, color=colr, ambient=0.3, diffuse=0.7, smooth_shading=False)
            hh = H if dx < 0 else H * (1 - eps)
            for zz in (-1.6, hh):
                plate = mesh(box(dx - W * 0.62, zz, dx + W * 0.62, zz + 1.6))
                p.add_mesh(plate, color="#8f8c84", ambient=0.3)
        # straight grey: original width; coloured curve: the deformed side profile
        zz = np.linspace(0, H, 60)
        prof = np.sin(np.pi * zz / H)
        for sx in (-1, 1):
            x0 = 0.75 * W + sx * W / 2
            p.add_mesh(pv.Line((x0, -2, 0), (x0, -2, H * (1 - eps))), color="#8f8c84", line_width=2)
            xs = 0.75 * W + sx * (W / 2 + 0.6) * (1 + nu * eps * prof)
            pts = np.c_[xs, np.full_like(zz, -2.0), zz * (1 - eps)]
            p.add_mesh(pv.lines_from_points(pts), color=c, line_width=5)
        p.camera_position = [(-150.0, -300, 170), (0.0, 0, 26), (0, 0, 1)]
        p.enable_parallel_projection()
        p.camera.parallel_scale = 50          # same scale for both panels
    p.screenshot(FIG / "barrel_3d.png")
    p.close()
    trim(FIG / "barrel_3d.png")
    print("  wrote barrel_3d.png")


def pipeline():
    """Flow diagram of the midsole optimisation."""
    steps = [("Bulk midsole\n(STEP)", "fill_midsole"), ("Lattice slices\n+ plates", "fill_midsole"),
             ("Fukuchi stance\n(force, COP, markers)", "foot_motion"),
             ("Posed\npseudo-foot", "foot_motion"), ("Beam on lattice\nfoundation", "foot_conform"),
             ("Sink from\nforce balance", "full_stance"), ("Dip, energy,\npeak strain", "grid search")]
    fig, ax = plt.subplots(figsize=(7.4, 1.6))
    ax.axis("off")
    xs = np.linspace(0, 1, len(steps))
    for i, (x, (txt, src)) in enumerate(zip(xs, steps)):
        ax.add_patch(plt.Rectangle((x - 0.062, 0.25), 0.124, 0.5, transform=ax.transAxes,
                                   fc=PAPER if i < 6 else "#f6e3df", ec=INK, lw=0.6))
        ax.text(x, 0.5, txt, ha="center", va="center", fontsize=6.6, transform=ax.transAxes)
        ax.text(x, 0.12, src + (".py" if "grid" not in src else ""), ha="center", fontsize=5.6,
                color=MUTED, family="monospace", transform=ax.transAxes)
        if i < len(steps) - 1:
            ax.annotate("", (xs[i + 1] - 0.065, 0.5), (x + 0.065, 0.5), xycoords="axes fraction",
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK))
    ax.set_xlim(-0.08, 1.08)
    save(fig, "pipeline.png")


# ---- grid evaluation: dip and energy per candidate -------------------------------------
GRID_CSV = OUTPUTS / "grid_eval.csv"
STEP = 4


def candidates():
    import zone_search as zs
    c = zs.candidates()
    c["tl"] = c.t / c.l
    # 24 spread over stiffness, keep a mix of angles
    idx = np.unique(np.linspace(0, len(c) - 1, 24).astype(int))
    return c.sort_values("E2").iloc[idx].reset_index(drop=True)


def grid_eval():
    """Uniform lattice over the whole midsole, every STEP-th stance row:
    dip = largest sink; energy = work done while loading (integral F d sink up to max)."""
    if GRID_CSV.exists():
        print("  grid_eval cached:", GRID_CSV.name)
        return pd.read_csv(GRID_CSV)
    from fill_midsole import read_volume
    from foot_conform import ES, foot_bodies
    from full_stance import F_MIN, sink_solve
    from lattice_law import x2_law
    mo = pd.read_csv(OUTPUTS / "foot_motion.csv")
    feet, (o, zt) = foot_bodies(), read_volume()
    rows = []
    for i, c in candidates().iterrows():
        try:
            law = x2_law(c.theta, c.h, c.l, c.t, ES)
            if len(law[0]) < 5:
                raise ValueError("short Ohtaki curve")
        except Exception as ex:
            print("  skip", dict(c), ex); continue
        F, S, P = [], [], []
        for k in range(0, len(mo), STEP):
            r = mo.loc[k]
            if r.F < F_MIN:
                continue
            s, tops = sink_solve(feet, o, zt, law, r)
            F.append(r.F); S.append(s); P.append(max(t["eps"].max() for t in tops))
        F, S = np.array(F), np.array(S)
        kmax = int(np.argmax(S))
        energy = np.trapezoid(F[:kmax + 1], S[:kmax + 1]) / 1000.0     # J (N*mm -> J)
        rows.append({**c.to_dict(), "dip": S.max(), "energy_J": energy, "peak_strain": max(P)})
        print(f"  {i:2d} th {c.theta:5.0f} l {c.l} h {c.h:4.1f} t {c.t:.1f}  E2 {c.E2:.2f}  "
              f"dip {S.max():5.2f} mm  energy {energy:.3f} J  strain {max(P) * 100:4.1f} %")
    d = pd.DataFrame(rows)
    d.to_csv(GRID_CSV, index=False)
    return d


def pareto_mask(d):
    """Non-dominated among feasible cells: lower dip and higher energy are both better."""
    m = ~d.bottom.to_numpy()
    for i in np.where(m)[0]:
        better = ~d.bottom & (d.dip <= d.dip[i]) & (d.energy_J >= d.energy_J[i]) & \
                 ((d.dip < d.dip[i]) | (d.energy_J > d.energy_J[i]))
        m[i] = not better.any()
    return m


def grid_plot():
    d = grid_eval()
    d["bottom"] = d.dip >= 19.99                        # hit the sink limit: bottoms out
    d["pareto"] = pareto_mask(d)
    ex = example_ids(d)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axs[0]
    ok = d[~d.bottom]
    sc = ax.scatter(ok.dip, ok.energy_J, c=ok.tl, cmap="viridis", s=26, edgecolor=PAPER, lw=0.6, zorder=3)
    bo = d[d.bottom]
    ax.scatter(bo.dip, bo.energy_J, marker="x", color=MUTED, s=18, lw=0.8, zorder=3)
    ax.text(19.5, bo.energy_J.min() - 0.1, "bottomed out (20 mm limit) ×", fontsize=6.5,
            color=MUTED, ha="right", va="center")
    f = d[d.pareto].sort_values("dip")
    ax.plot(f.dip, f.energy_J, color=INK, lw=0.8, zorder=2)
    for lab, i in ex.items():
        ax.annotate(lab, (d.dip[i], d.energy_J[i]), xytext=(5, -10), textcoords="offset points",
                    fontsize=8, weight="bold", color=ACCENT)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("t / l", fontsize=7); cb.outline.set_linewidth(0.3); cb.ax.tick_params(labelsize=6.5)
    ax.set_xlabel("dip: largest sink over the stance (mm)")
    ax.set_ylabel("energy absorbed in loading (J)")
    ax.set_title("Trade-off: dip vs energy", fontsize=9, loc="left")
    ax = axs[1]
    ax.scatter(ok.E2, ok.peak_strain * 100, c=ok.tl, cmap="viridis", s=20, edgecolor=PAPER, lw=0.5)
    ax.scatter(bo.E2, bo.peak_strain * 100, marker="x", color=MUTED, s=16, lw=0.8)
    ax.axhspan(15, 20, color=GRID, alpha=0.5, lw=0)
    ax.text(0.05, 17.5, "zoning target", fontsize=6.5, color=MUTED, va="center")
    for lab, i in ex.items():
        ax.annotate(lab, (d.E2[i], d.peak_strain[i] * 100), xytext=(4, 3), textcoords="offset points",
                    fontsize=8, weight="bold", color=ACCENT)
    ax.set_xscale("log")
    ax.set_xlabel("E₂ (M&E, MPa)"); ax.set_ylabel("peak column strain (%)")
    ax.set_title("Stiffness vs strain", fontsize=9, loc="left")
    fig.tight_layout()
    save(fig, "grid_pareto.png")


def example_ids(d):
    """Six labelled examples spread along the trade-off."""
    f = d[d.dip < 19.99] if "dip" in d else d
    o = f.sort_values("E2").index.to_numpy()
    pick = o[np.linspace(0, len(o) - 1, 6).astype(int)]
    return {chr(65 + k): int(i) for k, i in enumerate(pick)}


def grid_examples():
    """Rendered forefoot slice for six candidates, each with a zoomed lattice inset."""
    import pyvista as pv
    import trimesh
    from PIL import Image, ImageDraw, ImageFont
    from shapely.geometry import LineString, box
    from shapely.ops import unary_union
    import fill_midsole as fm
    pv.OFF_SCREEN = True
    d = grid_eval()
    ex = example_ids(d)
    o, z = fm.read_volume()
    x = 34.5
    cut = o.intersection(LineString([(x, -1e3), (x, 1e3)]))
    _, ya, _, yb = cut.bounds
    ys = np.linspace(ya, yb, 150)
    zs = z(np.full_like(ys, x), ys)
    from shapely.geometry import Polygon
    region = Polygon([(ya, 0), *zip(ys, zs), (yb, 0)]).buffer(0)
    tiles = []
    for lab, i in ex.items():
        c = d.loc[i]
        p = (c.theta, c.l, c.h, c.t)
        solid = unary_union([fm.lattice_walls(region, p), region.difference(region.buffer(-c.t))]).intersection(region)
        m = fm.extrude_yz(solid, 0, fm.B * 3)
        pl = pv.Plotter(off_screen=True, window_size=(1400, 900))
        pl.set_background(PAPER)
        pl.add_mesh(pv.wrap(m), color=CLAY, ambient=0.3, diffuse=0.75)
        pl.enable_parallel_projection()
        pl.camera_position = [(-260, (ya + yb) / 2 - 90, 60), (7, (ya + yb) / 2, 20), (0, 0, 1)]
        pl.camera.parallel_scale = 34
        f = FIG / f"_ex_{lab}.png"
        pl.screenshot(f); pl.close()
        im = Image.open(f).convert("RGB")
        # zoom: front face, centre
        W, H = im.size
        box_ = (int(W * 0.40), int(H * 0.38), int(W * 0.58), int(H * 0.62))
        zoom = im.crop(box_).resize((470, 627))
        ImageDraw.Draw(im).rectangle(box_, outline=ACCENT, width=4)
        from PIL import ImageChops
        bb = ImageChops.difference(im, Image.new("RGB", im.size, im.getpixel((2, 2)))).getbbox()
        body = im.crop((bb[0] - 20, bb[1] - 20, bb[2] + 20, bb[3] + 20))
        # tile: [label | render | zoom] above two text lines; nothing overlaps
        tw, th = body.width + 40 + zoom.width + 20, max(body.height, zoom.height) + 190
        tile = Image.new("RGB", (tw, th), PAPER)
        tile.paste(body, (0, 70 + (zoom.height - body.height) // 2))
        zx = body.width + 40
        tile.paste(zoom, (zx, 70))
        dr = ImageDraw.Draw(tile)
        dr.rectangle([zx, 70, zx + zoom.width - 1, 70 + zoom.height - 1], outline=ACCENT, width=4)
        try:
            fnt = ImageFont.truetype("arial.ttf", 48); fs = ImageFont.truetype("arial.ttf", 32)
        except OSError:
            fnt = fs = ImageFont.load_default()
        dr.text((10, 0), lab, fill=ACCENT, font=fnt)
        dr.text((10, th - 110), f"θ {c.theta:g}°   l {c.l:g}   h {c.h:g}   t {c.t:g} mm", fill=INK, font=fs)
        dr.text((10, th - 60), f"E2 {c.E2:.2f} MPa    dip {c.dip:.1f} mm    energy {c.energy_J:.2f} J",
                fill=MUTED, font=fs)
        tiles.append(tile); f.unlink()
    W = max(t.width for t in tiles); H = max(t.height for t in tiles)
    sheet = Image.new("RGB", (W * 2 + 80, H * 3 + 100), PAPER)       # 2 columns x 3 rows
    for k, im in enumerate(tiles):
        sheet.paste(im, ((k % 2) * (W + 80), (k // 2) * (H + 50)))
    sheet.save(FIG / "grid_examples.png")
    print("  wrote grid_examples.png")


ZONES = [("Forefoot", "fore", 0.0, 63.1, 34.5, "#c0392b"),
         ("Midfoot", "mid", 63.1, 186.6, 122.5, "#d98a73"),
         ("Heel", "heel", 186.6, 300.0, 226.5, "#e9c3b5")]


def zones():
    """Top view of the zoned midsole with each zone's 3D slice placed directly below
    its position (leader lines), plus a zoom of the lattice and the cell parameters."""
    import pyvista as pv
    from matplotlib.patches import ConnectionPatch, Rectangle
    from PIL import Image, ImageChops
    from shapely.geometry import LineString, Polygon, box
    from shapely.ops import unary_union
    import fill_midsole as fm
    pv.OFF_SCREEN = True
    pick = pd.read_csv(OUTPUTS / "zone_pick.csv", index_col=0)
    o, z = fm.read_volume()

    # ---- render one slice per zone (body + zoom crop)
    renders = {}
    for name, key, a, b, xs, col in ZONES:
        p = pick.loc[key]
        cut = o.intersection(LineString([(xs, -1e3), (xs, 1e3)]))
        _, ya, _, yb = cut.bounds
        ys = np.linspace(ya, yb, 150)
        zs = z(np.full_like(ys, xs), ys)
        region = Polygon([(ya, 0), *zip(ys, zs), (yb, 0)]).buffer(0)
        solid = unary_union([fm.lattice_walls(region, (p.theta, p.l, p.h, p.t)),
                             region.difference(region.buffer(-p.t))]).intersection(region)
        m = fm.extrude_yz(solid, 0, fm.B * 3)
        pl = pv.Plotter(off_screen=True, window_size=(1500, 1000))
        pl.set_background(PAPER)
        pl.add_mesh(pv.wrap(m), color=CLAY, ambient=0.3, diffuse=0.75)
        pl.enable_parallel_projection()
        zc = float(np.max(zs)) / 2
        pl.camera_position = [(-260, (ya + yb) / 2 - 90, zc + 40), (7, (ya + yb) / 2, zc), (0, 0, 1)]
        pl.camera.parallel_scale = 40
        img = np.asarray(pl.screenshot(None, return_img=True)); pl.close()
        im = Image.fromarray(img).convert("RGB")
        W, H = im.size
        bx = (int(W * 0.40), int(H * 0.36), int(W * 0.58), int(H * 0.62))
        zoom = np.asarray(im.crop(bx))
        bb = ImageChops.difference(im, Image.new("RGB", im.size, im.getpixel((2, 2)))).getbbox()
        crop = (bb[0] - 15, bb[1] - 15, bb[2] + 15, bb[3] + 15)
        body = np.asarray(im.crop(crop))
        rect = (bx[0] - crop[0], bx[1] - crop[1], bx[2] - bx[0], bx[3] - bx[1])   # zoom box in body px
        renders[key] = (body, zoom, rect)

    # ---- compose
    fig = plt.figure(figsize=(7.4, 6.6))
    top = fig.add_axes([0.06, 0.62, 0.88, 0.36])
    top.axis("off"); top.set_aspect("equal")
    for name, key, a, b, xs, col in ZONES:
        part = o.intersection(box(a, -1e3, b, 1e3))
        for g in getattr(part, "geoms", [part]):
            top.fill(*g.exterior.xy, color=col, alpha=0.35, lw=0)
        top.text((a + min(b, 270)) / 2, o.bounds[3] + 5, name.upper(), ha="center", fontsize=7,
                 color=INK, family="sans-serif")
    top.plot(*o.exterior.xy, color=INK, lw=0.7)
    for hx in (63.1, 186.6):
        top.plot([hx, hx], [o.bounds[1] - 2, o.bounds[3] + 2], color=INK, lw=0.5, ls=(0, (3, 2)))
    top.text(-2, o.bounds[1] - 9, "toe", fontsize=6.5, color=MUTED)
    top.text(272, o.bounds[1] - 9, "heel", fontsize=6.5, color=MUTED, ha="right")
    top.set_xlim(-8, 278); top.set_ylim(o.bounds[1] - 14, o.bounds[3] + 12)
    fig.canvas.draw()

    colw = 0.29
    centres = {"fore": 0.17, "mid": 0.50, "heel": 0.83}   # evenly spaced; leaders angle to them
    for name, key, a, b, xs, col in ZONES:
        p = pick.loc[key]
        body, zoom, rect = renders[key]
        cut = o.intersection(LineString([(xs, -1e3), (xs, 1e3)]))
        _, ya, _, yb = cut.bounds
        top.plot([xs, xs], [ya, yb], color=ACCENT, lw=1.6)
        top.plot(xs, ya, "o", ms=3, color=ACCENT)
        # column centred under the slice, kept on the page
        x0 = centres[key] - colw / 2
        axb = fig.add_axes([x0, 0.33, colw, 0.22])
        axb.imshow(body); axb.axis("off")
        axb.add_patch(Rectangle(rect[:2], rect[2], rect[3], fill=False, ec=ACCENT, lw=0.8))
        axb.set_title(f"{name}  ·  x {xs:g} mm", fontsize=7.5, color=INK, loc="left", pad=2)
        fig.add_artist(ConnectionPatch((xs, ya), (0.5, 1.12), "data", "axes fraction",
                                       axesA=top, axesB=axb, color=ACCENT, lw=0.6,
                                       ls=(0, (3, 2))))
        axz = fig.add_axes([x0, 0.04, colw * 0.52, 0.26])
        axz.imshow(zoom); axz.set_xticks([]); axz.set_yticks([])
        for sp in axz.spines.values():
            sp.set_visible(True); sp.set_color(ACCENT); sp.set_linewidth(1.0)
        fig.add_artist(ConnectionPatch((rect[0] + rect[2] / 2, rect[1] + rect[3]), (0.5, 1.0),
                                       "data", "axes fraction", axesA=axb, axesB=axz,
                                       color=ACCENT, lw=0.5))
        txt = "\n".join([f"θ  {p.theta:g}°", f"l  {p.l:g} mm", f"h  {p.h:g} mm", f"t  {p.t:g} mm",
                         f"t/l  {p.t / p.l:.2f}", "", f"E₂  {p.E2:.2f} MPa", f"ν₂₁  {p.nu21:.2f}",
                         f"peak  {p.peak_strain * 100:.1f} %"])
        fig.text(x0 + colw * 0.58, 0.29, txt, fontsize=6.6, va="top", color=INK, linespacing=1.45)
    fig.savefig(FIG / "zones_detail.png", dpi=260, facecolor=PAPER)
    plt.close(fig)
    print("  wrote zones_detail.png")


def foot_views(row=30):
    """Fig. 15: side cross-section (top) and a line-only isometric view (bottom)
    of the posed pseudo-foot on the midsole, same hairline style."""
    import trimesh
    import fill_midsole as fm
    import foot_anim as fa
    import foot_conform as fc
    mo = pd.read_csv(OUTPUTS / "foot_motion.csv")
    o, z = fm.read_volume()
    feet = fc.pose(fc.foot_bodies(), mo.loc[row])
    cz, sz = np.cos(np.radians(30)), np.sin(np.radians(30))

    def iso(p):
        p = np.atleast_2d(p)
        x, y, zz = p[:, 0], p[:, 1], p[:, 2]
        return np.c_[(x - y) * cz, (x + y) * sz * 0.55 + zz]   # flattened isometric

    fig = plt.figure(figsize=(7.4, 6.4))
    # ---- top: side section (reuse the first report's drawing routine, black lines)
    sh = fa.sections()
    ax = fig.add_axes([0.02, 0.56, 0.96, 0.40]); ax.axis("off"); ax.set_aspect("equal")
    for e in sh["sole"]:
        ax.plot(e[:, 0], e[:, 1], color=INK, lw=0.25)
    rear, toe, shank = fa.pose(mo, row)
    for nm, a in (("toe", toe), ("mid", rear), ("heel", rear)):
        for e in sh[nm]:
            q = fa.rot(e, fa.MTP, a)
            ax.plot(q[:, 0], q[:, 1], color=INK, lw=0.7)
    an = fa.rot(fa.ANKLE, fa.MTP, rear)
    d = np.array([-np.cos(np.radians(shank)), np.sin(np.radians(shank))])
    ax.plot([an[0], an[0] + 110 * d[0]], [an[1], an[1] + 110 * d[1]], color=INK, lw=0.7)
    for pnt, lab, off in ((an, "ankle", (40, 40)), (fa.rot(fa.HEEL_MID, fa.MTP, rear), "heel–midfoot", (10, 55)),
                          (fa.MTP, "MTP", (-10, 50))):
        ax.plot(*pnt, "o", ms=3.5, color=ACCENT)
        ax.annotate(lab, pnt, xytext=off, textcoords="offset points", fontsize=7, color=ACCENT, ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=ACCENT))
    ax.plot([-10, 285], [0, 0], color=INK, lw=0.6)
    ax.text(-8, 195, "SIDE SECTION  ·  y = 60 mm", fontsize=7, color=MUTED)

    # ---- bottom: isometric, lines only
    ax = fig.add_axes([0.02, 0.0, 0.96, 0.55]); ax.axis("off"); ax.set_aspect("equal")
    ring = np.asarray(o.segmentize(3.0).exterior.coords)
    ztop = z(ring[:, 0], ring[:, 1])
    bot, top = iso(np.c_[ring, np.zeros(len(ring))]), iso(np.c_[ring, ztop])
    ax.plot(*bot.T, color=INK, lw=0.5); ax.plot(*top.T, color=INK, lw=0.6)
    for i in range(0, len(ring), 6):                                   # vertical hairlines
        seg = iso(np.array([[*ring[i], 0], [*ring[i], ztop[i]]]))
        ax.plot(*seg.T, color=INK, lw=0.2)
    x0, _, x1, _ = o.bounds
    for xs in np.arange(x0 + fm.B / 2, x1, fm.B + fm.GAP)[::2]:        # slice top profiles
        from shapely.geometry import LineString
        cut = o.intersection(LineString([(xs, -1e3), (xs, 1e3)]))
        _, ya, _, yb = cut.bounds
        ys = np.linspace(ya, yb, 40)
        ax.plot(*iso(np.c_[np.full_like(ys, xs), ys, z(np.full_like(ys, xs), ys)]).T, color=INK, lw=0.2)
    for nm, m in feet.items():
        if nm == "heel":                                               # sphere: outline + equator
            c = m.vertices.mean(axis=0); r = (m.bounds[1, 2] - m.bounds[0, 2]) / 2
            t = np.linspace(0, 2 * np.pi, 120)
            cc = iso(c)[0]
            ax.fill(cc[0] + r * np.cos(t), cc[1] + r * np.sin(t), color=PAPER, lw=0, zorder=2)
            ax.plot(cc[0] + r * np.cos(t), cc[1] + r * np.sin(t), color=INK, lw=0.7, zorder=3)
            ax.plot(*iso(c + np.c_[r * np.cos(t), r * np.sin(t), 0 * t]).T, color=INK, lw=0.3, zorder=3)
            continue
        for f in m.faces:                                              # solid: hides the midsole lines
            ax.fill(*iso(m.vertices[f]).T, color=PAPER, lw=0, zorder=2)
        ang = m.face_adjacency_angles
        for a, b in m.face_adjacency_edges[ang > np.radians(20)]:
            ax.plot(*iso(m.vertices[[a, b]]).T, color=INK, lw=0.7, zorder=3)
    # shank line and hinges in 3D at y = 60
    for (px, pz) in (an, fa.rot(fa.HEEL_MID, fa.MTP, rear), fa.MTP):
        ax.plot(*iso([px, 60, pz])[0], "o", ms=3.5, color=ACCENT, zorder=4)
    s0, s1 = iso([an[0], 60, an[1]])[0], iso([an[0] + 110 * d[0], 60, an[1] + 110 * d[1]])[0]
    ax.plot([s0[0], s1[0]], [s0[1], s1[1]], color=INK, lw=0.7, zorder=3)
    ax.text(*iso([0, 0, 0])[0] + np.array([-20, -14]), "ISOMETRIC", fontsize=7, color=MUTED)
    fig.savefig(FIG / "foot_views.png", dpi=260, facecolor=PAPER)
    plt.close(fig)
    trim(FIG / "foot_views.png")
    print("  wrote foot_views.png")


def bulk():
    """3D render of the bulk midsole volume (midsole_base.step, scaled to 270 mm)."""
    import pyvista as pv
    import trimesh
    import fill_midsole as fm
    pv.OFF_SCREEN = True
    o, z = fm.read_volume()
    o = o.segmentize(2.0).simplify(0.05)
    # top surface: triangulate the outline densely, lift vertices to z_top
    from shapely.geometry import MultiPoint
    from shapely.ops import triangulate
    xs, ys = np.meshgrid(np.arange(o.bounds[0], o.bounds[2], 3.0), np.arange(o.bounds[1], o.bounds[3], 3.0))
    pts = [p for p in MultiPoint(np.c_[xs.ravel(), ys.ravel()]).geoms if o.buffer(-1).contains(p)]
    ring = np.asarray(o.exterior.coords)[:-1]
    cloud = MultiPoint([*pts, *map(tuple, ring)])
    tris = [t for t in triangulate(cloud) if o.buffer(0.01).contains(t)]
    V, F, idx = [], [], {}
    def vid(x, y, top):
        k = (round(x, 4), round(y, 4), top)
        if k not in idx:
            idx[k] = len(V)
            V.append((x, y, float(z(np.array([x]), np.array([y]))[0]) if top else 0.0))
        return idx[k]
    for t in tris:
        (a, b, c) = list(t.exterior.coords)[:3]
        F.append([vid(*a, True), vid(*b, True), vid(*c, True)])
        F.append([vid(*c, False), vid(*b, False), vid(*a, False)])
    for i in range(len(ring)):
        (x0, y0), (x1, y1) = ring[i], ring[(i + 1) % len(ring)]
        a, b = vid(x0, y0, True), vid(x1, y1, True)
        c, d = vid(x1, y1, False), vid(x0, y0, False)
        F += [[a, d, c], [a, c, b]]
    m = trimesh.Trimesh(np.array(V), np.array(F), process=True)
    m.fix_normals()
    p = pv.Plotter(off_screen=True, window_size=(2600, 1500))
    p.set_background(PAPER)
    p.add_mesh(pv.wrap(m), color=CLAY, ambient=0.3, diffuse=0.75, specular=0.05, smooth_shading=False)
    p.add_mesh(pv.wrap(m).extract_feature_edges(35), color="#8a7b66", line_width=2)
    p.enable_parallel_projection()
    p.camera_position = [(-260, -330, 260), (135, 54, 18), (0, 0, 1)]
    p.camera.parallel_scale = 105
    p.enable_eye_dome_lighting()
    p.screenshot(FIG / "bulk_render.png"); p.close()
    trim(FIG / "bulk_render.png")
    print("  wrote bulk_render.png")


def filled():
    """3D render of the filled (zoned) midsole, midsole_filled.stl."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    m = pv.read(OUTPUTS / "midsole_filled.stl")
    p = pv.Plotter(off_screen=True, window_size=(3000, 1700))
    p.set_background(PAPER)
    p.add_mesh(m, color=CLAY, ambient=0.3, diffuse=0.75, specular=0.05, smooth_shading=False)
    p.enable_parallel_projection()
    p.camera_position = [(-300, -330, 230), (135, 54, 18), (0, 0, 1)]   # from the toe, outer side
    p.camera.parallel_scale = 100
    p.enable_eye_dome_lighting()
    p.screenshot(FIG / "filled_render.png"); p.close()
    trim(FIG / "filled_render.png")
    print("  wrote filled_render.png")


def reuse():
    """Midsole figures from the first report."""
    src = FIG  # reused figures already live in report/fig
    for n in ("film_dark.png", "loadpath.png", "strain_map.png", "sink.png", "exploded_render.png",
              "section_render.png", "cover_side.png", "cover_stack.png", "cover_slice.png",
              "stance_dark.png", "grid_dark.png"):
        shutil.copy(src / n, FIG / n)
    print("  copied midsole figures")


ALL = {"cell": cell, "timeline": timeline, "models": models, "ansys": ansys_curve, "barrel": barrel,
       "pipeline": pipeline, "grid": grid_plot, "examples": grid_examples, "zones": zones, "bulk": bulk, "filled": filled, "footviews": foot_views, "reuse": reuse}

if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    style()
    for name in (sys.argv[1:] or ALL):
        print(name)
        ALL[name]()
