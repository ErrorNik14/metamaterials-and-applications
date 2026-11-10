"""All report images, in the portfolio style. Output: report/fig/

Run:  python make_figs.py            (everything, ~3-4 min)
      python make_figs.py film cover (only some, by name)
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
MS = OUTPUTS                      # data the scripts produced
FIG = HERE / "fig"
sys.path[:0] = [str(_ROOT / "midsole"), str(_ROOT / "models")]

import fill_midsole as fm                                    # noqa: E402
import foot_anim as fa                                       # noqa: E402

INK, MUTED, PAPER, DARK, ACCENT = "#1b1b1a", "#6b6a66", "#fcfcfb", "#141413", "#c0392b"
CLAY, WARM = "#d9c9b0", "#e9e5de"
PASTEL = {"ground": "#7fb7b2", "foot": "#efc6c0", "plate": "#ecd9b0", "lattice": "#cfc4e6",
          "base": "#e3ddd2"}


def serif():
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "savefig.transparent": False})


def trim(path, pad=40):
    from PIL import Image, ImageChops
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, im.getpixel((2, 2)))
    bb = ImageChops.difference(im, bg).getbbox()
    if bb:
        im.crop((max(bb[0] - pad, 0), max(bb[1] - pad, 0), min(bb[2] + pad, im.width),
                 min(bb[3] + pad, im.height))).save(path)


def save(fig, name, dark=False, transparent=False):
    fig.savefig(FIG / name, dpi=220, facecolor="none" if transparent else (DARK if dark else PAPER),
                transparent=transparent, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("  wrote", name)


def blank(w, h, dark=False):
    fig, ax = plt.subplots(figsize=(w, h))
    fig.patch.set_facecolor(DARK if dark else PAPER)
    ax.set_facecolor(DARK if dark else PAPER)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def outline(ax, geom, color, lw, fill=None, hole=None):
    for g in getattr(geom, "geoms", [geom]):
        if g.is_empty or not hasattr(g, "exterior"):
            continue
        if fill:
            ax.fill(*g.exterior.xy, color=fill, lw=0)
        ax.plot(*g.exterior.xy, color=color, lw=lw)
        for r in g.interiors:
            if fill:
                ax.fill(*r.xy, color=hole or ax.get_facecolor(), lw=0)
            ax.plot(*r.xy, color=color, lw=lw)


# ---- shared geometry ----------------------------------------------------------
_vol = None


def volume():
    global _vol
    if _vol is None:
        _vol = fm.read_volume()
    return _vol


def slice_xs():
    o, _ = volume()
    x0, _, x1, _ = o.bounds
    return np.arange(x0 + fm.B / 2, x1 - fm.B / 2 + 1e-9, fm.B + fm.GAP)


# ---- 1. cover linework (black, transparent) ------------------------------------
def cover():
    o, z = volume()
    # a: one forefoot slice, walls + dimension lines
    s, ya, yb = fm.slice_solid(o, z, 34.5)
    fig, ax = blank(5, 3)
    fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    outline(ax, s, INK, 0.35)
    zt = s.bounds[3]
    for yy, txt in ((ya, f"{yb - ya:.1f}"),):
        ax.annotate("", (ya, -3), (yb, -3), arrowprops=dict(arrowstyle="<->", lw=0.3, color=INK))
        ax.text((ya + yb) / 2, -6, txt, ha="center", fontsize=5, color=INK)
    ax.annotate("", (yb + 3, 0), (yb + 3, zt), arrowprops=dict(arrowstyle="<->", lw=0.3, color=INK))
    ax.text(yb + 5, zt / 2, f"{zt:.1f}", rotation=90, va="center", fontsize=5, color=INK)
    ax.text(ya, zt + 4, "SLICE 05 / x 34.5", fontsize=5, color=INK)
    save(fig, "cover_slice.png", transparent=True)

    # b: oblique stack of slice outlines (axonometric feel)
    fig, ax = blank(7, 4)
    fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    k = np.radians(30)
    for x in slice_xs()[::3]:
        s, _, _ = fm.slice_solid(o, z, x)
        for g in getattr(s, "geoms", [s]):
            if not hasattr(g, "exterior"):
                continue
            for ring, lw in [(g.exterior, 0.3)] + [(r, 0.12) for r in g.interiors]:
                yy, zz = np.asarray(ring.xy)
                ax.plot(x * np.cos(k) + yy * 0.55, x * np.sin(k) * 0.45 + zz, color=INK, lw=lw)
    save(fig, "cover_stack.png", transparent=True)

    # c: side elevation, midsole section + posed foot at mid-stance + hinges
    side(row=30, dark=False, name="cover_side.png", transparent=True)


def side(row, dark, name, transparent=False, w=8, h=3.2):
    mo = pd.read_csv(MS / "foot_motion.csv")
    sh = fa.sections()
    fig, ax = blank(w, h, dark)
    if transparent:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    col = "white" if dark else INK
    for e in sh["sole"]:
        ax.plot(e[:, 0], e[:, 1], color=col, lw=0.2)
    rear, toe, shank = fa.pose(mo, row)
    for nm, a in (("toe", toe), ("mid", rear), ("heel", rear)):
        for e in sh[nm]:
            p = fa.rot(e, fa.MTP, a)
            ax.plot(p[:, 0], p[:, 1], color=col, lw=0.5)
    an = fa.rot(fa.ANKLE, fa.MTP, rear)
    d = np.array([-np.cos(np.radians(shank)), np.sin(np.radians(shank))])
    ax.plot([an[0], an[0] + 120 * d[0]], [an[1], an[1] + 120 * d[1]], color=col, lw=0.5)
    for pnt in (an, fa.rot(fa.HEEL_MID, fa.MTP, rear), fa.MTP):
        ax.plot(*pnt, "o", ms=2.5, color=ACCENT)
    ax.plot([-10, 285], [0, 0], color=col, lw=0.6)
    for xx in range(0, 271, 30):
        ax.plot([xx, xx], [-2, 0], color=col, lw=0.3)
    save(fig, name, dark=dark, transparent=transparent)


# ---- 2. stance filmstrip (white on black) ----------------------------------------
def film():
    mo = pd.read_csv(MS / "foot_motion.csv")
    sh = fa.sections()
    rows = np.linspace(0, len(mo) - 1, 12).astype(int)
    fig, axs = plt.subplots(3, 4, figsize=(10, 4.6))
    fig.patch.set_facecolor(DARK)
    for ax, k in zip(axs.flat, rows):
        ax.set_facecolor(DARK); ax.set_aspect("equal"); ax.axis("off")
        for e in sh["sole"]:
            ax.plot(e[:, 0], e[:, 1], color="#8a8780", lw=0.15)
        rear, toe, shank = fa.pose(mo, k)
        for nm, a in (("toe", toe), ("mid", rear), ("heel", rear)):
            for e in sh[nm]:
                p = fa.rot(e, fa.MTP, a)
                ax.fill(p[:, 0], p[:, 1], color="white", alpha=0.08, lw=0)
                ax.plot(p[:, 0], p[:, 1], color="white", lw=0.5)
        an = fa.rot(fa.ANKLE, fa.MTP, rear)
        d = np.array([-np.cos(np.radians(shank)), np.sin(np.radians(shank))])
        ax.plot([an[0], an[0] + 110 * d[0]], [an[1], an[1] + 110 * d[1]], color="white", lw=0.6)
        ax.set_xlim(-15, 285); ax.set_ylim(-5, 230)
        ax.text(0, 215, f"{mo.pct[k]:.0f} %", color="white", fontsize=7)
        ax.text(0, 195, f"{mo.F[k]:.0f} N", color="#8a8780", fontsize=6)
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.99, 0.04, 0.04)
    save(fig, "film_dark.png", dark=True)
    side(30, True, "side_dark.png", w=9, h=3.4)


# ---- 3. lattice contact sheet (white on black) -----------------------------------
def cell_walls(th, l, h, t, nx=2, ny=2):
    from shapely.geometry import LineString, box
    from shapely.ops import unary_union
    c, s = np.cos(np.radians(th)), np.sin(np.radians(th))
    px, py = 2 * l * c, 2 * (h + l * s)
    cell = [((0, 0), (0, h)), ((0, h), (l * c, h + l * s)), ((l * c, h + l * s), (px, h)),
            ((l * c, h + l * s), (l * c, 2 * h + l * s)), ((l * c, 2 * h + l * s), (0, py)),
            ((l * c, 2 * h + l * s), (px, py))]
    lines = [LineString([(a[0] + i * px, a[1] + j * py), (b[0] + i * px, b[1] + j * py)])
             for i in range(nx) for j in range(ny) for a, b in cell]
    w = unary_union(lines).buffer(t / 2, quad_segs=4)
    return w.intersection(box(-t / 2, 0, nx * px + t / 2, ny * py - abs(l * s) - t))


def grid():
    import zone_search as zs
    cand = zs.candidates()
    pick = pd.read_csv(MS / "zone_pick.csv", index_col=0)
    sel = cand.sort_values("E2").iloc[np.linspace(0, len(cand) - 1, 38).astype(int)]
    for _, pr in pick.iterrows():
        hit = cand[(cand.theta == pr.theta) & (cand.l == pr.l) & (cand.h == pr.h) & (abs(cand.t - pr.t) < 1e-6)]
        sel = pd.concat([sel, hit])
    sel = sel[~sel.index.duplicated()].sort_values("E2")
    fig, axs = plt.subplots(5, 8, figsize=(10, 6.6))
    fig.patch.set_facecolor(DARK)
    for ax, (_, r) in zip(axs.flat, sel.iterrows()):
        ax.set_facecolor(DARK); ax.set_aspect("equal"); ax.axis("off")
        chosen = any((abs(pick.theta - r.theta) < 1e-6) & (abs(pick.l - r.l) < 1e-6) &
                     (abs(pick.h - r.h) < 1e-6) & (abs(pick.t - r.t) < 1e-6))
        outline(ax, cell_walls(r.theta, r.l, r.h, r.t), ACCENT if chosen else "white",
                0.6 if chosen else 0.35)
        ax.set_title(f"E₂ {r.E2:.2f}   ν {r.nu21:.2f}", fontsize=5.5,
                     color=ACCENT if chosen else "#b9b6ae", pad=2)
    for ax in axs.flat:
        ax.set_facecolor(DARK); ax.axis("off")
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.97, 0.12, 0.3)
    save(fig, "grid_dark.png", dark=True)


# ---- 4. strain map (zoned lattice, whole stance) ---------------------------------
def strain():
    from fill_midsole import read_volume
    from foot_conform import ES, foot_bodies
    from full_stance import F_MIN, sink_solve
    from lattice_law import x2_law
    import zone_search as zs
    f = MS / "full_stance_zoned.csv"
    if not f.exists():
        pick = pd.read_csv(MS / "zone_pick.csv", index_col=0)
        laws = {z: x2_law(p.theta, p.h, p.l, p.t, ES) for z, p in pick.iterrows()}
        mo = pd.read_csv(MS / "foot_motion.csv")
        feet, (o, zt) = foot_bodies(), read_volume()
        rows = []
        for k in range(len(mo)):
            r = mo.loc[k]
            if r.F < F_MIN:
                continue
            sink, tops = sink_solve(feet, o, zt, lambda x: laws[zs.zone_of(x)], r)
            rows += [{"pct": r.pct, "x": t["x"], "eps": t["eps"].max(), "force": t["force"],
                      "sink": sink} for t in tops]
        pd.DataFrame(rows).to_csv(f, index=False)
    d = pd.read_csv(f)
    m = d.pivot(index="x", columns="pct", values="eps") * 100
    fig, ax = plt.subplots(figsize=(8, 3.6))
    fig.patch.set_facecolor(PAPER); ax.set_facecolor(PAPER)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", [PAPER, "#f1d9cf", "#d98a73", "#a8402a", "#5e1a0f"])
    im = ax.pcolormesh(m.columns, m.index, m.values, cmap=cmap, shading="nearest", vmin=0, vmax=22)
    for xz, lab in ((63.1, "MTP hinge"), (186.6, "heel–midfoot hinge")):
        ax.axhline(xz, color=INK, lw=0.4, ls=(0, (2, 2)))
        ax.text(101, xz + 3, lab, fontsize=6.5, color=MUTED, ha="right")
    for z, y in (("forefoot", 30), ("midfoot", 125), ("heel", 228)):
        ax.text(1.5, y, z, fontsize=7, ha="left", va="center", style="italic", color=MUTED)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.outline.set_linewidth(0.3); cb.ax.tick_params(labelsize=6.5, width=0.3)
    cb.set_label("peak column strain (%)", fontsize=7)
    ax.set_xlabel("stance (%)", fontsize=7.5); ax.set_ylabel("slice x (mm)", fontsize=7.5)
    ax.tick_params(labelsize=6.5, width=0.3, color=MUTED)
    for sp in ax.spines.values():
        sp.set_linewidth(0.3); sp.set_color(MUTED)
    save(fig, "strain_map.png")
    s = d.groupby("pct").sink.first()
    fig, ax = plt.subplots(figsize=(4, 1.6))
    ax.plot(s.index, s.values, color=INK, lw=1)
    ax.fill_between(s.index, s.values, color=WARM)
    ax.set_ylabel("sink (mm)", fontsize=7); ax.set_xlabel("stance (%)", fontsize=7)
    ax.tick_params(labelsize=6); [sp.set_linewidth(0.3) for sp in ax.spines.values()]
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "sink.png")


# ---- 5. stance data chart ---------------------------------------------------------
def stance():
    mo = pd.read_csv(MS / "foot_motion.csv")
    fig, axs = plt.subplots(2, 1, figsize=(5, 3.6), sharex=True)
    for ax in axs:
        ax.set_facecolor(DARK); ax.tick_params(colors="#b9b6ae", labelsize=6.5, width=0.3)
        for sp in ax.spines.values():
            sp.set_color("#55534e"); sp.set_linewidth(0.4)
        ax.spines[["top", "right"]].set_visible(False)
    fig.patch.set_facecolor(DARK)
    axs[0].fill_between(mo.pct, mo.F_fore, color="white", alpha=0.85, lw=0)
    axs[0].fill_between(mo.pct, mo.F_fore, mo.F, color=ACCENT, lw=0)
    axs[0].text(40, 600, "forefoot", color=DARK, fontsize=7)
    i = int(np.argmin(abs(mo.pct - 22))); axs[0].annotate("heel", (mo.pct[i], (mo.F_fore[i] + mo.F[i]) / 2), xytext=(8, 1500), color=ACCENT, fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.6, color=ACCENT))
    axs[0].set_ylim(0, 1700)
    axs[0].set_ylabel("force (N)", color="#b9b6ae", fontsize=7)
    for c, ls in (("pitch", "-"), ("ankle", (0, (3, 2))), ("mtp", (0, (1, 1.5)))):
        axs[1].plot(mo.pct, mo[c], color="white", lw=0.9, ls=ls)
        axs[1].text(101, mo[c].iloc[-1], c, color="white", fontsize=6.5, va="center")
    axs[1].set_ylabel("angle (°)", color="#b9b6ae", fontsize=7)
    axs[1].set_xlabel("stance (%)", color="#b9b6ae", fontsize=7)
    save(fig, "stance_dark.png", dark=True)


# ---- 6. pyvista renders ----------------------------------------------------------
def renders():
    import pyvista as pv
    import trimesh
    pv.OFF_SCREEN = True
    sole = pv.read(MS / "midsole_filled.stl")
    # a: section cut, dark
    cut = sole.clip(normal="x", origin=(38.0, 0, 0), invert=False)
    p = pv.Plotter(off_screen=True, window_size=(3000, 1100))
    p.set_background(DARK)
    p.add_mesh(cut, color=CLAY, smooth_shading=False, ambient=0.25, diffuse=0.75, specular=0.05)
    p.enable_parallel_projection()
    p.camera_position = [(-420, -260, 170), (150, 50, 20), (0, 0, 1)]
    p.camera.parallel_scale = 78
    p.enable_eye_dome_lighting()
    p.screenshot(FIG / "section_render.png")
    p.close()
    trim(FIG / "section_render.png")
    print("  wrote section_render.png")
    # b: exploded, warm grey: bottom plate / slices / top plate / foot
    o, z = volume()
    parts = {"bottom": [trimesh.creation.extrude_polygon(o, fm.PLATE_BOT)], "slices": [], "top": []}
    for x in slice_xs():
        s, _, _ = fm.slice_solid(o, z, x)
        parts["slices"].append(fm.extrude_yz(s, x - fm.B / 2, x + fm.B / 2))
    x0, _, x1, _ = o.bounds
    for x in np.arange(x0 + 1, x1, 2.0):
        b = fm.top_band(o, z, x)
        if b.area > 1e-6:
            parts["top"].append(fm.extrude_yz(b, x - 1, x + 1))
    import foot_conform as fc
    mo = pd.read_csv(MS / "foot_motion.csv")
    foot = fc.pose(fc.foot_bodies(), mo.loc[30])
    lift = {"bottom": -35, "slices": 0, "top": 45, "foot": 95}
    p = pv.Plotter(off_screen=True, window_size=(3000, 2000))
    p.set_background(WARM)
    for k, ms in parts.items():
        m = pv.wrap(trimesh.util.concatenate(ms))
        m.translate((0, 0, lift[k]), inplace=True)
        p.add_mesh(m, color=CLAY if k == "slices" else "#efe9df", ambient=0.3, diffuse=0.7)
    for m in foot.values():
        mm = pv.wrap(m); mm.translate((0, 0, lift["foot"]), inplace=True)
        p.add_mesh(mm, color=ACCENT, opacity=0.35)
        ed = mm.extract_feature_edges(30)
        if ed.n_points:
            p.add_mesh(ed, color=ACCENT, line_width=2)
    p.enable_parallel_projection()
    p.camera_position = [(-380, -420, 440), (135, 55, 60), (0, 0, 1)]
    p.camera.parallel_scale = 190
    p.enable_eye_dome_lighting()
    p.screenshot(FIG / "exploded_render.png")
    p.close()
    trim(FIG / "exploded_render.png")
    print("  wrote exploded_render.png")


# ---- 7. load-path cut-away (pastel) ----------------------------------------------
def loadpath():
    o, z = volume()
    s, ya, yb = fm.slice_solid(o, z, 34.5)
    zt = s.bounds[3]
    fig, ax = plt.subplots(figsize=(7, 5.6))
    fig.patch.set_facecolor(PAPER); ax.set_facecolor(PAPER)
    ax.set_aspect("equal"); ax.axis("off")
    from shapely.geometry import box
    ax.fill([-30, 140, 140, -30], [-25, -25, 0, 0], color=PASTEL["ground"], lw=0)
    ax.plot([-30, 140], [0, 0], color=INK, lw=1)
    lat = s.intersection(box(ya, fm.PLATE_BOT, yb, zt - fm.PLATE_TOP))
    outline(ax, s, INK, 0.6, fill=PASTEL["lattice"], hole="white")
    ax.fill([ya, yb, yb, ya], [0, 0, fm.PLATE_BOT, fm.PLATE_BOT], color=PASTEL["plate"], lw=0.8, ec=INK)
    yy = np.linspace(ya, yb, 60)
    ztop = z(np.full_like(yy, 34.5), yy)
    ax.fill_between(yy, ztop - fm.PLATE_TOP, ztop, color=PASTEL["plate"], lw=0.8, ec=INK)
    foot_y = np.linspace(ya + 6, yb - 6, 50)
    ax.fill_between(foot_y, ztop[5:55].max() + 0.5 + 5 * (1 - np.cos(np.linspace(0, 2 * np.pi, 50))) / 2,
                    zt + 45, color=PASTEL["foot"], lw=1, ec=INK)
    for yv in np.linspace(ya + 15, yb - 15, 5):
        ax.annotate("", (yv, zt + 3), (yv, zt + 22),
                    arrowprops=dict(arrowstyle="-|>", lw=1.2, color=INK))
    for yv in np.linspace(ya + 10, yb - 10, 7):
        ax.annotate("", (yv, -14), (yv, -3), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK))
    tags = [(yb - 8, zt + 32, 1, "pseudo-foot"), (yb + 6, zt - 1, 2, "top plate"),
            (yb + 6, zt * 0.55, 3, "re-entrant slice"), (yb + 6, 1, 4, "bottom plate"),
            (yb + 6, -12, 5, "ground reaction")]
    for x, y, n, lab in tags:
        ax.add_patch(plt.Circle((x, y), 3.2, color="white", ec=INK, lw=0.8, zorder=5))
        ax.text(x, y, str(n), ha="center", va="center", fontsize=7, zorder=6)
        ax.text(x + 5, y, lab, va="center", fontsize=8, style="italic")
    ax.set_xlim(ya - 12, yb + 45); ax.set_ylim(-25, zt + 52)
    save(fig, "loadpath.png")


def models():
    import os, subprocess
    subprocess.run([sys.executable, str(_ROOT / "models" / "model_compare.py")], cwd=MS, check=True,
                   env={**os.environ, "REPORT": "1", "MPLBACKEND": "Agg"}, capture_output=True)
    for n in ("model_compare_moduli.png", "model_compare_curves.png"):
        shutil.copy(MS / n, FIG / n)
    print("  copied model figures")


ALL = {"cover": cover, "film": film, "grid": grid, "strain": strain, "stance": stance,
       "renders": renders, "loadpath": loadpath, "models": models}

if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    serif()
    for name in (sys.argv[1:] or ALL):
        print(name)
        ALL[name]()
