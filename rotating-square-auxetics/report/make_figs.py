"""Figures for the rotating-squares report (same style as polished_files/report_v2).

Run from this folder:  python make_figs.py            (all figures)
                       python make_figs.py film sweep (only some, by name)
Inputs : ../results/data/*.csv (MATLAB results) and data/fe_*.csv (export_fe_for_report.m)
Output : fig/*.png
The closed-form models are re-implemented here and checked against the MATLAB CSV first.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "results" / "data"
FIG = HERE / "fig"
FIG.mkdir(exist_ok=True)

INK, MUTED, PAPER, ACCENT, GRID, DARK, BEIGE = ("#1b1b1a", "#6b6a66", "#fcfcfb", "#c0392b",
                                               "#d9d8d3", "#141413", "#f4efe8")
SQ = "#9fb9d6"          # square fill (neutral blue-grey; red stays reserved for data)
GOOD, WARN = "#5a9e6f", "#d9a441"

# ---- design (same values as params_default.m) ------------------------------------
a, Lh, w, t = 10.0, 1.0, 0.8, 5.0
th0 = np.radians(15.0)
E, sy, SF, Kt = 3500.0, 50.0, 1.5, 1.0
N = 5


def style():
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 8.5,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.5, "axes.spines.top": False,
        "axes.spines.right": False, "axes.facecolor": PAPER, "figure.facecolor": PAPER,
        "savefig.facecolor": PAPER, "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.labelcolor": INK,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.4,
        "grid.linestyle": (0, (1, 3)), "legend.frameon": False})


def save(fig, name, **kw):
    fig.savefig(FIG / name, dpi=220, bbox_inches="tight", pad_inches=0.04, **kw)
    plt.close(fig)
    print("  wrote", name)


# ---- models (mirror the MATLAB functions) --------------------------------------------
def pitch(th, Lh_=Lh, a_=a):
    return a_ * (np.cos(th) + np.sin(th)) + Lh_ * np.cos(th - th0)


def dpitch(th, Lh_=Lh, a_=a):
    return a_ * (np.cos(th) - np.sin(th)) - Lh_ * np.sin(th - th0)


def props(a_=a, Lh_=Lh, w_=w):
    k_r = E * t * w_ ** 3 / 12 / Lh_
    b = a_ * (np.cos(th0) - np.sin(th0))
    rho = w_ ** 2 / (12 * b ** 2)
    return dict(k_r=k_r, b=b, rho=rho, nu=-1 / (1 + 8 * rho),
                E=8 * k_r / (t * b ** 2 * (1 + 8 * rho)), p0=a_ * (np.cos(th0) + np.sin(th0)) + Lh_)


def nl_path(th):
    C = props()
    p0 = pitch(th0)
    eps = pitch(th) / p0 - 1
    sig = 8 * C["k_r"] * (th - th0) / (p0 * t * dpitch(th))
    sh = Kt * (E * w * np.abs(th - th0) / Lh + np.abs(sig) * p0 / w)
    return eps, sig, sh


def theta_lock():
    from scipy.optimize import brentq
    return brentq(dpitch, th0, np.pi / 2)


ROT = lambda al: np.array([[np.cos(al), -np.sin(al)], [np.sin(al), np.cos(al)]])
CORNER = (a / 2) * np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]])
HINGE = {"right": (0, 3), "left": (2, 1), "top": (1, 0), "bottom": (3, 2)}   # 0-based, (A, B)


def layout(th, n=N, npts=14):
    """Squares (corners) and Hermite ligaments at rotation th, as in rs_layout.m."""
    p = pitch(th)
    sq, typ, al = {}, {}, {}
    for j in range(n):
        for i in range(n):
            ty = (i + j) % 2                       # 0 = A (-th), 1 = B (+th)
            typ[i, j], al[i, j] = ty, (-1 if ty == 0 else 1) * th
            sq[i, j] = np.array([i * p, j * p]) + CORNER @ ROT(al[i, j]).T
    s = np.linspace(0, 1, npts)[:, None]
    h00, h10, h01, h11 = 2*s**3 - 3*s**2 + 1, s**3 - 2*s**2 + s, -2*s**3 + 3*s**2, s**3 - s**2
    ligs = []
    for j in range(n):
        for i in range(n):
            for (i2, j2, s1, s2, e) in ((i + 1, j, "right", "left", (1, 0)), (i, j + 1, "top", "bottom", (0, 1))):
                if i2 >= n or j2 >= n:
                    continue
                P0 = sq[i, j][HINGE[s1][typ[i, j]]]
                P1 = sq[i2, j2][HINGE[s2][typ[i2, j2]]]
                dA = al[i, j] - (-1 if typ[i, j] == 0 else 1) * th0
                dB = al[i2, j2] - (-1 if typ[i2, j2] == 0 else 1) * th0
                t0 = Lh * ROT(dA) @ np.array(e)
                t1 = Lh * ROT(dB) @ np.array(e)
                ligs.append(h00 * P0 + h10 * t0 + h01 * P1 + h11 * t1)
    c = np.mean([np.array([i, j]) * p for i in range(n) for j in range(n)], axis=0)
    return [v - c for v in sq.values()], [l - c for l in ligs]


def draw(ax, th, fill=SQ, edge=INK, lig=INK, lw=0.6, ghost=None, n=N):
    from matplotlib.patches import Polygon
    if ghost is not None:
        for v in layout(th0, n)[0]:
            ax.add_patch(Polygon(v, closed=True, fill=False, ec=ghost, lw=0.5, ls=(0, (3, 2))))
    sq, ligs = layout(th, n)
    for v in sq:
        ax.add_patch(Polygon(v, closed=True, fc=fill, ec=edge, lw=lw))
    for l in ligs:
        ax.plot(l[:, 0], l[:, 1], color=lig, lw=lw * 2.2, solid_capstyle="round")
    ax.set_aspect("equal"); ax.axis("off")


# ---- 0. consistency with MATLAB ------------------------------------------------------
def check():
    M = pd.read_csv(DATA / "S6_matlab_results.csv").dropna(subset=["theta_deg"])
    th = np.radians(M.theta_deg.values)
    eps, sig, sh = nl_path(th)
    e1 = np.max(np.abs(eps - M.strain_x.values))
    e2 = np.max(np.abs(sig - M.sigma_lattice_MPa.values))
    e3 = np.max(np.abs(sh * 1.06496 - M.hinge_stress_MPa.values) / np.maximum(M.hinge_stress_MPa.values, 1e-9))
    C = props()
    print(f"  check vs MATLAB: strain {e1:.1e}, stress {e2:.1e} MPa, nu {C['nu']:.5f}, E {C['E']:.3f} MPa,"
          f" hinge (rel, edge factor approx) {e3:.1e}")
    # theta in the CSV comes from interp1 on the strain path, hence ~1e-7 in strain
    assert e1 < 1e-5 and e2 < 1e-6 and abs(C["nu"] + 0.991539) < 1e-5


# ---- figures -------------------------------------------------------------------------
def cover():
    """Faint technical line drawings for the cover (transparent background)."""
    for name, th, n in (("cover_open.png", np.radians(32), 6), ("cover_closed.png", th0, 6)):
        fig, ax = plt.subplots(figsize=(6, 6))
        draw(ax, th, fill="none", edge=INK, lig=INK, lw=0.45, n=n)
        sq, _ = layout(th, n)
        # a few dimension marks, like a drawing
        v = sq[0]
        ax.annotate("", xy=v[1], xytext=v[2], arrowprops=dict(arrowstyle="<->", lw=0.4, color=MUTED))
        ax.autoscale()
        save(fig, name, transparent=True)


def timeline():
    ev = [(1987, "Lakes — re-entrant\nauxetic foam"),
          (1991, "Evans et al. — the\nword 'auxetic'"),
          (2000, "Grima & Evans —\nrotating squares"),
          (2005, "Grima et al. —\nrotating rigid units"),
          (2010, "Bertoldi et al. — holes\nbuckle into squares"),
          (2015, "Shan et al. — planar\ncut patterns"),
          (2016, "Rafsanjani & Pasini —\nbistable rotating units"),
          (2018, "Coulais et al. — size\neffects, length scale"),
          (2022, "Czajkowski et al. —\nconformal elasticity")]
    fig, ax = plt.subplots(figsize=(7.4, 1.95))
    ax.axis("off")
    xs = np.arange(len(ev), dtype=float)
    ax.plot([-0.4, len(ev) - 0.6], [0, 0], color=INK, lw=0.8)
    for i, (x, (y, txt)) in enumerate(zip(xs, ev)):
        up = i % 2 == 0
        ax.plot([x, x], [0, 0.35 if up else -0.35], color=MUTED, lw=0.5)
        ax.plot(x, 0, "o", ms=4, color=ACCENT if y in (2000, 2022) else INK)
        ax.text(x, 0.42 if up else -0.42, f"{y}\n{txt}" if up else f"{txt}\n{y}", ha="center",
                va="bottom" if up else "top", fontsize=6.1, color=INK, linespacing=1.25)
    ax.set_xlim(-0.7, len(ev) - 0.3); ax.set_ylim(-1.3, 1.3)
    save(fig, "timeline.png")


def geometry():
    fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.15))
    lim = 0.5 * 3 * a * np.sqrt(2) + a * 0.9
    for ax, d in zip(axs, (0, 15, 30, 45)):
        th = np.radians(d)
        sq, _ = layout(th, 4)
        from matplotlib.patches import Polygon
        for v in sq:
            ax.add_patch(Polygon(v, closed=True, fc=SQ if d != 15 else "#c9d6e6", ec=INK, lw=0.5))
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.axis("off")
        p = a * (np.cos(th) + np.sin(th))
        ax.set_title(f"θ = {d}°" + ("  (as printed)" if d == 15 else ""), fontsize=7.5, color=INK, pad=2)
        ax.text(0, -lim - 0.5, f"pitch {p:.2f} mm\nporosity {100*np.sin(2*th)/(1+np.sin(2*th)):.0f} %",
                ha="center", va="top", fontsize=6.2, color=MUTED, linespacing=1.3)
    save(fig, "geometry.png")


def loop_diagram():
    """Beige page: vector loop closure through one hinge (ideal pins, Lh = 0)."""
    from matplotlib.patches import Polygon, Arc, FancyArrowPatch
    th = np.radians(20)
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    fig.patch.set_facecolor(BEIGE); ax.set_facecolor(BEIGE); ax.axis("off")
    p = a * (np.cos(th) + np.sin(th))
    cA, cB, cC = np.array([0, 0]), np.array([p, 0]), np.array([0, p])
    for c, s, lab in ((cA, -1, "A  (−θ)"), (cB, 1, "B  (+θ)"), (cC, 1, "B  (+θ)")):
        v = c + CORNER @ ROT(s * th).T
        ax.add_patch(Polygon(v, closed=True, fc="#e6dccd", ec=INK, lw=0.7))
        ax.add_patch(Polygon(c + CORNER, closed=True, fill=False, ec=MUTED, lw=0.4, ls=(0, (2, 2))))
        ax.plot(*c, "o", ms=3, color=INK)
        ax.text(c[0], c[1] + 1.0, lab, ha="center", va="bottom", fontsize=7.5, color=INK)
    hA = cA + ROT(-th) @ CORNER[0]
    arr = dict(arrowstyle="-|>", lw=0.9, mutation_scale=8)
    ax.add_patch(FancyArrowPatch(cA, hA, color=ACCENT, **arr))
    ax.add_patch(FancyArrowPatch(cB, hA, color=ACCENT, **arr))
    ax.text(cA[0] - 1.0, cA[1] - 3.4, r"$\mathbf{r}_A=R(-\theta)\,\frac{a}{2}(1,1)$", fontsize=7.2,
            color=ACCENT, ha="center")
    ax.text(cB[0] + 1.0, cB[1] - 3.4, r"$\mathbf{r}_B=R(\theta)\,\frac{a}{2}(-1,1)$", fontsize=7.2,
            color=ACCENT, ha="center")
    ax.plot(*hA, "o", ms=5, mfc=PAPER, mec=ACCENT, mew=1)
    ax.text(hA[0] - 0.3, hA[1] + 1.4, "hinge", fontsize=7, color=ACCENT, ha="right")
    ax.annotate("", xy=(p, -a * 1.05), xytext=(0, -a * 1.05),
                arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK))
    ax.text(p / 2, -a * 1.12, r"$p=a(\cos\theta+\sin\theta)$", ha="center", va="top", fontsize=8.5, color=INK)
    ax.annotate("", xy=(-a * 1.0, p), xytext=(-a * 1.0, 0), arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK))
    ax.text(-a * 1.07, p / 2, r"$p$", ha="right", va="center", fontsize=8.5, color=INK)
    ax.add_patch(Arc(cB, 6, 6, theta1=90, theta2=90 + np.degrees(th), color=INK, lw=0.6))
    ax.text(cB[0] - 1.6, cB[1] + 3.4, "θ", fontsize=8, color=INK)
    ax.text(p + 7.5, p * 0.95,
            "Loop closure through one hinge:\n"
            r"$\mathbf{c}_B=\mathbf{c}_A+\mathbf{r}_A-\mathbf{r}_B$" "\n\n"
            "x:  the pitch grows as $a(\\cos\\theta+\\sin\\theta)$\n"
            "y:  the terms cancel — B stays level\n\n"
            "Same result upwards (A to the square above),\n"
            "so $\\varepsilon_x=\\varepsilon_y$ and $\\nu=-1$.\n\n"
            "Dashed: the closed state, θ = 0.",
            fontsize=7.3, color=INK, va="top", linespacing=1.45)
    ax.set_xlim(-a * 1.35, p + 25); ax.set_ylim(-a * 1.35, p + a * 0.9); ax.set_aspect("equal")
    save(fig, "loop.png", facecolor=BEIGE)


def kinematics():
    th = np.radians(np.linspace(0, 60, 400))
    eps_ideal = (np.cos(th) + np.sin(th)) / (np.cos(th0) + np.sin(th0)) - 1
    eps_arm = pitch(th) / pitch(th0) - 1
    tl = theta_lock()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 2.5), gridspec_kw=dict(wspace=0.32))
    a1.plot(np.degrees(th), 100 * eps_ideal, color=INK, lw=1.1)
    a1.plot(np.degrees(th), 100 * eps_arm, color=MUTED, lw=1.1, ls=(0, (4, 2)))
    a1.text(51, 17.0, "ideal pins", fontsize=6.8, color=INK)
    a1.text(46.5, 5.5, f"with L$_h$ = {Lh:g} mm", fontsize=6.8, color=MUTED)
    for x in (45, np.degrees(tl)):
        a1.axvline(x, color=MUTED, lw=0.5)
    a1.text(45.8, -18.5, f"open locks:\n45° (pins)\n{np.degrees(tl):.1f}° (L$_h$)", fontsize=6.2,
            color=MUTED, va="bottom")
    a1.set_ylim(-20, 19.5)
    a1.plot(15, 0, "o", ms=3.5, color=INK); a1.text(16.5, -2.8, "as printed", fontsize=6.6, color=INK)
    a1.set_xlabel("square rotation θ (°)"); a1.set_ylabel("strain ε$_x$ = ε$_y$ (%)")
    a1.set_xlim(0, 60)
    a2.plot(np.degrees(th), 100 * np.sin(2 * th) / (1 + np.sin(2 * th)), color=INK, lw=1.1)
    a2.axvline(45, color=MUTED, lw=0.5)
    a2.set_xlabel("square rotation θ (°)"); a2.set_ylabel("porosity (%)"); a2.set_xlim(0, 60)
    a2.text(46, 8, "50 % at 45°", fontsize=6.6, color=MUTED)
    save(fig, "kinematics.png")


def nu_width():
    W = pd.read_csv(DATA / "S3_nu_vs_width.csv")
    ww = np.linspace(0.2, 3.0, 200)
    fig, ax = plt.subplots(figsize=(4.6, 2.75))
    for aa, ls in ((5, (0, (4, 2))), (10, "-"), (20, (0, (1, 1.5)))):
        nu = [props(a_=aa, w_=x)["nu"] for x in ww]
        ax.plot(ww, nu, color=INK, lw=1.0, ls=ls)
        ax.text(3.12, nu[-1], f"a = {aa} mm", fontsize=6.6, va="center", color=INK)
    ax.plot(W.w_mm, W.nu_FE_periodic, "D", ms=4.2, mfc="none", mec=INK, mew=0.8)
    ax.plot(W.w_mm, W.nu_FE_5x5, "o", ms=3.2, color=MUTED)
    ax.text(1.45, -0.925, "◇ own FE, periodic cell\n● own FE, 5 × 5 array", fontsize=6.5, color=MUTED)
    ax.axvline(0.8, color=MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.text(0.83, -0.69, "min. printable w", fontsize=6.4, color=MUTED)
    ax.axhline(-1, color=MUTED, lw=0.5)
    ax.set_xlabel("ligament width w (mm)"); ax.set_ylabel("Poisson's ratio ν")
    ax.set_xlim(0.2, 3.08); ax.set_ylim(-1.02, -0.65)
    save(fig, "nu_width.png")


def sweeps():
    base = props()
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.25), gridspec_kw=dict(wspace=0.35))
    for ax, (key, x, lab, dflt) in zip(axs, (
            ("w_", np.linspace(0.2, 3, 120), "ligament width w (mm)", w),
            ("Lh_", np.linspace(0.4, 4, 120), "ligament length L$_h$ (mm)", Lh),
            ("a_", np.linspace(5, 20, 120), "square side a (mm)", a))):
        y = [props(**{key: v})["E"] for v in x]
        ax.semilogy(x, y, color=INK, lw=1.0)
        ax.plot(dflt, base["E"], "o", ms=4, color=INK)
        ax.set_xlabel(lab)
    axs[0].set_ylabel("E$_{eff}$ (MPa)")
    axs[0].annotate(f"design: {base['E']:.1f} MPa", (w, base["E"]), (1.3, 6), fontsize=6.6, color=MUTED,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED))
    save(fig, "sweeps.png")


def stress():
    tl = theta_lock()
    th = np.linspace(0, tl - 0.05 * (tl - th0) * 0.3, 600)
    eps, sig, sh = nl_path(th)
    C = props()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 2.6), gridspec_kw=dict(wspace=0.3))
    a1.plot(100 * eps, sig, color=INK, lw=1.1)
    a1.plot(100 * eps, C["E"] * eps, color=MUTED, lw=0.9, ls=(0, (4, 2)))
    a1.text(-7.5, -3.9, f"linear, E$_{{eff}}$ = {C['E']:.1f} MPa", fontsize=6.6, color=MUTED)
    a1.text(-15, -0.9, "large rotation", fontsize=6.6, color=INK)
    a1.set_ylim(-4.5, 12); a1.set_xlabel("strain ε$_x$ (%)"); a1.set_ylabel("lattice stress σ$_x$ (MPa)")
    a1.axvline(100 * (pitch(tl) / pitch(th0) - 1), color=MUTED, lw=0.5)
    a1.text(100 * (pitch(tl) / pitch(th0) - 1) - 0.4, -4.2, "open\nlock", fontsize=6.3, color=MUTED,
            ha="right", va="bottom")
    m = np.abs(eps) < 0.03
    a2.plot(100 * eps[m], sh[m], color=INK, lw=1.1)
    a2.axhline(sy, color=MUTED, lw=0.6); a2.axhline(sy / SF, color=MUTED, lw=0.6, ls=(0, (4, 2)))
    a2.text(-2.9, sy + 2, "yield 50 MPa", fontsize=6.5, color=MUTED)
    a2.text(-2.9, sy / SF + 2, "yield / 1.5", fontsize=6.5, color=MUTED)
    L = pd.read_csv(DATA / "S4_limits.csv").set_index("Limit").strain_pct
    for v in (L["first yield (tension)"], L["first yield (compression)"]):
        a2.plot(v, sy, "o", ms=3.5, color=INK)
        a2.text(v + (0.15 if v > 0 else -0.15), sy - 9, f"{v:+.2f} %", fontsize=6.5, color=INK,
                ha="left" if v > 0 else "right")
    a2.set_xlim(-3, 3); a2.set_ylim(0, 110)
    a2.set_xlabel("strain ε$_x$ (%)"); a2.set_ylabel("peak ligament stress (MPa)")
    save(fig, "stress.png")


def limits():
    L = pd.read_csv(DATA / "S4_limits.csv").set_index("Limit").strain_pct
    rows = [("safe range, tension", L["SAFE max"]), ("safe range, compression", L["SAFE min"]),
            ("first yield, tension", L["first yield (tension)"]),
            ("first yield, compression", L["first yield (compression)"]),
            ("open lock", L["open lock"]), ("square contact", L["square contact"]),
            ("in-plane buckling (linear FE)", L["in-plane buckling"])]
    fig, ax = plt.subplots(figsize=(6.2, 2.3))
    for k, (lab, v) in enumerate(rows):
        y = len(rows) - k
        col = INK if k < 2 else MUTED
        ax.plot([0.1, abs(v)], [y, y], color=col, lw=1.2 if k < 2 else 0.8)
        ax.plot(abs(v), y, "o", ms=4, color=col)
        ax.text(abs(v) * 1.22, y, f"{v:+.2f} %", va="center", fontsize=6.6, color=col)
    ax.set_xscale("log"); ax.set_xlim(0.1, 2000)
    ax.set_yticks(range(len(rows), 0, -1)); ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_xlabel("|lattice strain| (%, log scale)")
    ax.grid(axis="y", visible=False)
    save(fig, "limits.png")


def fe():
    X = np.loadtxt(HERE / "data" / "fe_nodes.csv", delimiter=",")
    El = np.loadtxt(HERE / "data" / "fe_elems.csv", delimiter=",").astype(int)
    S = np.loadtxt(HERE / "data" / "fe_squares.csv", delimiter=",").astype(int)
    from matplotlib.patches import Polygon
    sc = 5.0
    c = X[:, :2].mean(axis=0)
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.3), gridspec_kw=dict(wspace=0.08))
    # left: the mesh (undeformed), right: deformed x5 over ghost
    for ax, deformed in ((axs[0], False), (axs[1], True)):
        Y = X[:, :2] + (sc * X[:, 2:4] if deformed else 0) - c
        for sq in S:
            idx = sq - 1
            if deformed:
                ax.add_patch(Polygon(X[idx, :2] - c, closed=True, fill=False, ec=MUTED, lw=0.4, ls=(0, (2, 2))))
            ax.add_patch(Polygon(Y[idx], closed=True, fc=SQ, ec=INK, lw=0.4))
        for n1, n2, lig in El:
            if lig or not deformed:
                ax.plot(Y[[n1 - 1, n2 - 1], 0], Y[[n1 - 1, n2 - 1], 1],
                        color=INK if lig else MUTED, lw=1.4 if lig else 0.35)
        if not deformed:
            ax.plot(Y[:, 0], Y[:, 1], ".", ms=1.2, color=INK)
        ax.set_aspect("equal"); ax.axis("off")
    axs[0].set_title("mesh: 4 beams per square (grey), 4 per ligament", fontsize=7, color=MUTED)
    axs[1].set_title("ε$_x$ = +1 %, displacements × 5, dashed: undeformed", fontsize=7, color=MUTED)
    save(fig, "fe.png")


def film():
    """Dark spread: 10 frames through one open/close cycle, ligaments by stress status."""
    tl = theta_lock()
    ths = np.concatenate([np.linspace(0.0, th0, 3, endpoint=False), np.linspace(th0, th0 + 0.97 * (tl - th0), 7)])
    L = pd.read_csv(DATA / "S4_limits.csv").set_index("Limit").strain_pct
    fig, axs = plt.subplots(2, 5, figsize=(7.4, 3.6), gridspec_kw=dict(wspace=0.05, hspace=0.18))
    fig.patch.set_facecolor(DARK)
    lim = 0.5 * (N - 1) * pitch(tl) + a
    for ax, th in zip(axs.flat, ths):
        ax.set_facecolor(DARK)
        eps, _, sh = nl_path(np.array([th]))
        s = sh[0] * 1.065
        col = GOOD if s <= sy / SF else (WARN if s <= sy else ACCENT)
        draw(ax, th, fill="#3a3a37", edge="#bdbab2", lig=col, lw=0.45, ghost="#5b5a55")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_title(f"θ {np.degrees(th):4.1f}°  ·  ε {100*eps[0]:+5.1f} %", fontsize=6.5, color="#d8d5ce", pad=2)
    save(fig, "film_dark.png", facecolor=DARK)


ALL = {"cover": cover, "timeline": timeline, "geometry": geometry, "loop": loop_diagram,
       "kin": kinematics, "nu": nu_width, "sweep": sweeps, "stress": stress, "limits": limits,
       "fe": fe, "film": film}

if __name__ == "__main__":
    style()
    check()
    for name in (sys.argv[1:] or ALL):
        ALL[name]()
