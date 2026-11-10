"""Compare the four lattice models with each other and (when available) ANSYS.

For each zone lattice in zone_pick.csv: E1, E2, nu12, nu21 from G&A, M&E, M&G
(formulas as the .m files) and Ohtaki's first step; plus the x2 compression
sigma-eps curve (Ohtaki large deflection vs the linear models).
If ansys/results/small_strain.csv / x2_compression.csv exist (forefoot
cell, see ansys/unit_cell_setup.md), they are added with % errors.

Run:  python model_compare.py
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lattice_law import oh_curve, x2_law

HERE = OUTPUTS
RES = ANSYS_RESULTS
ES, NUS = 19.3, 0.45


def main():
    zp = pd.read_csv(HERE / "zone_pick.csv", index_col=0)
    rows = []
    for zone, p in zp.iterrows():
        if zone == "heel" and (p[["theta", "l", "h", "t"]] == zp.loc["mid", ["theta", "l", "h", "t"]]).all():
            continue                                   # same cell as midfoot
        for name, f in MODELS.items():
            rows.append({"zone": zone, "model": name, **f(p.theta, p.h, p.l, p.t)})
    df = pd.DataFrame(rows)
    ans = ansys()
    if ans is not None:
        df = pd.concat([df, pd.DataFrame([{"zone": "fore", "model": "ANSYS", **ans}])])
    print(df.round(3).to_string(index=False))
    if ans is not None:
        errors(df)
    else:
        print("\n(no ANSYS results yet: put CSVs in ansys/results/)")
    df.to_csv(HERE / "model_compare.csv", index=False, float_format="%.4f")
    style()
    moduli_fig(df)
    curves(zp.loc["fore"])


# ---- the four models (same formulas as the .m files) -----------------------
def shorthand(th, h, l, t):
    return h / l, t / l, np.sin(np.radians(th)), np.cos(np.radians(th))


def compliance_form(th, h, l, t, Es, Cp, Ca, Cv, g_eff=None):
    a, g, s, c = shorthand(th, h, l, t)
    g = g if g_eff is None else g_eff
    D1 = s * s * Cp + c * c * Ca
    D2 = c * c * Cp + s * s * Ca + Cv
    return {"E1": Es * g**3 * c / ((a + s) * D1), "E2": Es * g**3 * (a + s) / (c * D2),
            "nu12": s * c * c * (Cp - Ca) / ((a + s) * D1), "nu21": s * (a + s) * (Cp - Ca) / D2}


def g_a(th, h, l, t):
    return compliance_form(th, h, l, t, ES, 1.0, 0.0, 0.0)


def m_e(th, h, l, t):
    a, g, _, _ = shorthand(th, h, l, t)
    return compliance_form(th, h, l, t, ES, 1.0, g * g, 2 * a * g * g)


def m_g(th, h, l, t):
    """Malek & Gibson: rigid joints shorten the bending lengths to lb, hb."""
    _, _, s, c = shorthand(th, h, l, t)
    lb = l - t / (2 * c)
    hb = h - t * (1 - s) / c
    gb = t / lb
    Cp = 1 + (2.4 + 1.5 * NUS) * gb**2
    # same algebra as malek_gibson.m: a = h/l, compliances in lb units
    a = h / l
    D1 = s * s * Cp + c * c * gb**2
    D2 = c * c * Cp + s * s * gb**2 + 2 * (hb / lb) * gb**2
    return {"E1": ES * gb**3 * c / ((a + s) * D1), "E2": ES * gb**3 * (a + s) / (c * D2),
            "nu12": s * c * c * (Cp - gb**2) / ((a + s) * D1),
            "nu21": s * (a + s) * (Cp - gb**2) / D2}


def oh(th, h, l, t):
    """Ohtaki first load step = its small-strain values."""
    e1, s1, n1 = oh_curve(th, h, l, t, ES, 1, 0.002, f_step=0.002)
    e2, s2, n2 = oh_curve(th, h, l, t, ES, 2, -0.002, f_step=0.002)
    return {"E1": s1[1] / e1[1], "E2": s2[1] / e2[1], "nu12": n1[1], "nu21": n2[1]}


MODELS = {"G&A": g_a, "M&E": m_e, "M&G": m_g, "Ohtaki": oh}


# ---- ANSYS ------------------------------------------------------------------
def ansys():
    f = RES / "small_strain.csv"
    if not f.exists():
        return None
    r = pd.read_csv(f).iloc[0]
    return {k: r[f"{k}_MPa" if k.startswith("E") else k] for k in ("E1", "E2", "nu12", "nu21")}


def errors(df):
    ref = df[df.model == "ANSYS"].iloc[0]
    fore = df[(df.zone == "fore") & (df.model != "ANSYS")].set_index("model")
    err = (fore[["E1", "E2", "nu12", "nu21"]] - ref[["E1", "E2", "nu12", "nu21"]]) / \
        ref[["E1", "E2", "nu12", "nu21"]].abs() * 100
    print("\nforefoot cell, % error vs ANSYS:\n", err.astype(float).round(1).to_string())


# ---- figures -----------------------------------------------------------------
# palette: reference categorical slots 1-3 for the linear models, ink for
# Ohtaki (the reference curve), red accent reserved for ANSYS. Models are also
# told apart by line style / marker and direct labels, never colour alone.
INK, MUTED, GRID, SURF = "#1b1b1a", "#6b6a66", "#d9d8d3", "#fcfcfb"
COL = {"G&A": "#2a78d6", "M&E": "#eb6834", "M&G": "#1baf7a", "Ohtaki": INK, "ANSYS": "#e34948"}
DASH = {"G&A": (0, (4, 2)), "M&E": (0, (1.5, 1.5)), "M&G": (0, (6, 2, 1.5, 2)), "Ohtaki": "-"}
MARK = {"G&A": "o", "M&E": "s", "M&G": "D", "Ohtaki": "^", "ANSYS": "*"}


def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Georgia", "DejaVu Serif"],
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 8.5,
        "axes.edgecolor": MUTED, "axes.linewidth": 0.5, "axes.labelcolor": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.facecolor": SURF, "figure.facecolor": SURF, "savefig.facecolor": SURF,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.4,
        "grid.linestyle": (0, (1, 3)), "legend.frameon": False, "legend.fontsize": 7.5,
        "text.color": INK})


def maybe_title(fig, *a, **k):
    if not os.environ.get("REPORT"):
        fig.suptitle(*a, **k)


def caption(fig, n, text):
    if os.environ.get("REPORT"):
        return
    fig.text(0.012, 0.015, f"{n:02d}", fontsize=7, color=COL["ANSYS"], va="bottom")
    fig.text(0.045, 0.015, text, fontsize=7, color=MUTED, style="italic", va="bottom")


def save(fig, name):
    fig.savefig(HERE / f"{name}.png", dpi=220)
    fig.savefig(HERE / f"{name}.svg")
    plt.close(fig)
    print(f"wrote {name}.png / .svg")


def moduli_fig(df):
    """Dot plot, one panel per property: models as rows, zones as dots."""
    from matplotlib.lines import Line2D
    props = [("E1", "E₁  (MPa)"), ("E2", "E₂  (MPa)"),
             ("nu12", "ν₁₂"), ("nu21", "ν₂₁")]
    models = ["G&A", "M&E", "M&G", "Ohtaki"]
    fig, axs = plt.subplots(1, 4, figsize=(10, 2.8), sharey=True)
    for ax, (k, lab) in zip(axs, props):
        for i, m in enumerate(models):
            v = df[df.model == m].set_index("zone")[k]
            ax.plot([v["mid"], v["fore"]], [i, i], color=GRID, lw=1.2, zorder=1)
            ax.scatter(v["fore"], i, s=42, marker=MARK[m], color=COL[m], zorder=3,
                       edgecolor=SURF, linewidth=1)
            ax.scatter(v["mid"], i, s=42, marker=MARK[m], facecolor=SURF, zorder=3,
                       edgecolor=COL[m], linewidth=1)
            same = abs(v["fore"] - v["mid"]) < 0.02
            ax.annotate(f"{v['fore']:.2f}" + ("  both cells" if same else ""), (v["fore"], i),
                        xytext=(0, 7), textcoords="offset points",
                        ha="left" if same else "center", fontsize=6.5, color=INK)
        a = df[df.model == "ANSYS"]
        if len(a):
            ax.axvline(a[k].iloc[0], color=COL["ANSYS"], lw=1, zorder=2)
            ax.annotate("ANSYS", (a[k].iloc[0], len(models) - 0.45), color=COL["ANSYS"],
                        fontsize=6.5, ha="center")
        ax.set_title(lab, loc="left", fontsize=9.5)
        ax.grid(axis="y", visible=False)
        ax.set_ylim(-0.6, len(models) - 0.3)
        if k.startswith("nu"):
            ax.axvline(0, color=MUTED, lw=0.5)
            ax.set_xlim(-1.12, 0.05)
        else:
            ax.set_xlim(0, None)
    axs[0].set_yticks(range(len(models)), models)
    axs[0].tick_params(axis="y", length=0, labelsize=8.5, labelcolor=INK)
    axs[0].invert_yaxis()
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=INK,
                               label="forefoot cell  (t 1.6, l 4)"),
                        Line2D([], [], marker="o", ls="", mfc=SURF, color=INK,
                               label="midfoot / heel cell  (t 1.0, l 5)")],
               loc="upper right", ncol=2, bbox_to_anchor=(0.99, 1.0))
    maybe_title(fig, "Four models, two cells", x=0.012, ha="left", fontsize=12, y=0.98)
    caption(fig, 1, "Small-strain in-plane properties. Bending-only models (G&A, Ohtaki) "
                    "coincide; stretching (M&E) and rigid joints (M&G) soften the cell.")
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    save(fig, "model_compare_moduli")


def draw_cell(ax, p, nx=3, ny=2):
    """Thickened re-entrant walls, nx x ny cells, as a line drawing."""
    from shapely.geometry import LineString
    from shapely.ops import unary_union
    th = np.radians(p.theta)
    c, s, l, h = np.cos(th), np.sin(th), p.l, p.h
    px, py = 2 * l * c, 2 * (h + l * s)
    cell = [((0, 0), (0, h)), ((0, h), (l * c, h + l * s)), ((l * c, h + l * s), (px, h)),
            ((l * c, h + l * s), (l * c, 2 * h + l * s)), ((l * c, 2 * h + l * s), (0, py)),
            ((l * c, 2 * h + l * s), (px, py))]
    lines = [LineString([(a[0] + i * px, a[1] + j * py), (b[0] + i * px, b[1] + j * py)])
             for i in range(nx) for j in range(ny) for a, b in cell]
    from shapely.geometry import box
    walls = unary_union(lines).buffer(p.t / 2, quad_segs=6)
    walls = walls.intersection(box(-p.t / 2, 0, nx * px + p.t / 2, ny * py - abs(l * s) - p.t))
    for g in getattr(walls, "geoms", [walls]):
        ax.fill(*g.exterior.xy, color="#ece9e2", lw=0)
        ax.plot(*g.exterior.xy, color=INK, lw=0.6)
        for r in g.interiors:
            ax.fill(*r.xy, color=SURF, lw=0)
            ax.plot(*r.xy, color=INK, lw=0.6)
    ax.annotate("", (0, -1.5), (px, -1.5), arrowprops=dict(arrowstyle="<->", lw=0.5, color=MUTED))
    ax.text(px / 2, -3.0, f"{px:.2f}", ha="center", fontsize=6.5, color=MUTED)
    ax.annotate("", (-1.5, 0), (-1.5, py), arrowprops=dict(arrowstyle="<->", lw=0.5, color=MUTED))
    ax.text(-2.2, py / 2, f"{py:.1f}", ha="right", va="center", fontsize=6.5, color=MUTED,
            rotation=90)
    ax.set_xlim(-4.5, nx * px + 1)
    ax.set_ylim(-4.5, ny * py + 1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"θ {p.theta:g}°   l {p.l:g}   h {p.h:g}   t {p.t:g}  mm",
                 fontsize=8, color=MUTED, loc="left")


def curves(p):
    """Cell drawing | x2 compression sigma-eps | Ohtaki secant nu vs strain."""
    f = RES / "x2_compression.csv"
    a = pd.read_csv(f) if f.exists() else None
    fig = plt.figure(figsize=(10, 3.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.35, 1])
    draw_cell(fig.add_subplot(gs[0]), p)

    ax = fig.add_subplot(gs[1])
    e = np.linspace(0, 0.3, 60)
    ax.axvspan(15, 20, color=GRID, alpha=0.35, lw=0)
    ax.text(17.5, 2.78, "target\npeak strain", ha="center", fontsize=6.5, color=MUTED)
    eo, so = x2_law(p.theta, p.h, p.l, p.t, ES)
    ymax = 3.0
    ends = {}
    for m in ("G&A", "M&E", "M&G"):
        y = MODELS[m](p.theta, p.h, p.l, p.t)["E2"] * e
        ax.plot(e * 100, y, color=COL[m], lw=1.4, ls=DASH[m])
        ends[m] = y[-1]
    ax.plot(eo * 100, so, color=INK, lw=2, clip_on=True)
    if a is not None:
        ax.plot(a.strain * 100, a.stress_MPa, ls="", marker="o", ms=3.5, color=COL["ANSYS"],
                mec=SURF, mew=0.6)
        ax.text(a.strain.iloc[-1] * 100, min(a.stress_MPa.iloc[-1], ymax * 0.95), "  ANSYS",
                color=COL["ANSYS"], fontsize=7, va="center")
    for m, y in ends.items():
        ax.text(30.4, y, m, color=INK, fontsize=7, va="center", clip_on=False)
    i = np.searchsorted(so, ymax) - 1
    ax.annotate(f"Ohtaki: {so[-1]:.1f} MPa at {eo[-1] * 100:.0f} %",
                (eo[i] * 100, so[i]), xytext=(-150, -60), textcoords="offset points",
                fontsize=7, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED))
    ax.set_xlim(0, 30)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("compression strain  ε₂  (%)")
    ax.set_ylabel("stress  σ₂  (MPa)")
    ax.set_title("Compression along x₂", loc="left")

    ax2 = fig.add_subplot(gs[2])
    eN, _, nN = oh_curve(p.theta, p.h, p.l, p.t, ES, 2, -0.3)
    ax2.plot(-eN[1:] * 100, nN[1:], color=INK, lw=2)
    for m in ("G&A", "M&E", "M&G"):
        v = MODELS[m](p.theta, p.h, p.l, p.t)["nu21"]
        ax2.axhline(v, color=COL[m], lw=1.2, ls=DASH[m])
        dy = {"M&E": 0.06, "M&G": -0.06}.get(m, 0)
        ax2.text(30.4, v + dy, m, fontsize=7, va="center", clip_on=False)
    ax2.text(-eN[-1] * 100 - 0.5, nN[-1], "Ohtaki  ", fontsize=7, va="bottom", ha="right")
    if a is not None and "lateral_strain" in a:
        ax2.plot(a.strain * 100, a.lateral_strain / a.strain.where(a.strain > 0), ls="", marker="o", ms=3.5,
                 color=COL["ANSYS"], mec=SURF, mew=0.6)
    ax2.axhline(0, color=MUTED, lw=0.5)
    ax2.set_xlim(0, 30)
    ax2.set_xlabel("compression strain  ε₂  (%)")
    ax2.set_ylabel("ν₂₁  (secant)")
    ax2.set_title("Auxetic response", loc="left")

    maybe_title(fig, "Forefoot cell under the foot", x=0.012, ha="left", fontsize=12, y=0.98)
    caption(fig, 2, "Linear models stay straight; Ohtaki's inextensible walls straighten and "
                    "lock past ~20 %. Shaded band: the zone-search target." +
                    ("" if a is not None else "  ANSYS points appear once results are added."))
    fig.tight_layout(rect=[0, 0.05, 0.97, 0.99])
    save(fig, "model_compare_curves")


if __name__ == "__main__":
    main()
