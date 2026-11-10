"""Compare the ANSYS plane-stress grid with the analytical models, then redo the
dip / energy search with ANSYS-calibrated stiffness.

  1. per cell: ANSYS E2, nu21 vs G&A, M&E, M&G, Ohtaki (first step)  -> compare_cells.csv
  2. calibration: Ohtaki curve scaled per cell so its initial E2 equals ANSYS E2
     (shape kept, stiffness corrected), full stance re-run for every grid cell
     -> grid_eval_ansys.csv (dip, energy, peak strain), same as the analytical search
  3. figures: grid_vs_models.png, grid_pareto_compare.png

Run:  python compare_grid.py        (after run_all.bat)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PF = HERE.parents[1]
sys.path[:0] = [str(PF), str(PF / "midsole"), str(PF / "models")]
from paths import OUTPUTS  # noqa: E402

INK, MUTED, GRID, PAPER, ACCENT = "#1b1b1a", "#6b6a66", "#d9d8d3", "#fcfcfb", "#c0392b"
COL = {"G&A": "#2a78d6", "M&E": "#eb6834", "M&G": "#1baf7a", "Ohtaki": INK}
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.linewidth": 0.5, "axes.facecolor": PAPER, "figure.facecolor": PAPER,
                     "savefig.facecolor": PAPER, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.4, "grid.linestyle": (0, (1, 3)), "legend.frameon": False})


def load():
    import model_compare as mc
    cells = pd.read_csv(HERE / "grid_plane_cells.csv")
    res = pd.read_csv(HERE / "grid_plane_results.csv", skipinitialspace=True); res["id"] = res.id.astype(int)
    d = cells.drop(columns=["E2", "nu21"]).merge(res[["id", "E2_MPa", "nu21", "E2_global_MPa"]], on="id")
    d = d.rename(columns={"E2_MPa": "E2_ansys", "nu21": "nu21_ansys"})
    for m in ("G&A", "M&E", "M&G", "Ohtaki"):
        f = mc.MODELS[m]
        if m == "Ohtaki":
            from lattice_law import x2_law
            v = [x2_law(r.theta, r.h, r.l, r.t, 19.3) for r in d.itertuples()]
            d["E2_Ohtaki"] = [s[1] / e[1] if len(e) > 1 else np.nan for e, s in v]
            d["nu21_Ohtaki"] = -1.0 * np.sign(d.theta) * 0 + np.nan        # nu from first step below
            for i, r in enumerate(d.itertuples()):
                d.loc[i, "nu21_Ohtaki"] = mc.oh(r.theta, r.h, r.l, r.t)["nu21"]
        else:
            vals = [f(r.theta, r.h, r.l, r.t) for r in d.itertuples()]
            d[f"E2_{m}"] = [v["E2"] for v in vals]
            d[f"nu21_{m}"] = [v["nu21"] for v in vals]
    d["tl"] = d.t / d.l
    return d


def fig_models(d):
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2))
    ax = axs[0]
    lim = [0.02, max(d.E2_ansys.max(), d["E2_G&A"].max()) * 1.4]
    ax.plot(lim, lim, color=MUTED, lw=0.6)
    for m, mk in (("G&A", "o"), ("M&E", "s"), ("M&G", "D"), ("Ohtaki", "^")):
        ax.scatter(d.E2_ansys, d[f"E2_{m}"], s=16, marker=mk, color=COL[m], label=m, edgecolor=PAPER, lw=0.4)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("ANSYS E₂ (MPa)"); ax.set_ylabel("model E₂ (MPa)")
    ax.set_title("Stiffness: model vs ANSYS", fontsize=9, loc="left"); ax.legend(fontsize=7)
    ax = axs[1]
    ax.plot([-1.1, 0.05], [-1.1, 0.05], color=MUTED, lw=0.6)
    for m, mk in (("G&A", "o"), ("M&E", "s"), ("M&G", "D")):
        ax.scatter(d.nu21_ansys, d[f"nu21_{m}"], s=16, marker=mk, color=COL[m], label=m, edgecolor=PAPER, lw=0.4)
    ax.set_xlim(-1.1, 0.05); ax.set_ylim(-1.1, 0.05)
    ax.set_xlabel("ANSYS ν₂₁"); ax.set_ylabel("model ν₂₁")
    ax.set_title("Poisson's ratio: model vs ANSYS", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(HERE / "grid_vs_models.png", dpi=220); plt.close(fig)


def recalibrated_search(d):
    """Full stance per cell with Ohtaki's curve scaled to the ANSYS stiffness."""
    f = HERE / "grid_eval_ansys.csv"
    if f.exists():
        return pd.read_csv(f)
    from fill_midsole import read_volume
    from foot_conform import ES, foot_bodies
    from full_stance import F_MIN, sink_solve
    from lattice_law import x2_law
    mo = pd.read_csv(OUTPUTS / "foot_motion.csv")
    feet, (o, zt) = foot_bodies(), read_volume()
    rows = []
    for r in d[d.zone.isna() | (d.zone == "")].itertuples():
        e, s = x2_law(r.theta, r.h, r.l, r.t, ES)
        k = r.E2_ansys / (s[1] / e[1])                  # stiffness correction, shape kept
        law = (e, s * k)
        F, S, P = [], [], []
        for i in range(0, len(mo), 4):
            m = mo.loc[i]
            if m.F < F_MIN:
                continue
            sk, tops = sink_solve(feet, o, zt, law, m)
            F.append(m.F); S.append(sk); P.append(max(t["eps"].max() for t in tops))
        F, S = np.array(F), np.array(S)
        km = int(np.argmax(S))
        rows.append({"id": r.id, "theta": r.theta, "l": r.l, "h": r.h, "t": r.t, "tl": r.tl,
                     "E2_ansys": r.E2_ansys, "k": k, "dip": S.max(),
                     "energy_J": np.trapezoid(F[:km + 1], S[:km + 1]) / 1000, "peak_strain": max(P)})
        print(f"  cell {r.id:2d}  k {k:4.2f}  dip {S.max():5.2f}  energy {rows[-1]['energy_J']:.2f}")
    out = pd.DataFrame(rows); out.to_csv(f, index=False)
    return out


def fig_pareto(d, g):
    a = pd.read_csv(OUTPUTS / "grid_eval.csv")
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    for df, lab, c, mk in ((a, "analytical (M&E-ranked, Ohtaki law)", MUTED, "o"),
                           (g, "ANSYS-calibrated", ACCENT, "s")):
        ok = df[df.dip < 19.99]
        ax.scatter(ok.dip, ok.energy_J, s=18, marker=mk, color=c, label=lab, edgecolor=PAPER, lw=0.4)
    ax.set_xlabel("dip (mm)"); ax.set_ylabel("energy absorbed (J)")
    ax.legend(fontsize=7)
    ax.set_title("Dip vs energy: analytical vs ANSYS-calibrated", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(HERE / "grid_pareto_compare.png", dpi=220); plt.close(fig)


def main():
    d = load()
    d.to_csv(HERE / "compare_cells.csv", index=False)
    show = ["id", "theta", "l", "h", "t", "zone", "E2_ansys", "E2_G&A", "E2_M&E", "E2_M&G", "E2_Ohtaki",
            "nu21_ansys", "nu21_G&A", "nu21_M&E", "nu21_M&G"]
    print(d[show].round(3).to_string(index=False))
    for m in ("G&A", "M&E", "M&G", "Ohtaki"):
        err = (d[f"E2_{m}"] / d.E2_ansys - 1) * 100
        print(f"{m:7s} E2 error vs ANSYS: median {err.median():+6.1f} %, range {err.min():+.0f} to {err.max():+.0f} %")
    fig_models(d)
    g = recalibrated_search(d)
    fig_pareto(d, g)
    print("wrote compare_cells.csv, grid_vs_models.png, grid_eval_ansys.csv, grid_pareto_compare.png")


if __name__ == "__main__":
    main()
