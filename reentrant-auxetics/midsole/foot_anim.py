"""Cross-section (y = Y_CUT) of the pseudo-foot posed from foot_motion.csv.

Run:  python foot_anim.py 40        -> frame for row 40, saved as foot_anim_row40.png
      python foot_anim.py 40 anim   -> also the whole stance as foot_anim.gif,
                                       row 40 highlighted in the force plot
Pose: heel+midfoot rigid, rotate by -pitch about the MTP hinge (forefoot striker,
the MTP stays planted); toes rotate a further mtp about the same hinge;
shank drawn from the ankle at its measured angle. Midsole is undeformed.
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import trimesh
from matplotlib.animation import FuncAnimation, PillowWriter

HERE = OUTPUTS
Y_CUT = 60.0
ANKLE = np.array([231.2, 82.9])
HEEL_MID = np.array([186.6, 74.2])
MTP = np.array([63.1, 52.0])
SHANK_LEN = 150.0


def main():
    row = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    mo = pd.read_csv(HERE / "foot_motion.csv")
    row = min(max(row, 0), len(mo) - 1)
    shapes = sections()
    fig, ax, art = setup(shapes, mo)
    draw(art, shapes, mo, row)
    out = HERE / f"foot_anim_row{row}.png"
    fig.savefig(out, dpi=110)
    print(f"row {row}: {mo.pct[row]:.1f} % stance -> {out.name}")
    if len(sys.argv) > 2 and sys.argv[2] == "anim":
        frames = range(0, len(mo), 2)
        anim = FuncAnimation(fig, lambda k: draw(art, shapes, mo, k), frames=frames)
        anim.save(HERE / "foot_anim.gif", writer=PillowWriter(fps=15))
        print("wrote foot_anim.gif")
    plt.show()


def sections():
    """x-z outlines of each body and of the midsole at y = Y_CUT."""
    t = (INPUTS / "pseudo_foot.step").read_text(encoding="latin-1")
    pts = {int(i): [float(v) for v in x.split(",")]
           for i, x in re.findall(r"#(\d+)=CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", t)}
    V = np.array([pts[int(i)] for i in re.findall(r"=VERTEX_POINT\('[^']*',#(\d+)\)", t)])
    toe = trimesh.convex.convex_hull(V[V[:, 0] < 62])
    mid = trimesh.convex.convex_hull(V[(V[:, 0] > 62) & (V[:, 0] < 195)])
    sph = trimesh.creation.icosphere(4, 37.5)
    sph.apply_translation([231.173, 69.357, 82.897])
    sole = trimesh.load(HERE / "midsole_filled.stl")

    def cut(m):
        s = m.section(plane_origin=[0, Y_CUT, 0], plane_normal=[0, 1, 0])
        return [e[:, [0, 2]] for e in s.discrete]

    return {"toe": cut(toe), "mid": cut(mid), "heel": cut(sph), "sole": cut(sole)}


def rot(p, c, a_deg):
    """Rotate points p about c by a_deg (counter-clockwise in the x-z view)."""
    a = np.radians(a_deg)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return (p - c) @ R.T + c


def pose(mo, k):
    rear = -mo.pitch[k]                  # heel lifts when pitch < 0
    toe = rear - mo.mtp[k]               # dorsiflexion turns the toes clockwise
    shank = 90 + mo.pitch[k] - mo.ankle[k]   # angle above forward (-x)
    return rear, toe, shank


def setup(shapes, mo):
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(12, 7), height_ratios=[3, 1])
    for e in shapes["sole"]:
        ax.plot(e[:, 0], e[:, 1], color="#BA7517", lw=0.4)
    art = {}
    for name, c in [("toe", "#534AB7"), ("mid", "#0F6E56"), ("heel", "#993C1D")]:
        art[name] = [ax.fill(e[:, 0], e[:, 1], color=c, alpha=0.3, ec=c, lw=1.2)[0]
                     for e in shapes[name]]
    art["shank"], = ax.plot([], [], color="#444", lw=6, solid_capstyle="round")
    art["hinges"], = ax.plot([], [], "ko", ms=5)
    art["title"] = ax.set_title("")
    ax.set_aspect("equal"); ax.set_xlim(-20, 330); ax.set_ylim(-5, 260)
    ax.set_xlabel("x (mm)  toe -> heel"); ax.set_ylabel("z (mm)")
    ax2.plot(mo.pct, mo.F_heel, label="heel"); ax2.plot(mo.pct, mo.F_fore, label="forefoot")
    ax2.set_xlabel("stance (%)"); ax2.set_ylabel("N"); ax2.legend(loc="upper right")
    ax2.grid(True)
    art["cursor"] = ax2.axvline(0, color="k", lw=1)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig, ax, art


def draw(art, shapes, mo, k):
    rear, toe, shank = pose(mo, k)
    for p, e in zip(art["toe"], shapes["toe"]):
        p.set_xy(rot(e, MTP, toe))
    for name in ("mid", "heel"):
        for p, e in zip(art[name], shapes[name]):
            p.set_xy(rot(e, MTP, rear))
    a = rot(ANKLE, MTP, rear)
    d = np.array([-np.cos(np.radians(shank)), np.sin(np.radians(shank))])
    art["shank"].set_data([a[0], a[0] + SHANK_LEN * d[0]], [a[1], a[1] + SHANK_LEN * d[1]])
    hm = rot(HEEL_MID, MTP, rear)
    art["hinges"].set_data([a[0], hm[0], MTP[0]], [a[1], hm[1], MTP[1]])
    art["cursor"].set_xdata([mo.pct[k]] * 2)
    art["title"].set_text(f"row {k}: {mo.pct[k]:.0f} % stance   F {mo.F[k]:.0f} N   "
                          f"pitch {mo.pitch[k]:.1f}  ankle {mo.ankle[k]:.1f}  mtp {mo.mtp[k]:.1f} deg")


if __name__ == "__main__":
    main()
