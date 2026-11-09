"""Turn one Fukuchi running stance (right foot, 3.5 m/s) into pseudo-foot motion.

Out per stance %:  F (N), u (COP position heel 0 -> toe 1), heel / forefoot force,
pitch (rearfoot+midfoot vs ground), ankle (shank vs foot), mtp (toe bend).
Angles in degrees, 0 = static standing trial. Writes foot_motion.csv and a plot.

Run:  python foot_motion.py
"""
import sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
_sys.path[:0] = [str(_ROOT), str(_ROOT / "midsole"), str(_ROOT / "models")]
from paths import ANSYS_RESULTS, INPUTS, OUTPUTS  # noqa: E402

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = INPUTS / "fukuchi"
OUT = OUTPUTS / "foot_motion.csv"
F_ON = 50.0          # N, force threshold for contact
SKIP = 1.0           # s, ignore the first second
MTP_MAX = 35.0       # deg at toe-off (Stefanyshyn & Nigg 1997, ~30-40 in running)
MTP_START = 0.60     # stance fraction where the MTP starts to dorsiflex
ALPHA = 0.3          # EMA factor, 0-1: smaller = smoother, 1 = no smoothing


def main():
    m = pd.read_csv(DATA / "RBDS001runT35markers.txt", sep="\t")
    f = pd.read_csv(DATA / "RBDS001runT35forces.txt", sep="\t")
    s = pd.read_csv(DATA / "RBDS001static.txt", sep="\t")
    f["t"] = (f.Time - 1) / 300.0                 # forces at 300 Hz

    i0, i1 = right_stance(m, f)
    st = f.iloc[i0:i1].reset_index(drop=True)
    mk = m.set_index("Time").reindex(st.t, method="nearest", tolerance=0.004)
    print(f"right stance {st.t.iloc[0]:.3f}-{st.t.iloc[-1]:.3f} s "
          f"({1e3 * (st.t.iloc[-1] - st.t.iloc[0]):.0f} ms), peak {st.Fy.max():.0f} N")

    pitch0, ankle0 = angles(s.mean(numeric_only=True).to_frame().T)
    pitch, ankle = angles(mk)
    u = ema(cop_u(mk, st))
    F = ema(st.Fy.values)
    pct = np.linspace(0, 100, len(st))
    out = pd.DataFrame({
        "pct": pct, "F": F, "u": u,
        "F_heel": (1 - u) * F, "F_fore": u * F,
        "pitch": ema(pitch - pitch0), "ankle": ema(ankle - ankle0),
        "mtp": mtp_curve(pct / 100)})
    out.to_csv(OUT, index=False, float_format="%.3f")
    print(out.iloc[::len(out) // 10].round(1).to_string(index=False))
    plot(out)


def ema(x, a=ALPHA):
    """Exponential moving average, run forward then backward so the
    smoothed curve does not lag behind the data."""
    def one_way(v):
        y = np.empty_like(v, dtype=float)
        y[0] = v[0]
        for i in range(1, len(v)):
            y[i] = a * v[i] + (1 - a) * y[i - 1]
        return y
    x = np.asarray(x, dtype=float)
    return one_way(one_way(x)[::-1])[::-1]


def right_stance(m, f):
    """First full stance after SKIP s where the right foot is on the belt
    (COP side, z, nearer the right MT markers than the left)."""
    on = (f.Fy > F_ON).astype(int).diff().fillna(0).values
    starts, ends = np.where(on == 1)[0], np.where(on == -1)[0]
    mk = m.set_index("Time")
    for a in starts:
        if f.t[a] < SKIP:
            continue
        b = ends[ends > a][0]
        mid = mk.iloc[mk.index.get_indexer([f.t[(a + b) // 2]], method="nearest")[0]]
        zc = f.COPz[(a + b) // 2]
        if abs(zc - mid["R.MT1Z"]) < abs(zc - mid["L.MT1Z"]):
            return a, b
    raise RuntimeError("no right stance found")


def pt(d, name):
    return np.c_[d[name + "X"].values, d[name + "Y"].values, d[name + "Z"].values]


def sag(v, fwd):
    """Angle (deg) of vector v in the sagittal plane, above the forward direction."""
    return np.degrees(np.arctan2(v[:, 1], v[:, 0] * fwd))


def angles(d):
    """pitch: heel->MT line vs ground (+ = toes up).
    ankle: shank vs foot (+ = dorsiflexion)."""
    heel = pt(d, "R.Heel.Bottom")
    mt = (pt(d, "R.MT1") + pt(d, "R.MT5")) / 2
    fwd = np.sign(np.nanmean(mt[:, 0] - heel[:, 0]))      # running along +x or -x
    sh_top = (pt(d, "R.Shank.Top.Lateral") + pt(d, "R.Shank.Top.Medial")) / 2
    sh_bot = (pt(d, "R.Shank.Bottom.Lateral") + pt(d, "R.Shank.Bottom.Medial")) / 2
    pitch = sag(mt - heel, fwd)
    shank = sag(sh_top - sh_bot, fwd)          # ~90 when upright
    ankle = (90 + pitch) - shank               # angle closing between shank and foot
    return pitch, ankle


def cop_u(mk, st):
    """COP projected on the heel->MT line in the ground plane (x, z): 0 heel, 1 MT."""
    heel = pt(mk, "R.Heel.Bottom")[:, [0, 2]]
    mt = ((pt(mk, "R.MT1") + pt(mk, "R.MT5")) / 2)[:, [0, 2]]
    cop = np.c_[st.COPx.values, st.COPz.values]
    d = mt - heel
    u = np.sum((cop - heel) * d, axis=1) / np.sum(d * d, axis=1)
    return np.clip(u, 0, 1)


def mtp_curve(s):
    """MTP dorsiflexion from the literature: ~0 until MTP_START, then a smooth
    rise to MTP_MAX at toe-off."""
    r = np.clip((s - MTP_START) / (1 - MTP_START), 0, 1)
    return MTP_MAX * (1 - np.cos(np.pi * r)) / 2


def plot(o):
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
    ax[0].plot(o.pct, o.F, "k", label="total")
    ax[0].plot(o.pct, o.F_heel, label="heel")
    ax[0].plot(o.pct, o.F_fore, label="forefoot")
    ax[0].set_ylabel("force (N)"); ax[0].legend(); ax[0].grid(True)
    for c in ("pitch", "ankle", "mtp"):
        ax[1].plot(o.pct, o[c], label=c)
    ax[1].set_ylabel("angle (deg)"); ax[1].set_xlabel("stance (%)")
    ax[1].legend(); ax[1].grid(True)
    fig.tight_layout()
    fig.savefig(OUT.with_suffix(".png"), dpi=120)
    plt.show()


if __name__ == "__main__":
    main()
