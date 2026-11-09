# Re-entrant auxetic midsole: scripts

Clean copy of the working scripts. The original working folder is in `../raw_files/`;
the hand-run ANSYS work stays in `../ansys_stuff/`.

## Layout

```
polished_files/
├── paths.py              folder locations used by every script
├── inputs/               things the scripts read (never written to)
│   ├── midsole_base.step     bulk midsole volume (scaled to 270 mm on read)
│   ├── pseudo_foot.step      toe block, midfoot block, heel sphere
│   └── fukuchi/              one subject of the Fukuchi running dataset (3.5 m/s + static)
├── models/               the cell models
│   ├── matlab/               G&A, M&E, M&G, Ohtaki (.m) + lattice generator
│   ├── lattice_law.py        Ohtaki large-deflection curve in Python (cached)
│   └── model_compare.py      4 models vs ANSYS, comparison figures
├── midsole/              the application, run in this order
│   ├── fill_midsole.py       1. fill the midsole with lattice slices -> STL
│   ├── foot_motion.py        2. Fukuchi stance -> foot angles + force split
│   ├── foot_anim.py             (view any frame / make the GIF)
│   ├── foot_conform.py       3. one frame: posed foot on the slices (beam on lattice)
│   ├── full_stance.py        4. every frame, sink solved from force balance
│   └── zone_search.py        5. pick a lattice per zone (writes zone_pick.csv)
├── outputs/              everything the scripts write (CSV, PNG, STL, GIF, cache)
└── report/               make_figs.py, full_report.html / .pdf, report.html / .pdf, fig/
```

## Running

Run from anywhere; paths resolve from `paths.py`.

```bash
python midsole/foot_motion.py
python midsole/fill_midsole.py          # ~3 min, opens a 3D preview at the end
python midsole/foot_conform.py 30 6.7   # row 30, sink 6.7 mm
python midsole/full_stance.py 1         # every row, ~1 min
python midsole/zone_search.py           # ~5 min; then re-run fill_midsole.py
python models/model_compare.py          # picks up ../ansys_stuff/results/*.csv
python report/make_figs.py              # all report figures, ~4 min
```

`fill_midsole.py` uses `outputs/zone_pick.csv` if it exists (one lattice per zone),
otherwise the single lattice set at the top of the file.

MATLAB scripts are self-contained: open one and press Run (test call at the top).

## Key numbers

- Solid: TPU, Es = 19.3 MPa, nu_s = 0.45 (to be replaced by the measured value).
- Forefoot cell: theta -30 deg, l 4, h 8, t 1.6 mm. Midfoot/heel cell: l 5, h 10, t 1.0 mm.
- ANSYS (forefoot cell): E1 3.45, E2 4.35 MPa, nu12 -0.14, nu21 -0.33; buckles at 0.36 MPa.

Python packages: numpy, scipy, pandas, matplotlib, shapely, trimesh, manifold3d, pyvista, Pillow.
