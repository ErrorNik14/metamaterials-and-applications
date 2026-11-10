# Rotating-squares auxetic: MATLAB analysis

📄 **Read the report:** [PDF](Rotating_Squares_Report.pdf) · [HTML](Rotating_Squares_Report.html)

**To run:** open `main_rotating_squares.m` and press F5. This folder must be the current folder in MATLAB.

**To change any input:** edit `params_default.m`. It holds every number in the project: geometry, material, printer, plot style and tolerances.

| Folder / file | What it is |
|---|---|
| `params_default.m` | All inputs, one struct `P`, units mm / N / MPa |
| `main_rotating_squares.m` | Runs the stages in order, one section per stage |
| `src/common/` | Figure style, saving, drawing helpers |
| `src/stage1_kinematics/` | Ideal rigid-square model and its checks |
| `src/stage2_compliant/` | Beam-ligament energy model, sweeps, nonlinear curve |
| `src/stage3_fe/` | Own frame FE solver: finite N×N arrays + periodic cell |
| `src/stage4_failure/` | Hinge stress, yield, locking/contact, FE buckling, safe range |
| `src/stage5_visual/` | Layout drawing, state-at-angle model, GIF/MP4 animation |
| `rs_explorer_app.m` | **Slider app**: run `rs_explorer_app` in the Command Window |
| `stage5_live_demo.m` | Open as Live Script, put a slider on `theta_deg` |
| `results/animation/` | `S5_rotating_squares.gif` and `.mp4` for the slides |

Every figure is saved three ways: PNG (300 dpi), PDF (vector, sharpest in
PowerPoint), and FIG (editable in MATLAB).
| `docs/` | Derivations (markdown + LaTeX) |
| `results/figures/` | PNG (300 dpi) + `.fig` files for the slides |
| `results/data/` | CSV exports (Stage 6) |

| `src/stage6_export/` | DXF/STL/SolidWorks exports, print check, results CSV |
| `results/cad/` | `S6_lattice_outline.dxf` (sketch), `S6_lattice.stl` (print), SolidWorks parameters/equations |
| `compare_with_ansys.m` | Overlay Ansys results; calibrates `root_length_factor` |
| `docs/theory_summary_one_page.md` | One-page theory |
| `docs/viva_questions.md` | 8 likely questions with answers |
| `docs/presentation_checklist.md` | What to show and say |

**SolidWorks:** File → Open → `results/cad/S6_lattice_outline.dxf` → *Import to a new part as 2D sketch* (units: mm) → Extruded Boss = t.
Or drive your own model from `S6_solidworks_equations.txt` (Tools → Equations → *Link to external file*).

**Ansys:** copy `results/data/ansys_results_TEMPLATE.csv` to `ansys_results.csv`, fill in one row per load step, run `compare_with_ansys`.

Status: all six stages done; `main_rotating_squares` prints "all stages passed = 1".
