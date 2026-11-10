# Presentation checklist (MATLAB part)

All files are in `results/figures/`. Each figure exists as PNG, PDF (vector) and FIG.

## Story, one slide each
- [ ] **What a rotating-squares auxetic is.** Show `S1_geometry_snapshots` and play `results/animation/S5_rotating_squares.mp4` on a loop.
- [ ] **Ideal theory, ν = −1.** Show the loop-closure equation p = a(cos θ + sin θ) with `S1_poisson_ratio` and `S1_strain_vs_theta` (the locks).
- [ ] **Real hinges.** Show ν = −1/(1+8ρ) and E_eff with `S2_nu_vs_hinge_width` and `S2_parametric_sweeps`.
- [ ] **Validation.** Show `S3_nu_comparison`, the S3 comparison table (`results/data/S3_comparison_table.csv`) and `S3_deformed_array`. Say "three independent methods agree to 1e-5".
- [ ] **Failure.** Show `S4_hinge_stress` and `S4_safe_window`, and state the safe range **±0.55 %**.
- [ ] **Design improvement.** Show `S4_design_chart` (longer ligaments mean a larger elastic range).
- [ ] **MATLAB vs Ansys.** Show `S6_ansys_comparison` after you run `compare_with_ansys`.
- [ ] **Printed part.** Show `S6_print_outline` and the print-check results next to a photo of the part.

## Before the talk
- [ ] Run `main_rotating_squares` once more and confirm "all stages passed = 1".
- [ ] Fill `results/data/ansys_results.csv` from your Ansys runs, using the same grips as Stage 3. Then run `compare_with_ansys`.
- [ ] If Ansys gives a different peak hinge stress, set `P.fail.Kt` = (Ansys peak) / (MATLAB nominal) and re-run Stage 4.
- [ ] Measure the printed part: mark the square centres, take photos at 0 % and a small strain, and compute ε_y/ε_x. Compare with −0.99.
- [ ] Test the live demo: `rs_explorer_app(P, S4.edge_factor)`.

## Things to say explicitly
- [ ] The assumptions: 2D plane stress, rigid squares, Euler–Bernoulli ligaments, linear PLA up to yield.
- [ ] The limitations: no out-of-plane buckling in MATLAB; K_t = 1 unless calibrated; FDM anisotropy (`print_knockdown`).
- [ ] What is verified and what is not. The beam models are verified against each other. Agreement with real parts still depends on Ansys and testing.
