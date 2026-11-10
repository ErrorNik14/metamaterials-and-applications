# ANSYS setup: re-entrant lattice block (for model comparison)

Goal: E1, E2, ν12, ν21 at small strain, and the σ–ε curve under x2 compression up to 30 %,
measured the same way the four models define them.

## 1. Geometry (mm): forefoot lattice (the stiffest, most loaded zone; from zone_pick.csv)
θ = −30°, l = 4, h = 8, t = 1.6 (t/l 0.40, h/l 2). Cell repeat: px = 2·l·cosθ = **6.928**, py = 2·(h + l·sinθ) = **12.0**.

One cell, wall centrelines (x = horizontal = x1, y = vertical = x2):

| Wall | From | To |
|---|---|---|
| lower vertical | (0, 0) | (0, 8) |
| lower-left inclined | (0, 8) | (3.464, 6) |
| lower-right inclined | (3.464, 6) | (6.928, 8) |
| middle vertical | (3.464, 6) | (3.464, 14) |
| upper-left inclined | (3.464, 14) | (0, 12) |
| upper-right inclined | (3.464, 14) | (6.928, 12) |

- Copy this cell 4 × 4: x step 6.928, y step 12.
- **Close the right edge:** add one more column of lower verticals (the (0,0)–(0,8) wall only) at x = 4 × 6.928 = 27.71, for each of the 4 rows (y = 0, 12, 24, 36).
- **Close the top (mirror of the bottom):** add verticals of length h = 8 from y = 48 to y = 56 at x = 0, 6.928, 13.86, 20.78, 27.71 (the valley joints of the top row).
- The block is 27.7 mm wide; the walls span y = 0 to 56 and the block is 180° symmetric.
- Thicken each centreline to t = 1.6 (in SpaceClaim: sketch the lines, offset ±0.8, then fill), or draw the walls as 1.6 mm wide rectangles and merge them.
- Add a 1.6 mm solid strip along the bottom (y = −1.6 to 0) and top (y = 56 to 57.6) so the loading is even. Leave the sides open. Use 56 mm as the block height for strain.
- Check: the clear gap inside each bow-tie is h + 2·l·sinθ − t = 2.4 mm.

## 2. Material and mesh
- Linear elastic, **E = 19.3 MPa, ν = 0.45** (for the large run, optional Neo-Hookean with μ = E/(2(1+ν)) = 6.66 MPa).
- 2D **plane stress, thickness 5 mm** (slice thickness B).
- Quad mesh with **at least 4 elements through the wall thickness** (size about 0.3 mm). Refine at the joints.

## 3. Boundary conditions
**Run A: x2, small strain**
- Bottom strip: UY = 0. Pin one node at bottom-centre (UX = 0) only.
- Top strip: UY = −0.56 mm (1 % of 56 mm). Sides free.

**Run B: x1, small strain**
- Same set-up rotated: left edge UX = 0 (one node UY = 0), right edge UX = +0.277 mm (1 % of 27.7). Top and bottom free.
- For run B, put the solid strips on the left and right edges instead of top and bottom.

**Run C: x2, large strain**
- As run A, but UY = −16.8 mm (30 %). **Large deflection ON**, at least 30 substeps (auto time stepping).
- Frictionless self-contact on the wall faces if you can. If not, note where the walls touch.

## 4. What to measure (middle 2 × 2 cells only, to stay away from edge effects)
Pick 4 nodes on the wall centrelines at the corners of the middle 2 × 2 cells (vertical-wall joints):
- ε2 = Δ(UY between top and bottom pair) ÷ their initial distance
- ε1 = Δ(UX between right and left pair) ÷ their initial distance
- σ2 = total top reaction FY ÷ (block width × 5)   (σ1 for run B: FX ÷ (block height × 5))
- E2 = σ2 / ε2, ν21 = −ε1 / ε2 (run A); E1 = σ1 / ε1, ν12 = −ε2 / ε1 (run B)

## 5. Export
Save in `ansys/results/`:
- `small_strain.csv` with one row:
  `E1_MPa,E2_MPa,nu12,nu21`
- `x2_compression.csv`, one row per substep (run C), compression positive:
  `strain,stress_MPa,lateral_strain`
- Also a screenshot of the deformed shape at 0 %, 15 % and 30 % (useful for the report).

## Sanity check (model predictions for this cell)
| | E1 | E2 | ν12 | ν21 |
|---|---|---|---|---|
| G&A | 2.85 | 2.85 | −1.00 | −1.00 |
| M&E | 1.93 | 1.50 | −0.57 | −0.44 |
| M&G | 2.37 | 1.99 | −0.59 | −0.50 |
| Ohtaki (1st step) | 2.86 | 2.85 | −1.00 | −1.00 |

(MPa; from `python model_scripts/model_compare.py`)

If ANSYS is far below all four models, check the mesh through the wall thickness first.
