# Tensegrity Structures

### Struts that never touch, held by cables that never slacken

A **tensegrity** is a network of compression struts held in place only by tension cables — the struts never
touch each other. It keeps its shape because the cables are pre-tensioned: this *self-stress* stiffens
mechanisms that would otherwise be free to move. Change the self-stress and the stiffness changes with it.

This strand builds on Gilewski et al. (2025), who turned the idea into **2D hexagonal tensegrity-inspired
cells** that can be 3D-printed, derived their elastic properties in closed form, and tested graded lattices
in compression.

<p align="center">
  <img src="images/Screenshot%202026-11-10%20204637.png" width="420" alt="CAD of the hexagonal tensegrity-inspired lattice">
  <br><sub>The project's hexagonal lattice: struts and cables in two parallel planes, tied by transverse stabilisers.</sub>
</p>

---

## Key ideas

- **Two parameters control a cell:** the cable-to-strut stiffness ratio *k* and the self-stress *σ*. At
  σ = √3·k/2 the cell loses all stiffness in one mode — a zero-energy *soft mode*.
- **Orientation is the strongest lever.** Rotating the same hexagonal module by 60° cuts its vertical stiffness
  to about a third and raises its Poisson's ratio about tenfold.
- **Grading places the damage.** In a layered lattice the softest layer fails first, so stiffening the top layers
  pushes failure downward and leaves a *mechanically stable region* above it.

## Test data

`M_hex_1_gra1_RawData_1.csv` is one compression record of the graded lattice with a rotated layer
(M.hex.1.gra1), from the paper's open dataset:

| | This specimen | Paper (mean of five) |
|---|---|---|
| Peak force | 47.2 N at 2.67 mm | 48.57 N at 2.41 mm |
| Linear stiffness | ≈ 25 N/mm | 31.0 N/mm |
| After the peak | plateau near 36 N | — |

The long plateau after the peak is the rotated layer shearing while the rest of the lattice stays intact —
a built-in load limiter.

---

## What's in here

| File | What it is |
|---|---|
| `Gilewski_2025_2D_Mater._12_045026.pdf` | Reference paper: graded 2D tensegrity-inspired lattices, model and tests |
| `book-1.pdf` | *Tensegrity Systems*, Skelton & de Oliveira — background theory |
| `M_hex_1_gra1_RawData_1.csv` | Raw force–displacement record (time, displacement, force; `;`-separated, decimal commas) |
| `images/` | CAD screenshots of the project's 2D lattice and 3D column |

**Reading the CSV:** `pandas.read_csv(path, sep=";", decimal=",", skiprows=[1], encoding="latin-1")` —
the headers are in Polish (*Czas* = time, *Przemieszczenie* = displacement, *Obciążenie* = load).

## Reference

Gilewski, W., Zawadzki, A., Thammavong, M., & Al Sabouni-Zawadzka, A. (2025). 2D graded tensegrity-inspired
metamaterials with hexagonal cells fabricated via additive manufacturing. *2D Materials*, 12, 045026.
[doi:10.1088/2053-1583/ae1620](https://doi.org/10.1088/2053-1583/ae1620)
