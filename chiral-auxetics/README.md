# Chiral Auxetics

### Rotating rings, winding ligaments, and a morphing airfoil

A **hexachiral honeycomb** is a lattice of rings, each joined to six neighbours by ligaments that leave it
*tangentially*. Squeeze it and every ring turns the same way, winding the ligaments up around it — so the
whole lattice contracts in every direction at once. In the ideal case its Poisson's ratio is exactly −1.

This strand asks two questions:

1. **How close does a real, finite lattice get to that ideal?** The unit cell is derived from first principles,
   compared with three published models, and checked against a 3D ANSYS compression specimen.
2. **What can the mechanism do in a structure?** A chiral core is mapped into an Eppler 420 airfoil section
   and compared with a straight-web core of the same mass, to see how easily the trailing edge can be morphed.

📄 **Read the report:** [Hexachiral_Honeycomb_Project_Report.docx](Hexachiral_Honeycomb_Project_Report.docx) ·
companion [airfoil application paper](Eppler420_Chiral_Honeycomb_Application_Paper.docx)

---

## Part 1 · The lattice

| Model | E/Eₛ (r 5, R 25, t 1 mm) | ν |
|---|---|---|
| Bending-only limit (derived here) | 1.80 × 10⁻³ | −1 (limit) |
| Spadoni & Ruzzene (axial + bending) | 1.74 × 10⁻³ | −0.954 |
| Prall & Lakes (printed formula) | 3.02 × 10⁻³ | −1 (limit) |
| **ANSYS finite coupon (R01)** | **2.89 × 10⁻³** | **−0.819** |

- The derivation reproduces the classic result: a stiffness that scales with the **cube of ligament thickness**,
  and ν → −1 in the slender limit.
- The finite 3D coupon is clearly auxetic in all seven ANSYS runs, but **less** auxetic and **stiffer** than the
  ideal lattice. Mesh refinement changes the modulus by at most 1.2 %, so the gap comes from the finite specimen
  (edges, end constraints, ring flexibility), not the mesh.
- Two literature formulas (Prall & Lakes, Bornengo et al.) are checked line by line, and inconsistencies in
  their prefactors are identified.

## Part 2 · The Eppler 420 airfoil

A 330 mm chord section with a chiral middle core, against a straight-web core of almost identical mass
(79.8 g vs 80.6 g), both actuated from inside the section:

| At 0.50 mm of trailing-edge travel | Chiral | Straight web |
|---|---|---|
| Actuator stroke needed | **0.47 mm** (−21 %) | 0.60 mm |
| Actuator force | 90.8 N | 28.4 N |
| Mechanical work | 0.021 J | 0.008 J |

The chiral core turns **less stroke into more tip motion**, but needs about 2.5 × the work to do it. Its peak
stresses stay below 6061-T6 yield at the comparison point and cross it at about 0.61 mm of stroke.

---

## What's in here

| File | What it is |
|---|---|
| `Hexachiral_Honeycomb_Project_Report.docx` | Main report: derivation, models, ANSYS coupon, airfoil study |
| `Eppler420_Chiral_Honeycomb_Application_Paper.docx` | Companion paper on the airfoil application |
| `chiral_models.m` | MATLAB function for the three closed-form lattice models |
| `Hexachiral_Honeycomb_ANSYS_Results.xlsx` | Coupon runs R01–R07: probe readings and calculated properties |
| `Eppler420_Aerofoil_ANSYS_Comparison.xlsx` | Airfoil ramps, core masses and matched-travel comparison |
| `Chiral Honeycomb sample.STEP` | CAD of the printed hexachiral sample |

## Running the models

In MATLAB, with this folder as the current folder:

```matlab
chiral_models()      % press Enter three times for r = 5, R = 25, t = 1 mm
```

It prints E/Eₛ and ν for the three models and plots both against ligament thickness.
