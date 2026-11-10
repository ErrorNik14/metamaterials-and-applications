# Metamaterials and Applications

**Team "Mettamaterials" · Mechanics of Solids Project 2026–27**

Most materials get thinner when you stretch them and fatter when you squash them. *Auxetic* structures do
the opposite: squeeze one and it pulls inward, stretch it and it widens. That behaviour comes entirely from
geometry — the shape of the cells, not the material they are printed in — which makes it something an
engineer can design.

This repository collects the team's work on four ways of building that geometry. Each strand takes one
cell type from first principles (analytical models), checks it against finite-element simulation, and then
puts it to work in a real application.

<p align="center">
  <img src="reentrant-auxetics/report/fig/barrel_3d.png" width="720"
       alt="A conventional honeycomb bulging outward under load next to a re-entrant honeycomb drawing inward">
  <br><sub>Same compression, opposite response: a conventional honeycomb barrels out (left), a re-entrant one draws in (right).</sub>
</p>

---

## The four strands

| | Strand | Cell | Application | Report |
|---|---|---|---|---|
| ⬡ | **[Re-entrant auxetics](reentrant-auxetics/)** | Bow-tie honeycomb whose walls fold inward | A running-shoe midsole, loaded by a measured stance | [Markdown](reentrant-auxetics/report/report.md) · [PDF](reentrant-auxetics/report/Re-entrant_Auxetics_Report.pdf) |
| ◎ | **[Chiral auxetics](chiral-auxetics/)** | Hexachiral lattice: rings joined by tangent ligaments that wind up as the rings turn | A morphing Eppler 420 airfoil trailing edge | [PDF](chiral-auxetics/Hexachiral_Honeycomb_Report.pdf) |
| ◇ | **[Rotating squares](rotating-square-auxetics/)** | Rigid squares linked at their corners, counter-rotating as a mechanism | Stage-by-stage MATLAB analysis, own FE solver and printable CAD | [HTML report](rotating-square-auxetics/Rotating_Squares_Report.html) |
| △ | **[Tensegrity structures](tensegrity-structures/)** | Struts held in a web of tension cables | — | [PDF](tensegrity-structures/Tensegrity_Lattices_Report.pdf) |

The **[experimental setup](experimental-setup/)** folder is reserved for the shared compression tests of the
printed samples; it will be filled in once printing is complete.

## How each strand is organised

Every strand follows the same arc, so the reports can be read side by side:

1. **Mechanism** — why the geometry produces a negative Poisson's ratio.
2. **Models** — closed-form predictions of stiffness and Poisson's ratio from the literature.
3. **Simulation** — ANSYS finite-element runs to test those predictions.
4. **Application** — the cell used in a design, and what it gains or costs there.

## Getting started

- **Just want the results?** Open a strand's report from the table above.
- **Want to run the code?** Each strand's own README lists its scripts in run order.
  The re-entrant and rotating-square strands are Python and MATLAB respectively; the chiral strand ships a
  MATLAB model file and the ANSYS result workbooks.
