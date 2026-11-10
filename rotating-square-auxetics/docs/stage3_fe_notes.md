# Stage 3: Frame finite-element check

## What the model is
- **Element:** a 2D Euler–Bernoulli frame element, the same kind as Ansys BEAM188 in 2D. Each node has 3 degrees of freedom: $u, v, \varphi$.
- **Local stiffness:** an axial bar plus a beam with cubic deflection,

$$k_l=\begin{bmatrix}\tfrac{EA}{L}&0&0&-\tfrac{EA}{L}&0&0\\ 0&\tfrac{12EI}{L^3}&\tfrac{6EI}{L^2}&0&-\tfrac{12EI}{L^3}&\tfrac{6EI}{L^2}\\ 0&\tfrac{6EI}{L^2}&\tfrac{4EI}{L}&0&-\tfrac{6EI}{L^2}&\tfrac{2EI}{L}\\ \cdots \end{bmatrix},\qquad k_g=T^Tk_lT$$

- **Squares:** each square is a stiff "star" of 4 beams running from its centre to its corners. They are $10^3\times$ stiffer than PLA, so effectively rigid. Results are unchanged from $10^2$ to $10^7$.
- **Ligaments:** each ligament is 4 PLA beam elements (w × t section), welded to the corner nodes. The results are the same with 1 element, because the cubic beam element is exact for end-loaded beams.
- **Two models:**
  1. **Finite N × N array with grips.** The centres of the left column have $u_x=0$ and those of the right column have $u_x=\varepsilon(N-1)p_0$. Everything else is free. Strain is measured by a least-squares fit to the centres of the interior squares, and $E_{eff}$ comes from the grip force.
  2. **Periodic cell.** A 2 × 2 cell whose right edge wraps onto its left edge and its top onto its bottom. The average strains $\varepsilon_x, \varepsilon_y$ are extra unknowns: the far end of a wrapped element moves by $u+\varepsilon L$. This model is an **infinite lattice**, the same thing Stage 2 describes.

## Results (default design)

| Model | ν | E_eff (MPa) |
|---|---|---|
| Stage 1, rigid squares + pins | −1 | — (mechanism) |
| Stage 2, energy method | −0.99154 | 23.691 |
| FE periodic (tension = compression) | −0.99153 | 23.691 |
| FE 5 × 5 array | −0.99036 | 23.682 |
| FE 13 × 13 array | −0.99003 | 23.684 |

## What this shows
1. **Verification** (is the solver right?):
   - the cantilever matches $PL^3/3EI$ to 1e-14;
   - a free 2 × 2 lattice has exactly 3 rigid-body modes;
   - the grip forces balance.
2. **Validation** (is the Stage 2 theory right?): the periodic FE agrees with Stage 2 to 1e-5 at your design, and to within 5e-4 at w = 3 mm. The two are completely independent: beam elements on one side, a hand derivation on the other.
3. **Finite size:**
   - Real arrays with free edges are about 0.1 % less auxetic than the infinite lattice. That is far too small to measure on a printed part.
   - The finite-array value does not settle quickly as N grows. A likely reason is that the rotation mechanism makes non-uniform deformations cheap, so boundary disturbances reach far into the array.
   - I have not proven this. What is certain is that the effect is about 0.1 %.
4. **Tension equals compression exactly**, because the solver is **linear** (small displacements). Real differences between tension and compression (locking, contact, buckling) are nonlinear. They are handled in Stage 4.

## Limitations
- Linear and small-strain. Valid for strains up to about 1–2 %. Beyond that, use the Stage 2 nonlinear curve or a nonlinear Ansys run.
- Squares are effectively rigid, and the ligament roots are ideal clamps. A 2D or 3D solid model in Ansys includes the root flexibility. Expect Ansys to give a **softer** $E_{eff}$, and possibly a ν a bit further from −1. That difference is real physics this beam model cannot capture, not an error in either model.
