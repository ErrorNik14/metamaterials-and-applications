# Stage 4: Failure analysis

## 1. Hinge stress
Each ligament carries two kinds of stress.

**Bending.** From Stage 2, the moment along the ligament is constant: $M = EI\,(2\Delta\theta)/L$. That gives

$$\sigma_b=\frac{M\,(w/2)}{I}=\frac{E\,w\,\Delta\theta}{L_h}$$

**Axial.** Cut the lattice vertically. The cut crosses one horizontal ligament per pitch $p_0$, so equilibrium alone (no model needed) gives the ligament force $N = \sigma_x p_0 t$. Then

$$\sigma_a=\frac{\sigma_x\,p_0}{w}$$

**Peak stress** at the outer fibre, on the worst side:

$$\sigma_{max}=K_t(\sigma_b+|\sigma_a|)$$

For small strain, Stage 2 equilibrium gives how much the squares rotate per unit strain:

$$\frac{\Delta\theta}{\varepsilon_x}=\frac{k_a\,b\,p_0(1-\nu)}{2k_ab^2+8k_r}=1.86\ \text{rad per unit strain}$$

For your design, $\sigma_{max}\approx 5590\,\varepsilon_x$ MPa. **93 % of it is bending.**

**First yield** happens at $\varepsilon_x \approx \pm0.89\,\%$. That is when the outer fibre of the ligament reaches the strain of PLA at yield, $\sigma_y/E = 1.43\,\%$. The ligament acts as a strain amplifier: the material in the hinge strains about 1.6 times more than the lattice as a whole.

**Validated** against the FE: the ligaments at the centre of a 13 × 13 array differ from the formula by 0.3 %. Near the **free edges** of your 5 × 5 part, the worst ligament carries **6.5 % more** stress. That edge factor is included in the safe range.

## 2. Locking
- **Open lock (tension):** $p'(\theta)=0$ at 43.1°, which is **+13.4 %** strain. This is a geometric limit.
- **Square contact (compression):** found with a separating-axis overlap test on each square and its 8 neighbours. The squares first touch at θ = −2.7° (the pinwheel has turned through zero), which is **−21 %** strain. The ligament is not included in this test.

## 3. Instability
- **In-plane buckling of the whole array:** found with an FE eigenvalue analysis, $(K+\lambda K_G)\phi=0$. The geometric stiffness was checked against the Euler formula $\pi^2EI/L^2$ to 2e-6.
  - The first mode needs λ = 272, i.e. a strain of −272 %. A strain that large is far outside what a linear analysis can describe. The honest reading is: **no in-plane bifurcation anywhere near the working range**.
  - This makes physical sense. The structure is already rotated ($\theta_0 = 15°$), so compression just carries on closing the mechanism smoothly. There is no symmetric state to branch away from.
- **Euler buckling of one ligament:** $P_{cr}=\pi^2EI/L^2=7369$ N, against about 14 N of axial load at yield, a ratio of 525. The ligament is also too stubby for Euler buckling ($L_h/w=1.25$). It yields long before it could buckle.
- **Out-of-plane buckling:** **not modelled** (the model is 2D).
  - The ligaments are 5 mm deep against 0.8 mm wide, so they are $(t/w)^2 \approx 39$ times stiffer bending out of plane than in plane.
  - That makes out-of-plane buckling unlikely to govern.
  - **Check it in Ansys** with an eigenvalue buckling analysis on a 3D model.

## 4. Safe operating range (SF 1.5 on stress, edge factor ×1.065)

$$\boxed{-0.56\,\% \le \varepsilon_x \le +0.55\,\%}$$

**Yield governs**, by a factor of more than 20 over every geometric or stability limit.

What this means in practice:
- The large auxetic motion seen in demonstration videos (5–10 % strain) is possible, but the PLA hinges then deform **plastically**. The part will not spring all the way back, and it can crack after a few cycles.
- For a reusable demo part:
  - make the ligaments longer, or use a more flexible material such as TPU or PETG;
  - see `S4_design_chart`: an L_h = 3 mm, w = 0.8 mm ligament roughly triples the elastic range.

## Assumptions and limitations
- PLA is treated as linear-elastic up to yield, with no rate or temperature effects.
- 3D-printed parts are anisotropic and weaker between layers. Use `print_knockdown` to account for this.
- $K_t = 1$ assumes no stress concentration at the ligament root. A sharp re-entrant corner raises it, and fillets lower it. **Take $K_t$ from your Ansys results:** the peak stress there divided by this model's nominal stress.
- The safe range is based on **yield**. PLA's strain at break is a few percent, so a part loaded beyond the safe range will deform permanently, not necessarily snap at once.
