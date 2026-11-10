# Rotating-squares auxetic: one-page theory summary

**Design:** PLA, a = 10 mm, L_h = 1 mm, w = 0.8 mm, t = 5 mm, θ₀ = 15°, 5 × 5 squares.

## 1. Kinematics (rigid squares, pin hinges)
Neighbouring squares rotate by −θ and +θ. Closing the vector loop through one hinge gives the centre-to-centre pitch:

$$p(\theta)=a(\cos\theta+\sin\theta)\quad\text{(the same in x and y)}$$

Because the cell grows equally in both directions, $\varepsilon_x=\varepsilon_y$ and **ν = −1 exactly**, for any θ and any size of strain.

**Locks:**
- **Open:** $dp/d\theta = 0$ at θ = 45°. From θ₀ = 15° this is a strain of $2/\sqrt3-1 = +15.5\,\%$.
- **Closed:** θ = 0, a strain of $\sqrt{2/3}-1 = -18.4\,\%$.

**Porosity:** $\sin2\theta/(1+\sin2\theta)$, which reaches 50 % at 45°.

## 2. Compliant hinges (energy method)
Each ligament acts as two springs:
- bending: $k_r=EI/L_h$, with $I = tw^3/12$;
- stretching: $k_a=EA/L_h$.

Linearising the hinge gap gives $\delta=p_0\varepsilon-b\,\Delta\theta$. Here $b=a(\cos\theta_0-\sin\theta_0)$ is the lever arm.

The cell energy is $U=4k_a(p_0e-b\Delta\theta)^2+4k_ap_0^2d^2+16k_r\Delta\theta^2$. Minimising it over Δθ (rotation and stretching act as springs in series) gives:

$$\nu=-\frac{1}{1+8\rho},\quad \rho=\frac{w^2}{12b^2}\qquad E_{eff}=\frac{8k_r}{t\,b^2(1+8\rho)}$$

**Your design:**
- ν = −0.9915;
- E_eff = 23.7 MPa, about 148 × softer than PLA;
- open lock with real ligaments: 43.1°, a strain of +13.4 %.

## 3. Independent check (own beam FE)
The FE uses 2D Euler–Bernoulli frame elements. The squares are rigid stars of beams, and each ligament is 4 PLA beam elements.

| Model | ν | E_eff |
|---|---|---|
| Periodic cell (infinite lattice) | −0.99153 | 23.691 MPa: matches theory to 1e-5 |
| Finite 5 × 5 array | −0.9904 | 23.68 MPa: about 0.1 % edge effect |

## 4. Failure
The peak ligament stress is $\sigma=K_t(E w\Delta\theta/L_h+\sigma_x p_0/w)$. At about 5600 MPa per unit strain, **93 % of it is bending**.

| Limit | Strain |
|---|---|
| First yield | ±0.9 % |
| Open lock | +13.4 % |
| Square contact | −21 % |
| In-plane buckling | far beyond any useful strain (linear eigenvalue gives −272 %) |
| **Safe range** (SF 1.5, edge factor 1.065) | **−0.56 % to +0.55 %** |

**Yield governs.** Large auxetic motion in PLA needs plastic hinges. Longer ligaments (L_h ≈ 3 mm) roughly triple the elastic range.

## 5. Assumptions and limitations
- 2D plane stress.
- Rigid squares.
- Linear-elastic beam ligaments with ideal clamped roots.
- PLA treated as isotropic.
- No out-of-plane buckling; check it in Ansys.
- $K_t = 1$.
- Ansys solid models include root flexibility, so expect a softer E_eff there. Calibrate with `root_length_factor` using `compare_with_ansys`.
