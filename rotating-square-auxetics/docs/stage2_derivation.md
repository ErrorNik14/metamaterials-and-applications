# Stage 2: Compliant hinges (energy method)

## Assumptions
1. The squares are still **rigid**. All the flexibility is in the ligaments.
2. Each ligament is a beam of length $L_h$ and width $w$ (the bending direction), with out-of-plane depth $t$. It is clamped into both squares, and its axis follows the line joining the two corners it connects.
3. The ligaments are **linear-elastic Euler–Bernoulli beams** with small rotations and plane stress, so we use $E$ (set `wide_beam` in the parameters to use $E/(1-\nu^2)$ instead).
4. Square **centres move with the average strain** (an affine displacement). Each square also rotates by an unknown extra angle $\Delta\theta$. The pattern has 4-fold symmetry, which rules out any other internal motion when the load is along x or y.
5. Small strain about the as-printed angle $\theta_0$. The only exception is the nonlinear curve in §6.

**Limitations:**
- Squares deform a little in reality. This model ignores that.
- The ligament root is not a perfect clamp, and real ligaments with $L_h/w \approx 1$ are stubby. Both effects make real hinges *softer* than this model predicts. The parameter `root_length_factor` exists to calibrate against Ansys.
- Contact between squares and yield are not checked here. They come in Stage 4.

## 1. Ligament springs
The second moment of area is $I = tw^3/12$ and the cross-section area is $A = tw$.

When the two squares rotate by $\mp\Delta\theta$ about the ligament mid-point, the two ends of the ligament rotate by $\mp\Delta\theta$ with no relative sideways movement. With the slope-deflection equations, the end moments are then

$$M_1=\frac{EI}{L}(4\varphi_1+2\varphi_2),\qquad \varphi_{1,2}=\mp\Delta\theta.$$

This gives $M_1=-M_2=2EI\Delta\theta/L$: a constant moment, so pure bending with no shear force. That means:

$$k_r=\frac{EI}{L_h}\ \text{(moment per relative rotation } 2\Delta\theta),\qquad k_a=\frac{EA}{L_h}\ \text{(axial)}$$

Because there is no shear force, Euler–Bernoulli is acceptable here even for a short ligament. The neglected root flexibility is the larger error.

## 2. Kinematics of one hinge
Take two horizontal neighbours with the ligament along x. Linearise the gap between the two arm tips (each tip is a corner plus half a ligament):

$$\delta_h = p_0\,\varepsilon_x - b\,\Delta\theta,\qquad \delta_v = p_0\,\varepsilon_y - b\,\Delta\theta$$

Here:
- $b=a(\cos\theta_0-\sin\theta_0)$ is the **lever arm**, equal to $dp/d\theta$.
- $p_0=a(\cos\theta_0+\sin\theta_0)+L_h$ is the as-printed pitch.

The half-ligament arms cancel out of $b$.

A pure mechanism has $\delta=0$. That forces $\varepsilon_x=\varepsilon_y=b\Delta\theta/p_0$, which is Stage 1's result $\nu=-1$.

## 3. Strain energy of a 2 × 2 cell (8 ligaments: 4 horizontal, 4 vertical)
$$U = 4\cdot\tfrac12k_a\delta_h^2 + 4\cdot\tfrac12k_a\delta_v^2 + 8\cdot\tfrac12k_r(2\Delta\theta)^2$$

Split the strain into a volumetric part $e=(\varepsilon_x+\varepsilon_y)/2$ and a deviatoric part $d=(\varepsilon_x-\varepsilon_y)/2$:

$$U = 4k_a(p_0e-b\Delta\theta)^2 + 16k_r\Delta\theta^2 + 4k_ap_0^2d^2$$

Now minimise over $\Delta\theta$, the internal unknown. The first two terms behave as springs in **series**:
- $K_1=4k_a$ (stretching);
- $K_2=16k_r/b^2$ (rotation).

$$U_{min} = \frac{K_1K_2}{K_1+K_2}p_0^2e^2 + 4k_ap_0^2d^2$$

## 4. Effective elastic constants
Divide by the cell volume $V=4p_0^2t$ and compare with $W=(C_{11}+C_{12})e^2+(C_{11}-C_{12})d^2$. Define

$$\rho=\frac{k_r}{k_ab^2}=\frac{w^2}{12\,b^2}.$$

Then:

$$\boxed{\nu_{xy}=\frac{C_{12}}{C_{11}}=-\frac{1}{1+8\rho}=-\frac{1}{1+\tfrac23(w/b)^2}}$$

$$\boxed{E_{eff}=\frac{8k_r}{t\,b^2(1+8\rho)}=\frac{2}{3}\,\frac{E\,w^3}{L_h\,a^2(\cos\theta_0-\sin\theta_0)^2}\cdot\frac{1}{1+8\rho}}$$

What the formulas say:
- **Thin hinges** ($\rho\to0$) give $\nu\to-1$, which recovers Stage 1.
- **Thick hinges** (larger $\rho$): stretching the ligaments competes with rotating them. Stretching alone gives $\nu=0$, so $\nu$ rises from −1 towards 0.
- $\nu$ depends only on $w/b$: not on $E$, $t$ or $L_h$.
- $E_{eff}\propto E\,w^3/(L_h a^2)$. It is independent of $t$, because stiffness and area both scale with $t$.
- As $\theta_0\to45°$, $b\to0$ and $E_{eff}\to\infty$: the lattice is locked.

**Defaults** ($a=10$, $L_h=1$, $w=0.8$ mm, $\theta_0=15°$, PLA):
- $\rho=1.07\times10^{-3}$, $\nu=-0.9915$
- $E_{eff}=23.7$ MPa, about 148 × softer than solid PLA

## 5. Kinematics with real ligament arms (large rotation)
$$p(\theta)=a(\cos\theta+\sin\theta)+L_h\cos(\theta-\theta_0)$$

The open lock is where $p'(\theta)=0$:
- **43.1°**, at a strain of **+13.4 %** (for the defaults).
- Stage 1 gave 45° and +15.5 %, because it had no ligament arms.

## 6. Nonlinear stress–strain (rotation mode only)
The energy of a cell is $U=16k_r(\theta-\theta_0)^2$. Only $\sigma_x$ does work, so $\sigma_x\,L_{y0}t\,dL_x=dU$. That gives

$$\sigma_x(\theta)=\frac{8k_r(\theta-\theta_0)}{p_0\,t\,p'(\theta)}$$

The stress rises without limit as $p'\to0$ at the lock. Its slope at $\theta_0$ equals $E_{eff}$ for thin hinges (the code checks this).

## How the code checks this
`rs_energy_check` builds the periodic 2 × 2 cell directly from rotation matrices, with the same hinge table as Stage 1. It computes each ligament's gap Jacobian **numerically** and assembles the 3 × 3 stiffness matrix. It then solves for uniaxial stress.

It never uses the algebra above. It also includes a transverse beam spring, $12EI/L^3$, which the closed form leaves out.

The two routes agree to about $10^{-10}$ relative, both for your design and for a deliberately thick-hinged design.
