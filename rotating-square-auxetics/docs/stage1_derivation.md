# Stage 1 — Ideal rotating squares: derivation

## Assumptions
1. **Rigid squares**: the squares don't stretch or bend. All deformation is rotation.
2. **Ideal pin hinges** at the corners: no rotational stiffness, no friction, no size. The hinge sits exactly at the corner.
3. **2D, in-plane** motion only. The structure is periodic: a 2 × 2 block of squares repeats.
4. **Pure geometry**: no forces appear. Stage 1 is a *mechanism* (one degree of freedom, θ), so it predicts shape, not stiffness.
5. Strains are measured from the **as-printed** angle θ₀.

**Limitation:** the real part has ligaments of finite length L_h and width w, so it is not a perfect pin. Stage 2 adds that.

## Setup
- Every square has side *a*.
- Neighbouring squares rotate in opposite directions, like a checkerboard:
  - type A squares rotate by −θ;
  - type B squares rotate by +θ.
- θ = 0: the squares touch edge-to-edge (closed). As θ grows, diamond-shaped holes open.
- Each corner is pinned to exactly one neighbour (a "pinwheel" pattern).

## 1. Pitch between neighbouring squares (x direction)
Take square A at the origin and its right neighbour B. They share one hinge:
- the hinge is at A's top-right corner, local position (a/2)(1, 1);
- the same hinge is at B's top-left corner, local position (a/2)(−1, 1).

The rotation matrix is

$$R(\alpha)=\begin{bmatrix}\cos\alpha & -\sin\alpha\\ \sin\alpha & \cos\alpha\end{bmatrix}$$

The hinge must be at the same point when reached from A and from B (vector loop closure):

$$\mathbf c_B=\mathbf c_A + R(-\theta)\,\tfrac a2\!\begin{bmatrix}1\\1\end{bmatrix} - R(+\theta)\,\tfrac a2\!\begin{bmatrix}-1\\1\end{bmatrix}$$

Multiply out:
- $R(-\theta)(1,1)^T=(\cos\theta+\sin\theta,\ \cos\theta-\sin\theta)$
- $R(\theta)(-1,1)^T=(-\cos\theta-\sin\theta,\ \cos\theta-\sin\theta)$

The y-parts cancel. That leaves

$$\boxed{p_x(\theta)=a(\cos\theta+\sin\theta)}$$

## 2. Pitch in the y direction
Do the same with the top neighbour:
- the hinge is at A's top-left corner, (a/2)(−1, 1);
- on the square above, it is the bottom-left corner, (a/2)(−1, −1).

This time the x-parts cancel, and

$$\boxed{p_y(\theta)=a(\cos\theta+\sin\theta)}$$

The unit cell (2 × 2 squares) measures $L_x = 2p_x$ by $L_y = 2p_y$.

A useful identity: $\cos\theta+\sin\theta=\sqrt2\,\sin(\theta+45^\circ)$.

## 3. Stretch ratio and strain
Measured from the as-printed angle θ₀:

$$\lambda_x=\lambda_y=\frac{\cos\theta+\sin\theta}{\cos\theta_0+\sin\theta_0},\qquad \varepsilon_x=\varepsilon_y=\lambda-1$$

## 4. Poisson's ratio
**Incremental (tangent) value:**

$$\nu_{xy}=-\frac{d\varepsilon_y}{d\varepsilon_x}=-\frac{dp_y/d\theta}{dp_x/d\theta}=-\frac{a(\cos\theta-\sin\theta)}{a(\cos\theta-\sin\theta)}=\boxed{-1}$$

**Secant value:** $-\varepsilon_y/\varepsilon_x=-1$ too, because $\varepsilon_x=\varepsilon_y$ exactly.

So ν = −1 at every angle and for any size of strain, not just small strains. It is a consequence of the geometry alone. It does not depend on a, θ₀ or the material.

## 5. Area and porosity
A 2 × 2 cell holds 4 squares of area a² each. The cell area is

$$A_{cell}=L_xL_y=4a^2(\cos\theta+\sin\theta)^2=4a^2(1+\sin2\theta)$$

From that:
- solid fraction $= 1/(1+\sin2\theta)$;
- porosity $= \sin2\theta/(1+\sin2\theta)$, which is 50 % at 45°;
- areal ratio $A/A_0=(1+\varepsilon)^2$.

Stretching the auxetic increases its area a lot. A normal material (ν > 0) barely changes area.

## 6. Locking angles
- **Open lock:** $dp/d\theta=a(\cos\theta-\sin\theta)=0$ gives θ = 45°. The pitch is at its maximum, $a\sqrt2$. Pulling further cannot lengthen the cell; past 45° the mechanism would close again. At exactly 45°, $d\varepsilon/d\theta=0$, so the incremental ν is 0/0 and undefined. That is why the plotted line has a gap there.
- **Closed lock:** at θ = 0 the squares touch edge-to-edge and the pitch is a. Compression cannot go further without crushing the squares.

Strain limits starting from θ₀:

$$\varepsilon_{lock}^{+}=\frac{\sqrt2}{\cos\theta_0+\sin\theta_0}-1,\qquad \varepsilon_{lock}^{-}=\frac{1}{\cos\theta_0+\sin\theta_0}-1$$

For θ₀ = 15°: $\cos15°+\sin15°=\sqrt6/2$. So
- $\varepsilon^{+}=2/\sqrt3-1=+15.47\,\%$
- $\varepsilon^{-}=\sqrt{2/3}-1=-18.35\,\%$

## 7. Inverse kinematics (used in Stages 4–5)

$$\theta(\varepsilon)=\arcsin\!\left[\frac{(1+\varepsilon)(\cos\theta_0+\sin\theta_0)}{\sqrt2}\right]-45^\circ$$

## How the code checks this
`rs_build_array` places every square by **chaining hinges** and never uses the formula. Rows are built left to right; every other square is built from the one below it. The horizontal hinges in rows 2 and above are therefore never imposed. All of them still close (the gap is about 1e-14 mm), which shows the pinwheel mechanism is compatible.

The pitch measured from that array matches $a(\cos\theta+\sin\theta)$ to about 1e-14 mm. The ν computed from it by finite differences is −1 to within about 1e-13.
