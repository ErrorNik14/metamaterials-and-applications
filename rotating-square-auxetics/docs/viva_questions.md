# Likely viva questions with short answers

**1. Why is the ideal Poisson's ratio exactly −1?**
Rigid squares joined at their corners have one degree of freedom, the rotation θ. The pitch in x and in y is the same function, a(cos θ + sin θ). So ε_x = ε_y always, and ν = −ε_y/ε_x = −1. This comes from geometry alone, so it holds for any material, any size and any amount of strain.

**2. Why does ν move away from −1 with real hinges?**
A ligament can bend (the squares rotate, giving ν = −1) or stretch (giving ν = 0). The two act like springs in series, and the mix is set by ρ = w²/12b². This gives ν = −1/(1+8ρ). Thicker ligaments, or a smaller lever arm b, mean more stretching, so ν moves towards 0.

**3. Where does the stiffness come from, and why is it so low?**
Almost all of it comes from bending the ligaments: E_eff ≈ 8k_r/(t b²), with k_r = EI/L_h ∝ w³. The lattice turns a small hinge moment into a large deformation through the lever arm b. That makes it about 150 × softer than solid PLA.

**4. How did you validate your results?**
In three independent ways:
- an analytical derivation;
- a numerical energy minimisation built directly from rotation matrices;
- my own beam finite-element code, with a periodic cell for the infinite lattice and finite arrays for the real part.

They agree to 1e-5 for the infinite lattice. The FE code was itself verified on a cantilever (PL³/3EI), on rigid-body modes, and on Euler buckling (π²EI/L²).

**5. Why does your FE model of a finite array differ slightly from the theory?**
Free edges. Squares at the edge have fewer ligaments, and the soft rotation mechanism lets that disturbance spread into the array. The periodic FE model removes the edges and matches the theory. The finite 5 × 5 array is about 0.1 % less auxetic, which is too small to measure.

**6. What limits the strain of your printed part?**
Yield of the ligaments, at about ±0.9 % lattice strain, not locking (+13 %) or buckling. The hinge works as a strain amplifier: its outer fibre reaches PLA's yield strain (σ_y/E ≈ 1.4 %) when the lattice is still under 1 %. With a safety factor of 1.5 and the edge factor, the safe range is about ±0.55 %.

**7. Why doesn't the lattice buckle in compression?**
It is printed already rotated (θ₀ = 15°), so compression simply continues the existing mechanism smoothly. There is no symmetric state for it to suddenly branch away from. The FE eigenvalue analysis confirms no in-plane bifurcation near the working range. Each ligament is also far too short (L/w = 1.25) to Euler-buckle; it yields first. Out-of-plane buckling was not modelled; that is checked in Ansys.

**8. What would you change to improve the design, and what does Ansys add?**
To improve it: use longer, thinner ligaments (L_h ≈ 3 mm with w = 0.8 mm roughly triples the elastic range), or a more flexible material such as TPU or PETG.

Ansys adds what the beam model leaves out:
- the flexibility of the ligament root;
- the real stress concentration K_t;
- deformation of the squares themselves;
- 3D effects.

We use it to calibrate `root_length_factor` and K_t.
