function [kg, T, kl] = fe_frame_k(E, A, I, x1, x2)
%FE_FRAME_K  2D Euler-Bernoulli frame element stiffness in global axes.
%
%   [kg, T, kl] = fe_frame_k(E, A, I, x1, x2)
%   DOFs per node: [u, v, rotation];  element DOFs [u1 v1 r1 u2 v2 r2].
%   x1, x2: node coordinates [x y] (mm).  E [MPa], A [mm^2], I [mm^4].
%
%   Local stiffness = bar (axial) + Euler-Bernoulli beam (cubic deflection):
%     axial  EA/L [1 -1; -1 1]
%     bending EI/L^3 [12 6L -12 6L; 6L 4L^2 -6L 2L^2; ...]
%   Global: kg = T' kl T, T rotates global -> local axes.
d = x2(:) - x1(:);
L = norm(d);
c = d(1) / L;  s = d(2) / L;

ka = E * A / L;
kb = E * I / L^3;
kl = zeros(6);
kl([1 4], [1 4]) = ka * [1 -1; -1 1];
kl([2 3 5 6], [2 3 5 6]) = kb * [ 12    6*L   -12    6*L;
                                  6*L  4*L^2 -6*L  2*L^2;
                                 -12   -6*L    12   -6*L;
                                  6*L  2*L^2 -6*L  4*L^2];
R = [c s 0; -s c 0; 0 0 1];
T = blkdiag(R, R);
kg = T.' * kl * T;
end
