function C = rs_compliant_props(a, Lh, w, t, E, theta0, root_factor)
%RS_COMPLIANT_PROPS  Small-strain properties of rotating squares with
%   compliant (beam) ligaments, from the energy method.
%
%   C = rs_compliant_props(a, Lh, w, t, E, theta0, root_factor)
%   Inputs may be arrays of equal size (or scalars) -> used for sweeps.
%     a [mm], Lh [mm], w [mm], t [mm], E [MPa], theta0 [rad]
%     root_factor : effective ligament length = Lh + root_factor*w
%
%   MODEL (see docs/stage2_derivation.md):
%     - squares rigid; each ligament = Euler-Bernoulli beam, clamped in
%       the two squares, with
%           rotational spring  k_r = E I / L     (I = t w^3/12)
%           axial spring       k_a = E A / L     (A = t w)
%     - per 2x2 cell: 8 ligaments (4 horizontal, 4 vertical)
%     - square centres follow the average strain (affine); the square
%       rotation dtheta is the one internal unknown, found by minimising
%       the strain energy.
%
%   RESULTS (closed form):
%       rho   = k_r / (k_a b^2) = w^2 / (12 b^2),   b = a(cos th0 - sin th0)
%       nu    = -1 / (1 + 8 rho)
%       E_eff = 8 k_r / ( t b^2 (1 + 8 rho) )

L   = Lh + root_factor .* w;          % effective bending length
I   = t .* w.^3 / 12;                 % second moment of the ligament section
A   = t .* w;                         % ligament cross-section area

C.k_r = E .* I ./ L;                  % [N mm / rad] moment per relative rotation
C.k_a = E .* A ./ L;                  % [N / mm]     force per axial stretch

% Lever arm: how much the hinge gap opens per radian of square rotation
C.b   = a .* (cos(theta0) - sin(theta0));          % = dp/dtheta at theta0
C.p0  = a .* (cos(theta0) + sin(theta0)) + Lh;     % as-printed pitch
C.Lcell = 2 * C.p0;                                % 2x2 cell side

% Competition between the two ways a hinge can deform:
%   rotation (auxetic, nu = -1)  vs  stretching (nu = 0)
C.rho   = C.k_r ./ (C.k_a .* C.b.^2);
C.nu    = -1 ./ (1 + 8 * C.rho);
C.E_eff = 8 * C.k_r ./ (t .* C.b.^2 .* (1 + 8 * C.rho));   % [MPa]
C.E_ratio = C.E_eff ./ E;

% Limits that the code is tested against
C.E_eff_rotation_only = 8 * C.k_r ./ (t .* C.b.^2);         % k_a -> infinity
end
