function R = rs_energy_check(a, Lh, w, t, E, theta0, root_factor)
%RS_ENERGY_CHECK  Independent numerical route to nu and E_eff (Stage 2 check).
%
%   Builds one periodic 2x2 cell from rotation matrices (no hand algebra),
%   connects the 8 ligaments as springs, and finds the stiffness matrix of
%   the cell numerically:
%       unknowns  q = [eps_x; eps_y; dtheta]
%       energy    U(q) = sum over ligaments of  1/2 k_a (axial gap)^2
%                                             + 1/2 k_t (transverse gap)^2
%                                             + 1/2 k_r (relative rotation)^2
%       stiffness H = sum k * J'J,  J = d(gap)/dq by central differences
%   Then uniaxial stress sigma_x (sigma_y = 0, no moment on dtheta):
%       H q = [sigma_x * V; 0; 0]   ->   nu = -eps_y/eps_x,  E = sigma/eps_x
%
%   Assumption shared with rs_compliant_props: square centres move with the
%   average strain (justified by the 4-fold symmetry of the pattern).

L   = Lh + root_factor * w;
I   = t * w^3 / 12;
k_r = E * I / L;
k_a = E * t * w / L;
k_t = 12 * E * I / L^3;      % clamped-clamped beam, transverse

% --- Geometry of the as-printed 2x2 cell ------------------------------------
G.a = a;  G.Lh = Lh;  G.theta0 = theta0;
G.CORNER = (a/2) * [1 1; -1 1; -1 -1; 1 -1];      % TR TL BL BR (local)
% Hinge corner table, same as rs_build_array (type A = 1, type B = 2)
G.HINGE = struct('right', [1 4], 'left', [3 2], 'top', [2 1], 'bottom', [4 3]);
G.sgn = [-1 +1];
p0 = a * (cos(theta0) + sin(theta0)) + Lh;        % pitch incl. ligament
[I2, J2] = meshgrid(1:2, 1:2);                     % column i, row j
cells  = [I2(:) J2(:)];
G.typ  = 1 + mod(cells(:,1) + cells(:,2), 2);
G.cen0 = (cells - 1) * p0;
Lcell  = 2 * p0;
V      = Lcell^2 * t;

% --- Ligament list (8 per cell, periodic) ------------------------------------
% each row: square, neighbour, side, other side, direction, periodic shift
lig = cell(0, 6);
for s = 1:4
    i = cells(s,1); j = cells(s,2);
    nb = find(cells(:,1) == mod(i,2)+1 & cells(:,2) == j);
    lig(end+1, :) = {s, nb, 'right', 'left', [1 0], [Lcell*(i==2) 0]}; %#ok<AGROW>
    nb = find(cells(:,1) == i & cells(:,2) == mod(j,2)+1);
    lig(end+1, :) = {s, nb, 'top', 'bottom', [0 1], [0 Lcell*(j==2)]}; %#ok<AGROW>
end

% --- Stiffness matrix from numerical Jacobians ---------------------------------
dq_step = 1e-6;
H = zeros(3);
max_gap = 0;
for k = 1:size(lig, 1)
    dir  = lig{k, 5};
    n_ax = dir(:);                 % ligament axis
    n_tr = [-dir(2); dir(1)];      % perpendicular in-plane
    Jax = zeros(1,3); Jtr = zeros(1,3); Jrot = zeros(1,3);
    for m = 1:3
        dq = zeros(3,1); dq(m) = dq_step;
        [gp, pp] = ligament_gap(G, lig(k,:),  dq);
        [gm, pm] = ligament_gap(G, lig(k,:), -dq);
        Jax(m)  = n_ax.' * (gp - gm) / (2*dq_step);
        Jtr(m)  = n_tr.' * (gp - gm) / (2*dq_step);
        Jrot(m) = (pp - pm) / (2*dq_step);
    end
    H = H + k_a*(Jax.'*Jax) + k_t*(Jtr.'*Jtr) + k_r*(Jrot.'*Jrot);
    max_gap = max(max_gap, norm(ligament_gap(G, lig(k,:), zeros(3,1))));
end

sigma = 1;                           % any value: the problem is linear
q = H \ [sigma * V; 0; 0];
R.nu    = -q(2) / q(1);
R.E_eff = sigma / q(1);
R.H     = H;
R.max_gap_as_printed = max_gap;      % must be ~0: geometry is consistent
end

% ============================================================================
function [g, phi] = ligament_gap(G, L, q)
% Gap between the two arm tips of one ligament, and its relative rotation.
[s, nb, side, other, dir, wrap] = L{:};
F  = diag([1 + q(1), 1 + q(2)]);              % average deformation
xA = arm_tip(G, s,  side,  dir,  q);
xB = arm_tip(G, nb, other, -dir, q) + F * wrap(:);
g  = xB - xA;
phi = (G.sgn(G.typ(nb)) - G.sgn(G.typ(s))) * q(3);
end

function x = arm_tip(G, s, side, dir, q)
% Corner of square s plus half a ligament, in the deformed state q.
rot = @(al) [cos(al) -sin(al); sin(al) cos(al)];
F   = diag([1 + q(1), 1 + q(2)]);
al0 = G.sgn(G.typ(s)) * G.theta0;              % as-printed rotation
al  = G.sgn(G.typ(s)) * (G.theta0 + q(3));     % current rotation
arm = rot(-al0) * dir(:) * (G.Lh/2);           % ligament points along 'dir' when printed
loc = G.CORNER(G.HINGE.(side)(G.typ(s)), :).' + arm;
x   = F * G.cen0(s, :).' + rot(al) * loc;
end
