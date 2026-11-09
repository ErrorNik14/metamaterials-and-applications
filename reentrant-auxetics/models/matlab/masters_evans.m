%% Test: run this file to see an example
clear; clc; close all;
me = m_e(-20, 1.3*7, 7, 1.1, 19.3, 0);

%% Function
function me = m_e(theta_deg, h, l, t, Es, rh)
% M_E  In-plane elastic properties of a re-entrant honeycomb,
% Masters & Evans model: walls bend, stretch, and (optionally) hinge.
%
%   me = m_e(theta_deg, h, l, t, Es, rh)
%
%   theta_deg : wall angle from x1 in degrees (negative = re-entrant)
%   h, l, t   : vertical wall, inclined wall, wall thickness (mm)
%   Es        : Young's modulus of the solid material (MPa)
%   rh        : hinging, Kf/Kh (0 = rigid joints, no hinging)
%
%   me.E1, me.E2   : effective Young's moduli along x1 and x2 (MPa)
%   me.nu12, nu21  : Poisson's ratios (nu12 = load along x1)
%
% Study sheet form, compliances in units of l^3/(Es b t^3):
%   Cp = 1 + rh  (bending + hinging), Ca = g^2 (stretching), Cv = 2 a g^2
% Derived from the paper's mechanisms, not copied from its equations.

% ---- shorthand -----------------------------------------------------------
a = h / l;
g = t / l;
s = sind(theta_deg);
c = cosd(theta_deg);

% ---- wall compliances -----------------------------------------------------
Cp = 1 + rh;         % perpendicular to the inclined wall
Ca = g^2;            % along the inclined wall
Cv = 2 * a * g^2;    % vertical wall, x2 load only

% ---- x1 loading -----------------------------------------------------------
D1   = s^2 * Cp + c^2 * Ca;
E1   = Es * g^3 * c / ((a + s) * D1);
nu12 = s * c^2 * (Cp - Ca) / ((a + s) * D1);

% ---- x2 loading -----------------------------------------------------------
D2   = c^2 * Cp + s^2 * Ca + Cv;
E2   = Es * g^3 * (a + s) / (c * D2);
nu21 = s * (a + s) * (Cp - Ca) / D2;

% ---- pack the results -----------------------------------------------------
me = struct('E1', E1, 'E2', E2, 'nu12', nu12, 'nu21', nu21);

fprintf('M&E  theta %g, h/l %.2f, t/l %.3f, Kf/Kh %g\n', theta_deg, a, g, rh);
fprintf('  E1 = %.4g MPa   E2 = %.4g MPa\n', E1, E2);
fprintf('  nu12 = %.3f   nu21 = %.3f\n', nu12, nu21);
end
