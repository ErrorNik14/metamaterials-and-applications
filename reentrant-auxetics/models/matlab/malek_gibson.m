%% Test: run this file to see an example
clear; clc; close all;
mg = m_g(-20, 1.3*7, 7, 1.1, 19.3, 0.45);

%% Function
function mg = m_g(theta_deg, h, l, t, Es, nus)
% M_G  In-plane elastic properties of a re-entrant honeycomb,
% Malek & Gibson (2015) model: bending + stretching + shear, with solid
% (rigid) joints, so only the shorter lengths l_b and h_b bend.
%
%   mg = m_g(theta_deg, h, l, t, Es, nus)
%
%   theta_deg : wall angle from x1 in degrees (negative = re-entrant)
%   h, l, t   : vertical wall, inclined wall, wall thickness (mm)
%   Es, nus   : Young's modulus (MPa) and Poisson's ratio of the solid
%
%   mg.E1, mg.E2   : effective Young's moduli along x1 and x2 (MPa)
%   mg.nu12, nu21  : Poisson's ratios (nu12 = load along x1)
%   mg.lb, mg.hb   : bending lengths of inclined and vertical walls (mm)
%
% Paper validated for theta > 0 only; theta < 0 is an extrapolation.

% ---- shorthand -----------------------------------------------------------
a = h / l;
s = sind(theta_deg);
c = cosd(theta_deg);

% ---- bending lengths (rigid joints) ---------------------------------------
lb = l - t / (2*c);
hb = h - t * (1 - s) / c;
gb = t / lb;                         % thickness ratio of the bending part

% ---- wall compliances, units of lb^3/(Es b t^3) ---------------------------
Cp = 1 + (2.4 + 1.5*nus) * gb^2;     % bending + shear
Ca = gb^2;                           % stretching along the inclined wall
Cv = 2 * (hb/lb) * gb^2;             % vertical wall stretching, x2 only

% ---- x1 loading -----------------------------------------------------------
D1   = s^2 * Cp + c^2 * Ca;
E1   = Es * gb^3 * c / ((a + s) * D1);
nu12 = s * c^2 * (Cp - Ca) / ((a + s) * D1);

% ---- x2 loading -----------------------------------------------------------
D2   = c^2 * Cp + s^2 * Ca + Cv;
E2   = Es * gb^3 * (a + s) / (c * D2);
nu21 = s * (a + s) * (Cp - Ca) / D2;

% ---- pack the results -----------------------------------------------------
mg = struct('E1', E1, 'E2', E2, 'nu12', nu12, 'nu21', nu21, 'lb', lb, 'hb', hb);

fprintf('M&G  theta %g, h/l %.2f, t/l %.3f, nu_s %.2f\n', theta_deg, a, t/l, nus);
fprintf('  l_b = %.2f mm   h_b = %.2f mm\n', lb, hb);
fprintf('  E1 = %.4g MPa   E2 = %.4g MPa\n', E1, E2);
fprintf('  nu12 = %.3f   nu21 = %.3f\n', nu12, nu21);
end
