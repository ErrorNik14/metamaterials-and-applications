%% Test: run this file to see an example
clear; clc; close all;
ga = g_a(-20, 1.3*7, 7, 1.1, 19.3);

%% Function
function ga = g_a(theta_deg, h, l, t, Es)
% GIBSON_ASHBY  In-plane elastic properties of a re-entrant honeycomb,
% Gibson & Ashby model (walls bend only; no stretching, hinging or shear).
%
%   ga = g_a(theta_deg, h, l, t, Es)
%   ga = g_a(-20, 1.3*7, 7, 1.1, 19.3)    % example, TPU 90A
%
%   theta_deg : wall angle from x1 in degrees (negative = re-entrant)
%   h, l, t   : vertical wall, inclined wall, wall thickness (mm)
%   Es        : Young's modulus of the solid material (MPa)
%
%   ga.E1, ga.E2   : effective Young's moduli along x1 and x2 (MPa)
%   ga.nu12, nu21  : Poisson's ratios (nu12 = load along x1)
%   ga.G12         : shear modulus (MPa)
%   ga.rho         : relative density (solid fraction of the area)
%
% Same formulas as the study sheet, with Cp = 1, Ca = 0, Cv = 0.

% ---- shorthand -----------------------------------------------------------
a = h / l;            % alpha, aspect ratio
g = t / l;            % thickness ratio
s = sind(theta_deg);
c = cosd(theta_deg);

% ---- Young's moduli -------------------------------------------------------
E1 = Es * g^3 * c / ((a + s) * s^2);
E2 = Es * g^3 * (a + s) / c^3;

% ---- Poisson's ratios -----------------------------------------------------
nu12 = c^2 / ((a + s) * s);
nu21 = (a + s) * s / c^2;

% ---- shear modulus (G&A book, eq. 4.16 form; from memory, not checked) ----
G12 = Es * g^3 * (a + s) / (a^2 * (1 + 2*a) * c);

% ---- relative density -----------------------------------------------------
% wall area in one cell / cell area (vertical wall counted as single thickness)
rho = g * (a + 2) / (2 * c * (a + s));

% ---- pack the results -----------------------------------------------------
ga = struct('E1', E1, 'E2', E2, 'nu12', nu12, 'nu21', nu21, ...
            'G12', G12, 'rho', rho);

fprintf('G&A  theta %g, h/l %.2f, t/l %.3f\n', theta_deg, a, g);
fprintf('  E1 = %.4g MPa   E2 = %.4g MPa   G12 = %.4g MPa\n', E1, E2, G12);
fprintf('  nu12 = %.3f   nu21 = %.3f   rho = %.3f\n', nu12, nu21, rho);
end
