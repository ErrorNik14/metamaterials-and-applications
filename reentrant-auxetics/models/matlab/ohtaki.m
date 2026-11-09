%% Test: run this file to see an example
clear; clc; close all;
oh = oh_curve(-20, 1.3*7, 7, 1.1, 19.3, 1, 0.05);   % x1 tension up to 5 %

figure;
subplot(1,2,1); plot(oh.eps, oh.sigma, 'LineWidth', 1.5); grid on
xlabel('\epsilon (load direction)'); ylabel('\sigma (MPa)'); title('Stress-strain')
subplot(1,2,2); plot(oh.eps(2:end), oh.nu(2:end), 'LineWidth', 1.5); grid on
xlabel('\epsilon (load direction)'); ylabel('\nu (secant)'); title('Poisson''s ratio')

%% Function
function oh = oh_curve(theta_deg, h, l, t, Es, dir, eps_end)
% OH_CURVE  Large-deflection response of a re-entrant honeycomb,
% Ohtaki et al. (2004) model: each half inclined wall is an elastica
% (bending only, but large rotations allowed), so nu changes with strain.
%
%   oh = oh_curve(theta_deg, h, l, t, Es, dir, eps_end)
%
%   theta_deg : wall angle from x1 in degrees (negative = re-entrant)
%   h, l, t   : vertical wall, inclined wall, wall thickness (mm)
%   Es        : Young's modulus of the solid material (MPa)
%   dir       : load direction, 1 = x1, 2 = x2
%   eps_end   : strain to stop at; positive = tension, negative = compression
%
%   oh.eps   : strain in the load direction
%   oh.sigma : stress in the load direction (MPa)
%   oh.nu    : secant Poisson's ratio (-lateral strain / load strain)
%
% Half wall, s from the joint (clamped at angle theta) to mid-wall (no moment):
%   phi'' = fx sin(phi) - fy cos(phi),  phi(0) = theta,  phi'(1) = 0
% Lengths in units of l/2, force f = F (l/2)^2 / (Es I).

% ---- shorthand -----------------------------------------------------------
a  = h / l;
g  = t / l;
th = deg2rad(theta_deg);
s  = sin(th);
c  = cos(th);
sgn = sign(eps_end);

% ---- load steps -----------------------------------------------------------
f_step = 0.05;            % non-dimensional force increment
k0 = 0;                   % guess for phi'(0), reused step to step
eps = 0; sigma = 0; nu = NaN;

f = 0;
while true
    f = f + sgn * f_step;
    F = [0 0];  F(dir) = f;

    % shooting: find phi'(0) that makes the mid-wall moment-free
    k0 = fzero(@(k) end_slope(F, th, k), k0, optimset('Display', 'off'));
    [dx, dy] = half_wall_end(F, th, k0);

    % cell strains (units of l/2: full wall spans 2c by 2s, h = 2a)
    e1 = (dx - c) / c;
    e2 = (2*dy - 2*s) / (2*a + 2*s);

    % stress: F = f Es b t^3 / (3 l^2), spread over the cell face
    if dir == 1
        sig = f * Es * g^3 / (3 * (a + s));
        e_load = e1;  e_lat = e2;
    else
        sig = f * Es * g^3 / (3 * c);
        e_load = e2;  e_lat = e1;
    end

    % walls cannot stretch: in tension they straighten and strain saturates
    % at 1/cos(theta) - 1 (x1). Stop if the solver jumps to another shape.
    if abs(e_load - eps(end)) > 0.01 || sign(e_load) ~= sgn
        fprintf('  stopped: wall nearly straight, strain limit reached\n');
        break
    end

    eps(end+1)   = e_load;           %#ok<AGROW>
    sigma(end+1) = sig;              %#ok<AGROW>
    nu(end+1)    = -e_lat / e_load;  %#ok<AGROW>

    if abs(e_load) >= abs(eps_end), break; end
end

oh = struct('eps', eps(:), 'sigma', sigma(:), 'nu', nu(:));

fprintf('Ohtaki  theta %g, h/l %.2f, t/l %.3f, load along x%d\n', theta_deg, a, g, dir);
fprintf('  first step:  eps = %+.4f   nu = %.3f   (small-strain, ~G&A)\n', eps(2), nu(2));
fprintf('  last step:   eps = %+.4f   nu = %.3f   sigma = %.4g MPa\n', eps(end), nu(end), sigma(end));
end

% ---- helpers: integrate the elastica along the half wall ------------------
function r = end_slope(F, th, k0)
    [~, Y] = ode45(@(s, y) elastica(y, F), [0 1], [th; k0; 0; 0], odeset('RelTol', 1e-9));
    r = Y(end, 2);                      % phi'(1), should be zero
end

function [dx, dy] = half_wall_end(F, th, k0)
    [~, Y] = ode45(@(s, y) elastica(y, F), [0 1], [th; k0; 0; 0], odeset('RelTol', 1e-9));
    dx = Y(end, 3);  dy = Y(end, 4);    % mid-wall position from the joint
end

function dy = elastica(y, F)
    % y = [phi; phi'; x; y]
    dy = [y(2); F(1)*sin(y(1)) - F(2)*cos(y(1)); cos(y(1)); sin(y(1))];
end
