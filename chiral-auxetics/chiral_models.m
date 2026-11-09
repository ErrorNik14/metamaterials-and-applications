%% CHIRAL HONEYCOMB — three closed-form models from the literature
%
% Implements, for the hexachiral (tangent-ligament) honeycomb of
% Prall and Lakes, the three closed-form models compared in Table 2 of
% the project write-up:
%
%   1) Gibson & Ashby cell-wall-bending method, bending-only limit of
%      this derivation (Eq. 14). nu = -1 exactly.
%   2) Spadoni & Ruzzene (2012), rigid-ring closed form with both axial
%      and bending ligament deformation (their Eqs. 45-46, Eq. 20/21
%      of the note).
%   3) Prall & Lakes (1997), the original paper's own prefactor
%      (Eq. 23 of the note) — about 1.7x higher than (1) and (2) at
%      R/r = 5; this derivation and the literature check in Section 8
%      of the note settle that in favour of (1)/(2).
%
% Units: E is reported as the dimensionless ratio E/Es (set Es = 1
% throughout); multiply by your measured Es to get the physical
% modulus. r, R, t are in consistent length units (mm).
%
% This file is a function file (not a plain script) so that its local
% helper functions below are guaranteed to work the same way in any
% MATLAB version.
%
% Usage:
%   chiral_models()   % always prompts for r, R, t at the command line

function chiral_models()

clc;

%% ---- geometry: always prompt for it ----
r0 = input('Node radius r (mm) [default 5]: ');
if isempty(r0), r0 = 5; end
R0 = input('Node spacing R (mm) [default 25]: ');
if isempty(R0), R0 = 25; end
t0 = input('Ligament thickness t (mm) [default 1]: ');
if isempty(t0), t0 = 1; end

if R0 <= 2*r0
    error('Geometry is not valid: need R > 2*r (so that L = sqrt(R^2-4r^2) is real).');
end

fprintf('\nBaseline r=%.3g R=%.3g t=%.3g\n', r0, R0, t0);
[E14, nu14] = eq14(r0, R0, t0);
[E23, nu23] = eq23_prallLakes(r0, R0, t0);
[Esr, nusr] = eq20_spadoniRuzzene(r0, R0, t0);

fprintf('  Gibson & Ashby (bending-only)   E/Es = %.6e   nu = %.4f\n', E14, nu14);
fprintf('  Prall & Lakes (1997)            E/Es = %.6e   nu = %.4f\n', E23, nu23);
fprintf('  Spadoni & Ruzzene (2012)        E/Es = %.6e   nu = %.4f\n\n', Esr, nusr);

%% ---- Table 2 style comparison: only the geometry you entered ----
fprintf('%-14s %14s %16s %14s\n', 'r,R,t', 'Gibson-Ashby', 'Spadoni-Ruzzene', 'Prall-Lakes');
fprintf('%-14s %14.4e %16.4e %14.4e\n', ...
    sprintf('%g,%g,%g',r0,R0,t0), E14, Esr, E23);

%% ---- sweep t and plot all three E/Es curves (r, R fixed at input) ----
r = r0; Rfix = R0;
tvals = linspace(max(0.05, 0.3*t0), 3.0*t0, 60);
E14v = zeros(size(tvals)); Esrv = E14v; E23v = E14v;
nu14v = E14v; nusrv = E14v; nu23v = E14v;
for k = 1:numel(tvals)
    [E14v(k), nu14v(k)] = eq14(r, Rfix, tvals(k));
    [Esrv(k), nusrv(k)] = eq20_spadoniRuzzene(r, Rfix, tvals(k));
    [E23v(k), nu23v(k)] = eq23_prallLakes(r, Rfix, tvals(k));
end

figure('Name','E/Es vs t');
plot(tvals, E14v, '--', 'LineWidth', 1.5); hold on;
plot(tvals, Esrv, '-',  'LineWidth', 1.5);
plot(tvals, E23v, ':',  'LineWidth', 1.5);
xlabel('ligament thickness t (mm)'); ylabel('E / E_s');
legend('Gibson & Ashby (bending-only)','Spadoni & Ruzzene (2012)','Prall & Lakes (1997)','Location','northwest');
title(sprintf('r = %g mm, R = %g mm', r, Rfix)); grid on;
% --- overlay your own Ansys / physical-test points here, e.g.:
% ansys_t  = [1.0];
% ansys_E  = [ ... your measured E/Es values ... ];
% plot(ansys_t, ansys_E, 'ks', 'MarkerFaceColor','k', 'DisplayName','Ansys');

figure('Name','nu vs t');
plot(tvals, nu14v, '--', 'LineWidth', 1.5); hold on;
plot(tvals, nusrv, '-',  'LineWidth', 1.5);
plot(tvals, nu23v, ':',  'LineWidth', 1.5);
xlabel('ligament thickness t (mm)'); ylabel('\nu');
legend('Gibson & Ashby (bending-only)','Spadoni & Ruzzene (2012)','Prall & Lakes (1997)','Location','southeast');
title(sprintf('r = %g mm, R = %g mm', r, Rfix)); grid on;

end % chiral_models

%% ======================= local functions =======================

function [E, nu] = eq14(r, R, t)
% Gibson & Ashby cell-wall-bending method, bending-only closed form
% (this derivation), Eq. (14). nu is exactly -1.
    L = sqrt(R^2 - 4*r^2);
    E = (sqrt(3)/2) * (t/L)^3 * (R/r)^2;
    nu = -1;
end

function [E, nu] = eq23_prallLakes(r, R, t)
% Original prefactor of Prall & Lakes (1997), Eq. (23). Same t^3 scaling
% as Eq. (14) but about 1.7x higher at R/r = 5 (see Table 1 of the note).
    L = sqrt(R^2 - 4*r^2);
    E = sqrt(3) * (t/L)^3 * (L/r)^2;
    nu = -1;
end

function [E, nu] = eq20_spadoniRuzzene(r, R, t)
% Spadoni & Ruzzene (2012), rigid-ring closed form with axial + bending
% ligament deformation, their Eqs. 45-46 (Eq. (20)/(21) of the note).
    L = sqrt(R^2 - 4*r^2);
    x = t / L;
    c2 = (L^2) / (R^2);
    s2 = 1 - c2;
    E = (4*sqrt(3)*(1 + x^2)) / (2*x^4*c2 + 2*s2 + 6*x^2) * x^3;
    nu = (4*x^2) / (x^4*c2 + s2 + 3*x^2) - 1;
end
