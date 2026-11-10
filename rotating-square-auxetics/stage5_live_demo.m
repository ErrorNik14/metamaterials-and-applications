%% Rotating-squares auxetic: interactive demo
% Open this file in MATLAB, right-click its tab and choose
% "Open as Live Script". Then click on the number in the line
%     theta_deg = 15;
% and use  Live Editor > Insert > Control > Slider  (min 0, max 43, step 0.5).
% Moving the slider re-runs the section and redraws the lattice.

%% Setup (run once)
here = fileparts(which('stage5_live_demo'));
addpath(here); addpath(genpath(fullfile(here, 'src')));
P = params_default();

%% Choose the square rotation angle
theta_deg = 15;   % <-- put the slider on this value

th0 = deg2rad(P.geom.theta0_deg);
St  = rs_state_at_theta(P, deg2rad(theta_deg));
L   = rs_layout(deg2rad(theta_deg), P.geom.a, P.geom.Lh, th0, P.anim.N, P.anim.lig_points);
L0  = rs_layout(th0, P.geom.a, P.geom.Lh, th0, P.anim.N, P.anim.lig_points);

fh = new_figure(P, 12, 12);
ax = axes(fh);
draw_lattice(ax, P, L, L0, St.rgb);
title(ax, sprintf('\\theta = %.1f\\circ,  \\epsilon = %+.2f %%  [%s]', ...
    theta_deg, 100*St.eps, St.status));
fprintf('strain %+.2f %% | nu (ideal) %.0f | nu (beam hinges) %.4f | hinge stress %s | %s\n', ...
    100*St.eps, St.nu_ideal, St.nu_compliant, St.sig_text, St.status);
