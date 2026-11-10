function S1 = stage1_run(P)
%STAGE1_RUN  Stage 1: ideal rotating-squares kinematics (rigid squares, pins).
%   S1 = stage1_run(P) computes the closed-form model, checks it against an
%   independently built array, prints a summary, saves 4 figures.
fprintf('\n=== Stage 1: ideal kinematics (rigid squares, pin hinges) ===\n');

a      = P.geom.a;
theta0 = deg2rad(P.geom.theta0_deg);
theta  = deg2rad(linspace(P.stage1.theta_range_deg(1), ...
                          P.stage1.theta_range_deg(2), P.stage1.n_theta));

K  = rs_cell_kinematics(theta, a, theta0);                    % formula
Nm = rs_numeric_poisson(theta, a, P.stage1.n_check, theta0);  % geometry

% --- Summary numbers ---------------------------------------------------------
K0 = rs_cell_kinematics(theta0, a, theta0);
fprintf('  a = %g mm, theta0 = %g deg\n', a, P.geom.theta0_deg);
fprintf('  as-printed cell (2x2 squares): %.3f x %.3f mm, porosity %.1f %%\n', ...
    K0.Lx, K0.Ly, 100*K0.porosity);
fprintf('  open lock    theta = 45 deg  -> max tensile strain     %+.2f %%\n', ...
    100*K.eps_lock_tension);
fprintf('  closed lock  theta =  0 deg  -> max compressive strain %+.2f %%\n', ...
    100*K.eps_lock_compression);

pass = stage1_checks(P, K, Nm);

% --- Figures -------------------------------------------------------------------
figs  = stage1_plots(P, K, Nm);
names = {'S1_geometry_snapshots', 'S1_cell_size_porosity', ...
         'S1_strain_vs_theta', 'S1_poisson_ratio'};
for k = 1:numel(figs)
    save_figure(figs(k), names{k}, P);
end

S1.K = K; S1.numeric = Nm; S1.pass = pass;
if pass
    fprintf('  Stage 1: ALL CHECKS PASSED\n');
else
    warning('Stage 1: at least one check FAILED - see table above.');
end
end
