function S3 = stage3_run(P, S2)
%STAGE3_RUN  Stage 3: own beam-element FE check of the N x N array.
fprintf('\n=== Stage 3: frame finite-element check ===\n');
g  = fe_geometry_from_params(P);
N  = P.fe.N_cells;
ea = P.fe.applied_strain;
th0 = deg2rad(P.geom.theta0_deg);
rf  = P.model.root_length_factor;

% --- default design, tension and compression --------------------------------
M  = rs_fe_mesh(g, N);
T  = rs_fe_uniaxial(M,  ea);
Tc = rs_fe_uniaxial(M, -ea);
fprintf('  mesh: %d nodes, %d elements (%d per ligament)\n', ...
    size(M.X,1), size(M.conn,1), P.fe.n_el_lig);

% --- periodic cell (infinite lattice): like-for-like with Stage 2 ---------------
sig = S2.C0.E_eff * ea;                     % stress giving about eps_applied
Mp  = rs_fe_mesh(g, 2, true);
Pt  = rs_fe_periodic(Mp,  sig);
Pc  = rs_fe_periodic(Mp, -sig);

% --- convergence with array size ------------------------------------------------
Conv.N = P.fe.N_convergence(:);
Conv.nu = zeros(size(Conv.N));  Conv.E = Conv.nu;
for k = 1:numel(Conv.N)
    Rk = rs_fe_uniaxial(rs_fe_mesh(g, Conv.N(k)), ea);
    Conv.nu(k) = Rk.nu;  Conv.E(k) = Rk.E_eff;
end
Conv.nu_S2 = S2.C0.nu;  Conv.E_S2 = S2.C0.E_eff;
Conv.nu_per = Pt.nu;    Conv.E_per = Pt.E_eff;

% --- comparison over ligament width ---------------------------------------------
Wcmp.w = P.fe.w_compare(:);
Wcmp.N_big = P.fe.N_convergence(end);
Wcmp.nu_FE_N = zeros(size(Wcmp.w));  Wcmp.nu_FE_big = Wcmp.nu_FE_N;
Wcmp.E_FE_big = Wcmp.nu_FE_N;  Wcmp.nu_FE_per = Wcmp.nu_FE_N;  Wcmp.E_FE_per = Wcmp.nu_FE_N;
for k = 1:numel(Wcmp.w)
    gk = fe_geometry_from_params(P, Wcmp.w(k));
    Rp = rs_fe_periodic(rs_fe_mesh(gk, 2, true), sig);
    Wcmp.nu_FE_per(k) = Rp.nu;  Wcmp.E_FE_per(k) = Rp.E_eff;
    Wcmp.nu_FE_N(k) = rs_fe_uniaxial(rs_fe_mesh(gk, N), ea).nu;
    Rb = rs_fe_uniaxial(rs_fe_mesh(gk, Wcmp.N_big), ea);
    Wcmp.nu_FE_big(k) = Rb.nu;  Wcmp.E_FE_big(k) = Rb.E_eff;
end
Wcmp.w_fine = linspace(min(Wcmp.w), max(Wcmp.w), 100).';
Cf = rs_compliant_props(P.geom.a, P.geom.Lh, Wcmp.w_fine, P.geom.t, g.E, th0, rf);
Wcmp.nu_S2 = Cf.nu;
Cw = rs_compliant_props(P.geom.a, P.geom.Lh, Wcmp.w, P.geom.t, g.E, th0, rf);
Wcmp.nu_S2_at_w = Cw.nu;  Wcmp.E_S2_at_w = Cw.E_eff;

% --- one comparison table ---------------------------------------------------------
Method = {'Stage 1: rigid + pins (infinite)'; 'Stage 2: energy method (infinite)'; ...
          'Stage 3: FE periodic, tension'; 'Stage 3: FE periodic, compression'; ...
          sprintf('Stage 3: FE %dx%d array, tension', N, N); ...
          sprintf('Stage 3: FE %dx%d array, compression', N, N); ...
          sprintf('Stage 3: FE %dx%d array, tension', Wcmp.N_big, Wcmp.N_big)};
nu    = [-1; S2.C0.nu; Pt.nu; Pc.nu; T.nu; Tc.nu; Conv.nu(end)];
E_eff = [NaN; S2.C0.E_eff; Pt.E_eff; Pc.E_eff; T.E_eff; Tc.E_eff; Conv.E(end)];
dnu_pct = 100 * (nu - S2.C0.nu) / abs(S2.C0.nu);
dE_pct  = 100 * (E_eff - S2.C0.E_eff) / S2.C0.E_eff;
Tab = table(Method, nu, E_eff, dnu_pct, dE_pct, 'VariableNames', ...
    {'Method', 'nu', 'E_eff_MPa', 'nu_diff_vs_S2_pct', 'E_diff_vs_S2_pct'});
fprintf('\n'); disp(Tab);
writetable(Tab, fullfile(P.paths.data, 'S3_comparison_table.csv'));
Wtab = table(Wcmp.w, Wcmp.nu_S2_at_w, Wcmp.nu_FE_per, Wcmp.nu_FE_N, Wcmp.nu_FE_big, ...
    Wcmp.E_S2_at_w, Wcmp.E_FE_per, Wcmp.E_FE_big, 'VariableNames', {'w_mm', 'nu_S2', ...
    'nu_FE_periodic', sprintf('nu_FE_%dx%d', N, N), sprintf('nu_FE_%dx%d', Wcmp.N_big, Wcmp.N_big), ...
    'E_S2_MPa', 'E_FE_periodic_MPa', sprintf('E_FE_%dx%d_MPa', Wcmp.N_big, Wcmp.N_big)});
disp(Wtab);
writetable(Wtab, fullfile(P.paths.data, 'S3_nu_vs_width.csv'));

pass = stage3_checks(P, T, Tc, Pt, Pc, Wcmp, S2.C0);

figs  = stage3_plots(P, M, T, Conv, Wcmp);
names = {'S3_deformed_array', 'S3_nu_comparison', 'S3_convergence'};
for k = 1:numel(figs)
    save_figure(figs(k), names{k}, P);
end

S3.mesh = M; S3.tension = T; S3.compression = Tc; S3.conv = Conv;
S3.periodic = Pt;
S3.width = Wcmp; S3.table = Tab; S3.pass = pass;
if pass
    fprintf('  Stage 3: ALL CHECKS PASSED\n');
else
    warning('Stage 3: at least one check FAILED - see table above.');
end
end
