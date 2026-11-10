function S4 = stage4_run(P, S2)
%STAGE4_RUN  Stage 4: hinge stress, first yield, locking, buckling, safe range.
fprintf('\n=== Stage 4: failure analysis ===\n');
g   = P.geom;
E   = hinge_modulus(P);
sy  = P.mat.sigma_y * P.mat.print_knockdown;      % yield strength used
SF  = P.fail.safety_factor;
Kt  = P.fail.Kt;
th0 = deg2rad(g.theta0_deg);
rf  = P.model.root_length_factor;
NL  = S2.NL;

% --- 1. hinge stress vs applied strain ---------------------------------------
eps_lin = linspace(min(NL.eps), max(NL.eps), 400).';
H = rs_hinge_stress(g.a, g.Lh, g.w, g.t, E, th0, rf, Kt, eps_lin, NL);
eps_y_lin = sy / H.slope;                          % small-strain first yield
ycross    = first_crossing(H.NL.eps, H.NL.sig, sy);       % large rotation
ycross_SF = first_crossing(H.NL.eps, H.NL.sig, sy / SF);

% --- 2. locking ------------------------------------------------------------------
lock_open = NL.eps_lock;                           % dp/dtheta = 0 (Stage 2)
Ct = rs_contact_limit(g.a, g.Lh, th0);             % squares touch in compression

% --- 3. in-plane buckling of the printed array (FE eigenvalue) ---------------------
M = rs_fe_mesh(fe_geometry_from_params(P), P.fe.N_cells);
B = rs_fe_buckling(M, P.fail.buckle_ref_strain, P.fail.n_modes);

% Euler check of one ligament (sway-permitted, clamped ends: K = 1)
I_l  = g.t * g.w^3 / 12;
Pcr_lig = pi^2 * E * I_l / (g.Lh + rf*g.w)^2;
N_lig_at_yield = abs(S2.C0.E_eff * eps_y_lin) * S2.C0.p0 * g.t;

% --- edge amplification in the printed N x N array (FE) ------------------------------
Rp = rs_fe_uniaxial(M, P.fe.applied_strain);
[s_any, ~] = fe_ligament_stress(M, Rp.U, Kt);
H1 = rs_hinge_stress(g.a, g.Lh, g.w, g.t, E, th0, rf, Kt, P.fe.applied_strain);
edge_factor = max(1, s_any / H1.sig);   % >1: an edge ligament is worse than theory

% --- 4. safe operating range --------------------------------------------------------
% yield limits are divided by the edge factor (stress is ~linear in strain here)
lim_t = [ycross_SF(2) / edge_factor, lock_open / SF];
lim_c = [ycross_SF(1) / edge_factor, Ct.eps / SF, B.eps_cr(1) / SF];
safe = [max(lim_c, [], 'omitnan'), min(lim_t, [], 'omitnan')];

fprintf('  yield strength used            %.1f MPa,  SF = %.2f,  Kt = %.2f\n', sy, SF, Kt);
fprintf('  hinge stress per strain        %.0f MPa per unit strain (linear)\n', H.slope);
fprintf('     of which bending            %.0f %%\n', 100*H.sig_b(end)/H.sig(end));
fprintf('  first yield, small strain      +/- %.3f %%\n', 100*eps_y_lin);
fprintf('  first yield, large rotation    %+.3f %% / %+.3f %%\n', 100*ycross);
fprintf('  open lock (tension)            %+.2f %%\n', 100*lock_open);
fprintf('  square contact (compression)   %+.2f %%  (theta = %.1f deg)\n', ...
    100*Ct.eps, rad2deg(Ct.theta));
fprintf('  in-plane buckling (FE %dx%d)    %+.2f %%  (lambda = %.1f)\n', ...
    M.N, M.N, 100*B.eps_cr(1), B.lambda(1));
fprintf('  ligament Euler load %.0f N vs axial load at yield %.1f N (ratio %.0f)\n', ...
    Pcr_lig, N_lig_at_yield, Pcr_lig / N_lig_at_yield);
fprintf('  edge amplification, %dx%d FE   x%.3f (worst ligament / theory)\n', ...
    M.N, M.N, edge_factor);
fprintf('  >> SAFE OPERATING STRAIN       %+.2f %%  to  %+.2f %%\n', 100*safe);

% --- design chart: safe strain vs ligament length, for printable widths ---------------
Lh_v = linspace(P.fail.Lh_design(1), P.fail.Lh_design(2), 60).';
Dc.Lh = Lh_v;  Dc.w = P.fail.w_design;
Dc.eps_safe = zeros(numel(Lh_v), numel(Dc.w));
for k = 1:numel(Dc.w)
    for m = 1:numel(Lh_v)
        Hk = rs_hinge_stress(g.a, Lh_v(m), Dc.w(k), g.t, E, th0, rf, Kt, 0.01);
        Dc.eps_safe(m, k) = sy / SF / Hk.slope;
    end
end

S4.H = H; S4.eps_y_lin = eps_y_lin; S4.ycross = ycross; S4.ycross_SF = ycross_SF;
S4.lock_open = lock_open; S4.contact = Ct; S4.buckle = B; S4.mesh = M;
S4.safe = safe; S4.sy = sy; S4.edge_factor = edge_factor; S4.design = Dc; S4.Pcr_lig = Pcr_lig;

S4.pass = stage4_checks(P, S4, S2);

figs  = stage4_plots(P, S4);
names = {'S4_hinge_stress', 'S4_safe_window', 'S4_buckling_mode', 'S4_design_chart'};
for k = 1:numel(figs)
    save_figure(figs(k), names{k}, P);
end

T = table({'first yield (tension)'; 'first yield (compression)'; 'open lock'; ...
           'square contact'; 'in-plane buckling'; 'SAFE min'; 'SAFE max'}, ...
          100*[ycross(2); ycross(1); lock_open; Ct.eps; B.eps_cr(1); safe(1); safe(2)], ...
          'VariableNames', {'Limit', 'strain_pct'});
writetable(T, fullfile(P.paths.data, 'S4_limits.csv'));
if S4.pass
    fprintf('  Stage 4: ALL CHECKS PASSED\n');
else
    warning('Stage 4: at least one check FAILED - see table above.');
end
end
