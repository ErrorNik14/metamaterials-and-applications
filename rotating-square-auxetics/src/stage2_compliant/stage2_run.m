function S2 = stage2_run(P)
%STAGE2_RUN  Stage 2: compliant (beam) hinges - stiffness and nu deviation.
fprintf('\n=== Stage 2: compliant hinges (energy method) ===\n');

E   = hinge_modulus(P);
g   = P.geom;
th0 = deg2rad(g.theta0_deg);
rf  = P.model.root_length_factor;

C0 = rs_compliant_props(g.a, g.Lh, g.w, g.t, E, th0, rf);   % closed form
R0 = rs_energy_check   (g.a, g.Lh, g.w, g.t, E, th0, rf);   % numerical cell
SW = stage2_sweeps(P);
NL = rs_nonlinear_response(g.a, g.Lh, g.w, g.t, E, th0, rf, P.stage2.lock_margin);

fprintf('  hinge modulus used          E   = %.0f MPa (%s)\n', E, P.model.hinge_E_mode);
fprintf('  ligament bending stiffness  k_r = %.1f N.mm/rad\n', C0.k_r);
fprintf('  ligament axial stiffness    k_a = %.0f N/mm\n', C0.k_a);
fprintf('  lever arm b = a(cos-sin)        = %.3f mm\n', C0.b);
fprintf('  as-printed cell (2x2)           = %.3f mm square\n', C0.Lcell);
fprintf('  rho = w^2/(12 b^2)              = %.3e\n', C0.rho);
fprintf('  Poisson''s ratio  nu             = %.4f   (numerical cell %.4f)\n', C0.nu, R0.nu);
fprintf('  effective modulus E_eff         = %.2f MPa (numerical cell %.2f)\n', C0.E_eff, R0.E_eff);
fprintf('  E_eff / E_PLA                   = %.4f  (%.0fx softer)\n', C0.E_ratio, 1/C0.E_ratio);
fprintf('  open lock with ligament arms    = %.2f deg, strain %+.2f %%\n', ...
    rad2deg(NL.theta_lock), 100*NL.eps_lock);

pass = stage2_checks(P, C0, R0, SW, NL);

figs  = stage2_plots(P, C0, SW, NL);
names = {'S2_nu_vs_hinge_width', 'S2_parametric_sweeps', 'S2_nonlinear_stress_strain'};
for k = 1:numel(figs)
    save_figure(figs(k), names{k}, P);
end

S2.C0 = C0; S2.numeric = R0; S2.sweeps = SW; S2.NL = NL; S2.pass = pass;
if pass
    fprintf('  Stage 2: ALL CHECKS PASSED\n');
else
    warning('Stage 2: at least one check FAILED - see table above.');
end
end
