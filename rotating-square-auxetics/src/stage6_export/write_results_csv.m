function file = write_results_csv(P, S2, S3, S4, folder)
%WRITE_RESULTS_CSV  MATLAB predictions vs strain, ready to overlay on Ansys.
%   Columns: strain_x, theta_deg, nu_ideal, nu_beam_theory, nu_FE_NxN,
%   sigma_lattice_MPa (large rotation), sigma_linear_MPa (E_eff*eps),
%   hinge_stress_MPa (elastic estimate, incl. edge factor), status
eps_g = P.export.strain_grid(:);
NL = S2.NL;
% invert the (monotonic) strain-angle path of Stage 2
inside = eps_g >= min(NL.eps) & eps_g <= max(NL.eps);
th = nan(size(eps_g));
th(inside) = interp1(NL.eps, NL.theta, eps_g(inside));
n = numel(eps_g);
sig = nan(n,1); sh = nan(n,1); status = repmat({'out of range'}, n, 1);
for k = find(inside).'
    St = rs_state_at_theta(P, th(k), S4.edge_factor);
    sig(k) = St.sigma;  sh(k) = St.sig_hinge;  status{k} = St.status;
end
T = table(eps_g, rad2deg(th), -ones(n,1), S2.C0.nu*ones(n,1), ...
    S3.tension.nu*ones(n,1), sig, S2.C0.E_eff*eps_g, sh, status, ...
    'VariableNames', {'strain_x', 'theta_deg', 'nu_ideal', 'nu_beam_theory', ...
    sprintf('nu_FE_%dx%d', P.fe.N_cells, P.fe.N_cells), 'sigma_lattice_MPa', ...
    'sigma_linear_MPa', 'hinge_stress_MPa', 'status'});
file = fullfile(folder, 'S6_matlab_results.csv');
writetable(T, file);
end
