function CMP = compare_with_ansys(ansys_csv)
%COMPARE_WITH_ANSYS  Overlay your Ansys results on the MATLAB predictions.
%
%   compare_with_ansys                       reads results/data/ansys_results.csv
%   compare_with_ansys('my_file.csv')        reads another file
%
%   HOW TO FILL THE FILE  (copy results/data/ansys_results_TEMPLATE.csv to
%   results/data/ansys_results.csv and add one row per Ansys load step):
%     strain_x   applied nominal strain = (edge displacement)/(gauge length)
%                gauge length = distance between the centres of the left
%                and right square columns, (N-1)*pitch
%     nu         -eps_y/eps_x, with eps_y from the average y-displacement
%                of the top and bottom rows of square CENTRES (use interior
%                rows if you can, like Stage 3)
%     stress_MPa reaction force on the gripped edge / (N*pitch * t)
%   Leave a cell empty if you did not extract that quantity.
%   Use the SAME boundary conditions as Stage 3 (edge centres move in x,
%   free in y and rotation) so the comparison is like-for-like.
%
%   Output: figure S6_ansys_comparison (+ printed differences), and a
%   calibrated root_length_factor that makes beam theory match Ansys E.
here = fileparts(mfilename('fullpath'));
addpath(here);  addpath(genpath(fullfile(here, 'src')));
P = params_default();
if nargin < 1
    ansys_csv = fullfile(P.paths.data, 'ansys_results.csv');
end
if ~isfile(ansys_csv)
    error(['No Ansys file found at\n  %s\nCopy ansys_results_TEMPLATE.csv to ' ...
           'that name and fill it in (see help compare_with_ansys).'], ansys_csv);
end
A = readtable(ansys_csv);
for v = {'strain_x', 'nu', 'stress_MPa'}
    if ~ismember(v{1}, A.Properties.VariableNames)
        A.(v{1}) = nan(height(A), 1);
    end
end

M = readtable(fullfile(P.paths.data, 'S6_matlab_results.csv'));
nuFEcol = M.Properties.VariableNames{5};
E   = hinge_modulus(P);  g = P.geom;  th0 = deg2rad(g.theta0_deg);
C0  = rs_compliant_props(g.a, g.Lh, g.w, g.t, E, th0, P.model.root_length_factor);

% --- Ansys modulus: least-squares slope through the origin, small strains only
ok = ~isnan(A.stress_MPa) & abs(A.strain_x) <= 0.01 & A.strain_x ~= 0;
CMP.E_ansys = NaN;  CMP.root_factor = NaN;
if any(ok)
    CMP.E_ansys = (A.strain_x(ok).' * A.stress_MPa(ok)) / (A.strain_x(ok).' * A.strain_x(ok));
    % calibrate the effective ligament length so beam theory gives E_ansys
    f = @(rf) rs_compliant_props(g.a, g.Lh, g.w, g.t, E, th0, rf).E_eff - CMP.E_ansys;
    try
        CMP.root_factor = fzero(f, [-0.9*g.Lh/g.w, 20]);
    catch
        CMP.root_factor = NaN;   % Ansys outside what any ligament length gives
    end
end
nu_ok = ~isnan(A.nu);
CMP.nu_ansys = mean(A.nu(nu_ok));

fprintf('\n  MATLAB vs Ansys\n');
fprintf('  %-28s %12s %12s %10s\n', 'quantity', 'MATLAB', 'Ansys', 'diff %');
fprintf('  %-28s %12.4f %12.4f %10.2f\n', 'nu (beam theory)', C0.nu, CMP.nu_ansys, ...
    100*(CMP.nu_ansys - C0.nu)/abs(C0.nu));
fprintf('  %-28s %12.2f %12.2f %10.2f\n', 'E_eff (MPa)', C0.E_eff, CMP.E_ansys, ...
    100*(CMP.E_ansys - C0.E_eff)/C0.E_eff);
if ~isnan(CMP.root_factor)
    fprintf(['  -> setting P.model.root_length_factor = %.3f makes beam theory match\n' ...
             '     Ansys stiffness (effective ligament length %.3f mm instead of %.3f mm)\n'], ...
        CMP.root_factor, g.Lh + CMP.root_factor*g.w, g.Lh);
end

% --- figure -------------------------------------------------------------------
make_output_folders(P);
fh = new_figure(P);
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, 100*M.strain_x, M.sigma_lattice_MPa, 'Color', P.plot.c_series1, ...
    'DisplayName', 'MATLAB large rotation');
plot(ax, 100*M.strain_x, M.sigma_linear_MPa, '--', 'Color', P.plot.c_series1, ...
    'DisplayName', 'MATLAB linear');
plot(ax, 100*A.strain_x, A.stress_MPa, 'o', 'Color', P.plot.c_series2, ...
    'MarkerFaceColor', P.plot.c_series2, 'LineStyle', 'none', 'DisplayName', 'Ansys');
xlabel(ax, 'Strain \epsilon_x (%)'); ylabel(ax, 'Nominal stress (MPa)');
legend(ax, 'Location', 'northwest'); grid(ax, 'on'); ax.GridAlpha = 0.12;
ax = nexttile(tl); hold(ax, 'on');
yline(ax, -1, '-', 'Color', P.plot.c_muted, 'DisplayName', 'ideal (Stage 1)');
plot(ax, 100*M.strain_x, M.nu_beam_theory, 'Color', P.plot.c_series1, ...
    'DisplayName', 'beam theory (Stage 2)');
plot(ax, 100*M.strain_x, M.(nuFEcol), ':', 'Color', P.plot.c_series3, ...
    'DisplayName', 'own FE (Stage 3)');
plot(ax, 100*A.strain_x, A.nu, 'o', 'Color', P.plot.c_series2, ...
    'MarkerFaceColor', P.plot.c_series2, 'LineStyle', 'none', 'DisplayName', 'Ansys');
xlabel(ax, 'Strain \epsilon_x (%)'); ylabel(ax, '\nu_{xy} (-)');
legend(ax, 'Location', 'best'); grid(ax, 'on'); ax.GridAlpha = 0.12;
title(tl, 'MATLAB predictions vs Ansys', 'FontSize', P.plot.title_size);
save_figure(fh, 'S6_ansys_comparison', P);
CMP.table = A;
end
