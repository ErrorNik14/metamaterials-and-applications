%EXPORT_FE_FOR_REPORT  Write the 5x5 FE solution (1 % tension) to CSV for the report figures.
%   Run from the Rotating_Squares_MATLAB folder:  run report/export_fe_for_report.m
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);  addpath(genpath(fullfile(root, 'src')));
P = params_default();
M = rs_fe_mesh(fe_geometry_from_params(P), P.fe.N_cells);
R = rs_fe_uniaxial(M, P.fe.applied_strain);
U = reshape(R.U, 3, []).';
out = fullfile(root, 'report', 'data');
if ~exist(out, 'dir'), mkdir(out); end
writematrix([M.X U], fullfile(out, 'fe_nodes.csv'));              % x y ux uy rot
writematrix([M.conn double(M.is_lig)], fullfile(out, 'fe_elems.csv'));
C = reshape(permute(M.corner, [3 1 2]), 4, []).';                 % one row per square
writematrix(C, fullfile(out, 'fe_squares.csv'));
fprintf('FE export done: eps_x = %.5f, nu = %.5f\n', R.eps_x, R.nu);
