function S6 = stage6_run(P, S2, S3, S4)
%STAGE6_RUN  Stage 6: CAD exports, print check, results CSV, Ansys template.
fprintf('\n=== Stage 6: exports ===\n');
g = P.geom;
cad = fullfile(P.paths.results, 'cad');
if ~exist(cad, 'dir'), mkdir(cad); end

% --- CAD -----------------------------------------------------------------------
[ps, parts] = rs_cad_polyshape(P);
dxf = fullfile(cad, 'S6_lattice_outline.dxf');
stl = fullfile(cad, 'S6_lattice.stl');
write_dxf_polyshape(ps, dxf);
TR = write_stl_extrude(ps, g.t, stl);
SW = write_solidworks_files(P, S2, cad);
fprintf('  saved %s\n  saved %s\n  saved %s\n  saved %s\n  saved %s\n', ...
    dxf, stl, SW.csv, SW.eq, SW.squares);

% --- print feasibility -----------------------------------------------------------
PC = rs_print_check(P, ps);

% --- results CSV + Ansys template ----------------------------------------------------
res = write_results_csv(P, S2, S3, S4, P.paths.data);
tpl = fullfile(P.paths.data, 'ansys_results_TEMPLATE.csv');
fid = fopen(tpl, 'w');
fprintf(fid, 'strain_x,nu,stress_MPa\n');      % header only: fill with YOUR Ansys data
fclose(fid);
fprintf('  saved %s\n  saved %s (header only - fill in your Ansys data)\n', res, tpl);

% --- outline figure for the slides ------------------------------------------------------
fh = new_figure(P, P.plot.fig_width_cm*0.6, P.plot.fig_width_cm*0.6);
ax = axes(fh);
plot(ax, ps, 'FaceColor', P.plot.c_square, 'EdgeColor', P.plot.c_series1, 'FaceAlpha', 1);
axis(ax, 'equal'); grid(ax, 'on'); ax.GridAlpha = 0.12;
xlabel(ax, 'x (mm)'); ylabel(ax, 'y (mm)');
title(ax, sprintf('Print outline: %dx%d, a = %g, L_h = %g, w = %g mm', ...
    P.cad.N, P.cad.N, g.a, g.Lh, g.w), 'FontSize', P.plot.title_size);
save_figure(fh, 'S6_print_outline', P);

% --- checks ------------------------------------------------------------------------------
tests = {};
% 1. STL volume (divergence theorem) = outline area x t
V = TR.Points;  F = TR.ConnectivityList;
vol = sum(dot(V(F(:,1),:), cross(V(F(:,2),:), V(F(:,3),:), 2), 2)) / 6;
tests(end+1,:) = {'STL volume vs area*t (rel)', abs(vol - area(ps)*g.t)/(area(ps)*g.t), 1e-9};
% 2. every STL edge is shared by exactly 2 triangles (watertight)
E = sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2);
[~, ~, ic] = unique(E, 'rows');
tests(end+1,:) = {'STL watertight (# bad edges)', sum(accumarray(ic, 1) ~= 2), 0};
% 3. outline area = squares + exposed ligament strips (union removes overlap)
A_sq = numel(parts.squares) * g.a^2;
tests(end+1,:) = {'outline area >= squares area', double(area(ps) < A_sq - 1e-9), 0};
% 4. DXF has one LINE per boundary segment
txt = fileread(dxf);
nseg = 0;
for b = 1:numboundaries(ps), [x, ~] = boundary(ps, b); nseg = nseg + numel(x) - 1; end
tests(end+1,:) = {'DXF LINE count', abs(numel(strfind(txt, sprintf('\nLINE\n'))) - nseg), 0};
% 5. SolidWorks CSV reads back the parameters exactly
Tsw = readtable(SW.csv, 'Delimiter', ',', 'ReadVariableNames', true, ...
    'VariableNamingRule', 'preserve', 'TextType', 'string');
tests(end+1,:) = {'SolidWorks CSV round trip (a, Lh, w, t)', ...
    max(abs(Tsw.value(1:4) - [g.a; g.Lh; g.w; g.t])), 1e-12};
% 6. part is one piece with (N-1)^2 pores
tests(end+1,:) = {'one piece with (N-1)^2 pores', double(~PC.items{end,3}), 0};

fprintf('\n  Stage 6 checkpoint tests\n');
fprintf('  %-40s %12s %10s  %s\n', 'test', 'error', 'tol', 'result');
pass = true;
for k = 1:size(tests, 1)
    ok = tests{k,2} <= tests{k,3};
    pass = pass && ok;
    if ok, r = 'PASS'; else, r = 'FAIL'; end
    fprintf('  %-40s %12.3e %10.1e  %s\n', tests{k,1}, tests{k,2}, tests{k,3}, r);
end
S6.ps = ps;  S6.print = PC;  S6.pass = pass;  S6.files = SW;
if pass
    fprintf('  Stage 6: ALL CHECKS PASSED\n');
else
    warning('Stage 6: at least one check FAILED - see table above.');
end
end
