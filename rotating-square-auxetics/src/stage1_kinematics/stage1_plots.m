function figs = stage1_plots(P, K, Nm)
%STAGE1_PLOTS  The four Stage-1 figures. K = closed form, Nm = numerical check.
a      = P.geom.a;
theta0 = deg2rad(P.geom.theta0_deg);
thd    = rad2deg(K.theta);
lockd  = rad2deg(K.theta_lock_open);
figs   = gobjects(1, 4);

%% Fig 1.1 - geometry snapshots, undeformed (as-printed) outline overlaid
snap = P.stage1.snapshot_deg;
n    = P.stage1.n_vis;
figs(1) = new_figure(P, P.plot.fig_width_cm, P.plot.fig_width_cm / numel(snap) + 2);
tl = tiledlayout(figs(1), 1, numel(snap), 'TileSpacing', 'compact', 'Padding', 'compact');
G0 = rs_build_array(theta0, a, n, n);
shift0 = -mean(reshape(G0.centre, [], 2), 1);
% common limits: the largest array is the fully open one (45 deg)
Gmax = rs_build_array(K.theta_lock_open, a, n, n);
half = max(abs(cell2mat(Gmax.corners(:)) - mean(reshape(Gmax.centre, [], 2), 1)), [], 'all') * 1.05;
for s = 1:numel(snap)
    ax = nexttile(tl);
    G  = rs_build_array(deg2rad(snap(s)), a, n, n);
    shift = -mean(reshape(G.centre, [], 2), 1);
    draw_squares(ax, G, shift, P.plot.c_square, P.plot.c_series1, '-', true);
    draw_squares(ax, G0, shift0, 'none', P.plot.c_ghost, '--');
    axis(ax, 'equal'); axis(ax, 'off');
    xlim(ax, [-half half]); ylim(ax, [-half half]);
    title(ax, sprintf('\\theta = %g\\circ', snap(s)), 'FontSize', P.plot.title_size);
end
title(tl, sprintf(['Rotating squares, a = %g mm  (dashed: as-printed, ' ...
    '\\theta_0 = %g\\circ)'], a, P.geom.theta0_deg), 'FontSize', P.plot.title_size);

%% Fig 1.2 - cell size and area vs theta
figs(2) = new_figure(P);
tl = tiledlayout(figs(2), 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, thd, K.Lx / (2*a), 'Color', P.plot.c_series1);
mark_locks(ax, P, lockd);
xlabel(ax, 'Square rotation \theta (deg)');
ylabel(ax, 'Cell size  L / 2a  (-)');
title(ax, 'Unit-cell size (x and y identical)');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, thd, 100 * K.porosity, 'Color', P.plot.c_series1);
mark_locks(ax, P, lockd);
xlabel(ax, 'Square rotation \theta (deg)');
ylabel(ax, 'Porosity (%)');
title(ax, 'Open area of the lattice');

%% Fig 1.3 - strain vs theta (relative to as-printed state)
figs(3) = new_figure(P);
ax = axes(figs(3)); hold(ax, 'on');
plot(ax, thd, 100 * K.eps_x, 'Color', P.plot.c_series1, 'DisplayName', '\epsilon_x');
plot(ax, thd, 100 * K.eps_y, '--', 'Color', P.plot.c_series2, 'DisplayName', '\epsilon_y');
yline(ax, 0, '-', 'Color', P.plot.c_muted, 'LineWidth', 0.5, 'HandleVisibility', 'off');
mark_locks(ax, P, lockd);
plot(ax, rad2deg(theta0), 0, 'o', 'Color', P.plot.c_ink, ...
    'MarkerFaceColor', P.plot.c_ink, 'HandleVisibility', 'off');
text(ax, rad2deg(theta0) + 1, -2, 'as printed', 'Color', P.plot.c_muted);
text(ax, lockd - 1, 100*K.eps_lock_tension + 2, ...
    sprintf('tension lock %+.1f %%', 100*K.eps_lock_tension), ...
    'Color', P.plot.c_ink, 'HorizontalAlignment', 'right');
text(ax, 3, 100*K.eps_lock_compression, ...
    sprintf('compression lock %+.1f %%', 100*K.eps_lock_compression), ...
    'Color', P.plot.c_ink, 'VerticalAlignment', 'middle');
ylim(ax, [floor(100*K.eps_lock_compression/5)*5, ceil(100*K.eps_lock_tension/5)*5 + 5]);
xlabel(ax, 'Square rotation \theta (deg)');
ylabel(ax, 'Engineering strain (%)');
title(ax, sprintf('Strain relative to \\theta_0 = %g\\circ', P.geom.theta0_deg));
legend(ax, 'Location', 'southeast');

%% Fig 1.4 - Poisson's ratio: closed form vs numerical array
figs(4) = new_figure(P);
ax = axes(figs(4)); hold(ax, 'on');
plot(ax, thd, K.nu_incremental, 'Color', P.plot.c_series1, ...
    'DisplayName', 'Closed form (incremental)');
step = max(1, round(numel(Nm.theta) / 25));     % thin out markers
plot(ax, rad2deg(Nm.theta(1:step:end)), Nm.nu_incremental(1:step:end), 'o', ...
    'Color', P.plot.c_series3, 'MarkerFaceColor', P.plot.c_series3, 'LineStyle', 'none', ...
    'DisplayName', sprintf('Numerical, %dx%d array', P.stage1.n_check, P.stage1.n_check));
yline(ax, P.mat.nu, ':', sprintf('bulk %s, \\nu = %.2f', P.mat.name, P.mat.nu), ...
    'Color', P.plot.c_muted, 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');
mark_locks(ax, P, lockd);
ylim(ax, [-2 1]);
xlabel(ax, 'Square rotation \theta (deg)');
ylabel(ax, 'Poisson''s ratio \nu_{xy} (-)');
title(ax, 'Ideal rotating squares: \nu = -1 at every angle');
legend(ax, 'Location', 'south');
end

function mark_locks(ax, P, lockd)
% Vertical lines at the two kinematic locks.
xline(ax, 0, '-', 'closed', 'Color', P.plot.c_lock, 'LineWidth', 1, ...
    'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal', ...
    'HandleVisibility', 'off');
xline(ax, lockd, '-', 'open lock', 'Color', P.plot.c_lock, 'LineWidth', 1, ...
    'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal', ...
    'HandleVisibility', 'off');
xlim(ax, P.stage1.theta_range_deg);
grid(ax, 'on'); ax.GridAlpha = 0.12;
end
