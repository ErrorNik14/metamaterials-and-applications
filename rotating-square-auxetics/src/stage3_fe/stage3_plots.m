function figs = stage3_plots(P, M, T, Conv, Wcmp)
%STAGE3_PLOTS  Deformed shape, nu comparison vs w, convergence with N.
figs = gobjects(1, 3);
sc = P.fe.plot_scale;

%% Fig 3.1 - deformed array (magnified) over the undeformed outline
figs(1) = new_figure(P, P.plot.fig_width_cm * 0.6, P.plot.fig_width_cm * 0.6);
ax = axes(figs(1)); hold(ax, 'on');
U  = reshape(T.U, 3, []).';
Xd = M.X + sc * U(:, 1:2);
for j = 1:M.N
    for i = 1:M.N
        c = squeeze(M.corner(j, i, :));
        patch(ax, M.X(c,1), M.X(c,2), 'w', 'FaceColor', 'none', ...
            'EdgeColor', P.plot.c_ghost, 'LineStyle', '--', 'LineWidth', 0.8);
        patch(ax, Xd(c,1), Xd(c,2), 'w', 'FaceColor', P.plot.c_square, ...
            'EdgeColor', P.plot.c_series1, 'LineWidth', 0.8);
    end
end
lig = M.conn(M.is_lig, :);
for e = 1:size(lig, 1)
    plot(ax, Xd(lig(e,:),1), Xd(lig(e,:),2), '-', 'Color', P.plot.c_series2, 'LineWidth', 1.5);
end
axis(ax, 'equal'); axis(ax, 'off');
title(ax, sprintf('%dx%d FE array, \\epsilon_x = %+.0f %% (x%g)', M.N, M.N, ...
    100*P.fe.applied_strain, sc), 'FontSize', P.plot.title_size);

%% Fig 3.2 - Poisson's ratio vs ligament width: Stage 1, Stage 2, FE
figs(2) = new_figure(P);
ax = axes(figs(2)); hold(ax, 'on');
yline(ax, -1, '-', 'Color', P.plot.c_muted, 'LineWidth', 1, ...
    'DisplayName', 'Stage 1: rigid squares, pins');
plot(ax, Wcmp.w_fine, Wcmp.nu_S2, 'Color', P.plot.c_series1, ...
    'DisplayName', 'Stage 2: energy method (periodic)');
plot(ax, Wcmp.w, Wcmp.nu_FE_per, 'd', 'Color', P.plot.c_ink, 'MarkerSize', 9, ...
    'LineStyle', 'none', 'LineWidth', 1.2, 'DisplayName', 'Stage 3: FE periodic cell');
plot(ax, Wcmp.w, Wcmp.nu_FE_N, 'o', 'Color', P.plot.c_series2, ...
    'MarkerFaceColor', P.plot.c_series2, 'LineStyle', 'none', ...
    'DisplayName', sprintf('Stage 3: FE %dx%d', P.fe.N_cells, P.fe.N_cells));
plot(ax, Wcmp.w, Wcmp.nu_FE_big, 's', 'Color', P.plot.c_series3, ...
    'MarkerFaceColor', P.plot.c_series3, 'LineStyle', 'none', ...
    'DisplayName', sprintf('Stage 3: FE %dx%d', Wcmp.N_big, Wcmp.N_big));
xline(ax, P.geom.w, ':', 'your design', 'Color', P.plot.c_ink, 'LineWidth', 1, ...
    'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal', ...
    'HandleVisibility', 'off');
xlabel(ax, 'Ligament width w (mm)');
ylabel(ax, 'Poisson''s ratio \nu_{xy} (-)');
title(ax, 'Three independent models of the same lattice');
legend(ax, 'Location', 'southeast');
grid(ax, 'on'); ax.GridAlpha = 0.12;

%% Fig 3.3 - convergence with array size (edge effects)
figs(3) = new_figure(P);
tl = tiledlayout(figs(3), 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, Conv.N, Conv.nu, 'o-', 'Color', P.plot.c_series2, ...
    'MarkerFaceColor', P.plot.c_series2, 'DisplayName', 'FE finite array');
yline(ax, Conv.nu_S2, '--', 'Color', P.plot.c_series1, 'DisplayName', 'Stage 2');
yline(ax, Conv.nu_per, ':', 'Color', P.plot.c_ink, 'LineWidth', 1.5, ...
    'DisplayName', 'FE periodic');
xlabel(ax, 'Array size N (squares per side)');
ylabel(ax, '\nu_{xy} (-)');
grid(ax, 'on'); ax.GridAlpha = 0.12;
ax = nexttile(tl); hold(ax, 'on');
plot(ax, Conv.N, Conv.E, 'o-', 'Color', P.plot.c_series2, ...
    'MarkerFaceColor', P.plot.c_series2, 'DisplayName', 'FE finite array');
yline(ax, Conv.E_S2, '--', 'Color', P.plot.c_series1, 'DisplayName', 'Stage 2');
yline(ax, Conv.E_per, ':', 'Color', P.plot.c_ink, 'LineWidth', 1.5, ...
    'DisplayName', 'FE periodic');
xlabel(ax, 'Array size N (squares per side)');
ylabel(ax, 'E_{eff} (MPa)');
grid(ax, 'on'); ax.GridAlpha = 0.12;
lgd = legend(ax, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';                  % one shared legend under both panels
title(tl, 'Finite arrays vs the infinite lattice (note the tiny y-ranges)', ...
    'FontSize', P.plot.title_size);
end
