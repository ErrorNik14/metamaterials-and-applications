function figs = stage2_plots(P, C0, SW, NL)
%STAGE2_PLOTS  Stage-2 figures: nu vs hinge width, parametric sweeps,
%   nonlinear stress-strain. C0 = properties at the default design.
g = P.geom;
w_min_print = P.print.min_perims * P.print.nozzle;
cols = {P.plot.c_series1, P.plot.c_series2, P.plot.c_series3};
figs = gobjects(1, 3);

%% Fig 2.1 - Poisson's ratio vs ligament width, for several square sizes
figs(1) = new_figure(P);
ax = axes(figs(1)); hold(ax, 'on');
for k = 1:numel(P.stage2.a_family)
    plot(ax, SW.w, SW.nu_family(:, k), 'Color', cols{k}, ...
        'DisplayName', sprintf('a = %g mm', P.stage2.a_family(k)));
end
yline(ax, -1, ':', 'ideal pin hinges (Stage 1)', 'Color', P.plot.c_muted, 'LineWidth', 1, ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', ...
    'HandleVisibility', 'off');
xline(ax, w_min_print, '--', sprintf('min. printable w = %.1f mm', w_min_print), ...
    'Color', P.plot.c_muted, 'LineWidth', 1, 'LabelVerticalAlignment', 'top', ...
    'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
ylim(ax, [-1.03 -0.6]);
plot(ax, g.w, C0.nu, 'o', 'Color', P.plot.c_ink, 'MarkerFaceColor', P.plot.c_ink, ...
    'DisplayName', sprintf('your design: \\nu = %.3f', C0.nu));
xlabel(ax, 'Ligament width w (mm)');
ylabel(ax, 'Poisson''s ratio \nu_{xy} (-)');
title(ax, sprintf('Compliant hinges: \\nu = -1/(1 + 2w^2/3b^2),  \\theta_0 = %g\\circ', ...
    g.theta0_deg));
legend(ax, 'Location', 'west');
grid(ax, 'on'); ax.GridAlpha = 0.12;

%% Fig 2.2 - parametric sweeps: E_eff (top) and nu (bottom)
figs(2) = new_figure(P, P.plot.fig_width_cm, 12);
tl = tiledlayout(figs(2), 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
xs   = {SW.w, SW.Lh, SW.a};
Cs   = {SW.vs_w, SW.vs_Lh, SW.vs_a};
xdef = [g.w, g.Lh, g.a];
xlab = {'Ligament width w (mm)', 'Ligament length L_h (mm)', 'Square side a (mm)'};
for k = 1:3
    ax = nexttile(tl, k); hold(ax, 'on');
    semilogy(ax, xs{k}, Cs{k}.E_eff, 'Color', P.plot.c_series1);
    semilogy(ax, xdef(k), C0.E_eff, 'o', 'Color', P.plot.c_ink, 'MarkerFaceColor', P.plot.c_ink);
    set(ax, 'YScale', 'log'); grid(ax, 'on'); ax.GridAlpha = 0.12;
    if k == 1, ylabel(ax, 'E_{eff} (MPa)'); end
end
for k = 1:3
    ax = nexttile(tl, k + 3); hold(ax, 'on');
    plot(ax, xs{k}, Cs{k}.nu, 'Color', P.plot.c_series2);
    plot(ax, xdef(k), C0.nu, 'o', 'Color', P.plot.c_ink, 'MarkerFaceColor', P.plot.c_ink);
    ylim(ax, [-1.02 -0.6]); grid(ax, 'on'); ax.GridAlpha = 0.12;
    xlabel(ax, xlab{k});
    if k == 1, ylabel(ax, '\nu_{xy} (-)'); end
end
title(tl, sprintf('One parameter at a time  (dot: a = %g, L_h = %g, w = %g mm)', ...
    g.a, g.Lh, g.w), 'FontSize', P.plot.title_size);

%% Fig 2.3 - nonlinear stress-strain along the mechanism path
figs(3) = new_figure(P);
ax = axes(figs(3)); hold(ax, 'on');
plot(ax, 100*NL.eps, NL.sigma, 'Color', P.plot.c_series1, ...
    'DisplayName', 'large rotation (rotation mode only)');
plot(ax, 100*NL.eps, C0.E_eff * NL.eps, '--', 'Color', P.plot.c_series2, ...
    'DisplayName', sprintf('linear, E_{eff} = %.1f MPa', C0.E_eff));
xline(ax, 100*NL.eps_lock, '-', sprintf('open lock %+.1f %%', 100*NL.eps_lock), ...
    'Color', P.plot.c_lock, 'LineWidth', 1, 'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal', ...
    'HandleVisibility', 'off');
xline(ax, 0, '-', 'Color', P.plot.c_muted, 'LineWidth', 0.5, 'HandleVisibility', 'off');
% stress rises without bound at the lock: show the useful range only
ylim(ax, [1.2*min(NL.sigma), 4*C0.E_eff*NL.eps_lock]);
yline(ax, 0, '-', 'Color', P.plot.c_muted, 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel(ax, 'Nominal strain \epsilon_x (%)');
ylabel(ax, 'Nominal stress \sigma_x (MPa)');
title(ax, 'Lattice stiffens towards the open lock (contact/yield: Stage 4)');
legend(ax, 'Location', 'northwest');
grid(ax, 'on'); ax.GridAlpha = 0.12;
end
