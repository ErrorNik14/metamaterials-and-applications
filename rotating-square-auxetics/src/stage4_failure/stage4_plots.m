function figs = stage4_plots(P, S4)
%STAGE4_PLOTS  Hinge stress, safe window, buckling mode, design chart.
SF = P.fail.safety_factor;  sy = S4.sy;  H = S4.H;
figs = gobjects(1, 4);

%% Fig 4.1 - peak hinge stress vs applied strain
figs(1) = new_figure(P);
ax = axes(figs(1)); hold(ax, 'on');
plot(ax, 100*H.NL.eps, H.NL.sig, 'Color', P.plot.c_series1, ...
    'DisplayName', 'total (large rotation)');
plot(ax, 100*H.NL.eps, H.NL.sig_b, '--', 'Color', P.plot.c_series2, ...
    'DisplayName', 'bending part');
plot(ax, 100*H.NL.eps, H.NL.sig_a, ':', 'Color', P.plot.c_series3, ...
    'DisplayName', 'axial part');
yline(ax, sy, '-', sprintf('yield %.0f MPa', sy), 'Color', P.plot.c_lock, ...
    'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(ax, sy/SF, '--', sprintf('yield / SF (%.1f)', SF), 'Color', P.plot.c_lock, ...
    'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
plot(ax, 100*S4.ycross, [sy sy], 'o', 'Color', P.plot.c_ink, ...
    'MarkerFaceColor', P.plot.c_ink, 'HandleVisibility', 'off');
text(ax, 100*S4.ycross(2) + 0.3, sy, sprintf('%+.2f %%', 100*S4.ycross(2)), ...
    'VerticalAlignment', 'top');
text(ax, 100*S4.ycross(1) - 0.3, sy, sprintf('%+.2f %%', 100*S4.ycross(1)), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
xlim(ax, 100 * 4 * [-1 1] * S4.eps_y_lin);
ylim(ax, [0 2*sy]);
xlabel(ax, 'Applied lattice strain \epsilon_x (%)');
ylabel(ax, 'Peak ligament stress (MPa)');
title(ax, 'The ligaments yield at about 1 % lattice strain');
legend(ax, 'Location', 'north');
grid(ax, 'on'); ax.GridAlpha = 0.12;

%% Fig 4.2 - all limits on one strain axis
% Limits span 0.5 % to >100 %, so plot |strain| on a log axis, one row each.
figs(2) = new_figure(P, P.plot.fig_width_cm, 10);
ax = axes(figs(2)); hold(ax, 'on');
names = {'SAFE, tension', 'SAFE, compression', 'first yield, tension', ...
         'first yield, compression', 'open lock (tension)', ...
         'square contact (compression)', 'in-plane buckling (compression)'};
vals  = 100*[S4.safe(2), S4.safe(1), S4.ycross(2), S4.ycross(1), ...
             S4.lock_open, S4.contact.eps, S4.buckle.eps_cr(1)];
nrow = numel(vals);
for k = 1:nrow
    y = nrow + 1 - k;                                     % top row first
    if k <= 2, col = P.plot.c_series3; else, col = P.plot.c_lock; end
    plot(ax, [0.1 abs(vals(k))], [y y], '-', 'Color', col, 'LineWidth', 2);
    plot(ax, abs(vals(k)), y, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 8);
    text(ax, abs(vals(k)) * 1.25, y, sprintf('%+.2f %%', vals(k)), ...
        'VerticalAlignment', 'middle');
end
set(ax, 'XScale', 'log', 'YTick', 1:nrow, 'YTickLabel', flip(names));
xlim(ax, [0.1, 10^ceil(log10(max(abs(vals))) + 0.3)]);
ylim(ax, [0.4 nrow + 0.6]);
xlabel(ax, '|Applied lattice strain| (%, log scale)');
title(ax, sprintf('Limits of the printed part (SF = %.1f)', SF));
grid(ax, 'on'); ax.GridAlpha = 0.12;

%% Fig 4.3 - first in-plane buckling mode
figs(3) = new_figure(P, P.plot.fig_width_cm*0.6, P.plot.fig_width_cm*0.6);
ax = axes(figs(3)); hold(ax, 'on');
M = S4.mesh;  phi = reshape(S4.buckle.modes(:,1), 3, []).';
scale = 0.15 * M.p0 / max(abs(phi(:,1:2)), [], 'all');
Xd = M.X + scale * phi(:,1:2);
for j = 1:M.N
    for i = 1:M.N
        c = squeeze(M.corner(j,i,:));
        patch(ax, M.X(c,1), M.X(c,2), 'w', 'FaceColor', 'none', ...
            'EdgeColor', P.plot.c_ghost, 'LineStyle', '--', 'LineWidth', 0.8);
        patch(ax, Xd(c,1), Xd(c,2), 'w', 'FaceColor', P.plot.c_square, ...
            'EdgeColor', P.plot.c_series1, 'LineWidth', 0.8);
    end
end
axis(ax, 'equal'); axis(ax, 'off');
title(ax, sprintf('Buckling mode 1, \\epsilon_{cr} = %+.1f %%', 100*S4.buckle.eps_cr(1)), ...
    'FontSize', P.plot.title_size);

%% Fig 4.4 - design chart: safe strain vs ligament length
figs(4) = new_figure(P);
ax = axes(figs(4)); hold(ax, 'on');
cols = {P.plot.c_series1, P.plot.c_series2, P.plot.c_series3};
Dc = S4.design;
for k = 1:numel(Dc.w)
    plot(ax, Dc.Lh, 100*Dc.eps_safe(:,k), 'Color', cols{k}, ...
        'DisplayName', sprintf('w = %.1f mm', Dc.w(k)));
end
plot(ax, P.geom.Lh, 100*S4.safe(2), 'o', 'Color', P.plot.c_ink, ...
    'MarkerFaceColor', P.plot.c_ink, ...
    'DisplayName', sprintf('your design (incl. x%.2f edge factor)', S4.edge_factor));
xlabel(ax, 'Ligament length L_h (mm)');
ylabel(ax, 'Safe strain, small-strain estimate (%)');
title(ax, 'Longer, thinner ligaments = larger elastic range');
legend(ax, 'Location', 'northwest');
grid(ax, 'on'); ax.GridAlpha = 0.12;
end
