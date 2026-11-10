function A = stage5_animation(P, S2, S4)
%STAGE5_ANIMATION  GIF + MP4 of one open/close cycle of the printed array.
%   Path: theta0 -> near the open lock -> fully closed -> back to theta0,
%   with a smooth (sinusoidal) speed profile. The dashed outline is the
%   as-printed shape; ligament colour = stress status (green/amber/red).
g   = P.geom;
th0 = deg2rad(g.theta0_deg);
th_hi = th0 + P.anim.open_frac * (S2.NL.theta_lock - th0);
th_lo = deg2rad(P.anim.close_deg);
nf  = P.anim.n_frames;
s   = (0:nf-1).' / nf;
wave = sin(2*pi*s);
theta = th0 + (wave > 0) .* wave * (th_hi - th0) + (wave < 0) .* wave * (th0 - th_lo);

N  = P.anim.N;
L0 = rs_layout(th0, g.a, g.Lh, th0, N, P.anim.lig_points);
Lmax = rs_layout(th_hi, g.a, g.Lh, th0, N, P.anim.lig_points);
half = 0.5 * (N - 1) * Lmax.p + g.a;               % fixed axis limits

gif_path = fullfile(P.paths.anim, 'S5_rotating_squares.gif');
mp4_path = fullfile(P.paths.anim, 'S5_rotating_squares.mp4');
if exist(gif_path, 'file'), delete(gif_path); end
vw = VideoWriter(mp4_path, 'MPEG-4');
vw.FrameRate = P.anim.fps;  vw.Quality = 95;
open(vw);

fh = new_figure(P, P.anim.size_cm(1), P.anim.size_cm(2));
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axL = nexttile(tl);  axR = nexttile(tl);
eps_path = arrayfun(@(t) rs_state_at_theta(P, t).eps, linspace(th_lo, th_hi, 200));
img_size = [];
for k = 1:nf
    St = rs_state_at_theta(P, theta(k), S4.edge_factor);
    L  = rs_layout(theta(k), g.a, g.Lh, th0, N, P.anim.lig_points);
    draw_lattice(axL, P, L, L0, St.rgb);
    xlim(axL, [-half half]);  ylim(axL, [-half half]);
    title(axL, sprintf('\\theta = %5.1f\\circ    \\epsilon = %+5.1f %%', ...
        St.theta_deg, 100*St.eps), 'FontSize', P.plot.title_size);

    cla(axR); hold(axR, 'on');
    plot(axR, rad2deg(linspace(th_lo, th_hi, 200)), 100*eps_path, ...
        'Color', P.plot.c_series1);
    patch(axR, [0 90 90 0], 100*[S4.safe(1) S4.safe(1) S4.safe(2) S4.safe(2)], ...
        P.plot.c_good, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(axR, St.theta_deg, 100*St.eps, 'o', 'MarkerSize', 9, ...
        'Color', P.plot.c_ink, 'MarkerFaceColor', St.rgb);
    xlim(axR, rad2deg([th_lo th_hi])); ylim(axR, 100*[min(eps_path) max(eps_path)]*1.1);
    xlabel(axR, '\theta (deg)');  ylabel(axR, 'Strain \epsilon_x = \epsilon_y (%)');
    title(axR, sprintf('hinge: %s', St.sig_text), 'FontSize', P.plot.title_size - 1);
    grid(axR, 'on'); axR.GridAlpha = 0.12;
    text(axR, rad2deg(th_lo) + 1, 100*S4.safe(2), 'safe band', ...
        'Color', P.plot.c_muted, 'VerticalAlignment', 'bottom');

    drawnow;
    exportgraphics(fh, gif_path, 'Append', true, 'Resolution', 110);
    img = print(fh, '-RGBImage', '-r110');
    if isempty(img_size), img_size = size(img); end
    img = img(1:img_size(1), 1:img_size(2), :);      % keep size constant
    writeVideo(vw, img);
end
close(vw);
close(fh);
fprintf('  saved %s\n  saved %s\n', gif_path, mp4_path);
A.gif = gif_path;  A.mp4 = mp4_path;  A.theta = theta;  A.n_frames = nf;
end
