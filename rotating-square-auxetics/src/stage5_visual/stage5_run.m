function S5 = stage5_run(P, S2, S4)
%STAGE5_RUN  Stage 5: animation (GIF/MP4), layout checks. The slider app is
%   rs_explorer_app.m (run it yourself - it needs a display).
fprintf('\n=== Stage 5: visual outputs ===\n');
g   = P.geom;
th0 = deg2rad(g.theta0_deg);
tests = {};

% 1. drawn layout at theta0 = the FE mesh (same squares, same corners)
M  = rs_fe_mesh(fe_geometry_from_params(P), P.anim.N);
L0 = rs_layout(th0, g.a, g.Lh, th0, P.anim.N, P.anim.lig_points);
err = 0;
for j = 1:P.anim.N
    for i = 1:P.anim.N
        err = max(err, max(abs(L0.corners{j,i} - M.X(squeeze(M.corner(j,i,:)), :)), [], 'all'));
    end
end
tests(end+1,:) = {'layout at theta0 = FE mesh [mm]', err, P.tol.geom};

% 2. corner-to-corner ligament chord = Lh cos(turn/2) (rigid-arm kinematics)
th = th0 + deg2rad(10);
L  = rs_layout(th, g.a, g.Lh, th0, P.anim.N, P.anim.lig_points);
tests(end+1,:) = {'ligament chord = Lh cos(dphi/2) [mm]', ...
    max(abs(L.lig_chord - g.Lh*cos(L.lig_turn/2))), P.tol.geom};

% 3. state at theta0: zero strain, zero stress
S0 = rs_state_at_theta(P, th0);
tests(end+1,:) = {'state at theta0: |eps| + |sigma|', abs(S0.eps) + abs(S0.sigma), 1e-12};

% 4. app/animation strain matches Stage 2 path
k  = round(numel(S2.NL.theta) * 0.7);
Sk = rs_state_at_theta(P, S2.NL.theta(k));
tests(end+1,:) = {'state strain vs Stage 2 path', abs(Sk.eps - S2.NL.eps(k)), 1e-12};

% 5. animation files written
A = stage5_animation(P, S2, S4);
info = imfinfo(A.gif);
tests(end+1,:) = {'GIF frame count - n_frames', abs(numel(info) - A.n_frames), 0};
vr = VideoReader(A.mp4);
tests(end+1,:) = {'MP4 duration - expected [s]', abs(vr.Duration - A.n_frames/P.anim.fps), 0.1};

% 6. a still of the slider view for the slides (as-printed and near lock)
fh = new_figure(P, P.plot.fig_width_cm, P.plot.fig_width_cm/2 + 1);
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for th_show = [th0, th0 + P.anim.open_frac*(S2.NL.theta_lock - th0)]
    ax = nexttile(tl);
    St = rs_state_at_theta(P, th_show, S4.edge_factor);
    draw_lattice(ax, P, rs_layout(th_show, g.a, g.Lh, th0, P.anim.N, P.anim.lig_points), L0, St.rgb);
    title(ax, sprintf('\\theta = %.1f\\circ, \\epsilon = %+.1f %% [%s]', ...
        rad2deg(th_show), 100*St.eps, St.status), 'FontSize', P.plot.title_size);
end
save_figure(fh, 'S5_printed_vs_open', P);

fprintf('\n  Stage 5 checkpoint tests\n');
fprintf('  %-38s %12s %10s  %s\n', 'test', 'error', 'tol', 'result');
pass = true;
for k = 1:size(tests, 1)
    ok = tests{k,2} <= tests{k,3};
    pass = pass && ok;
    if ok, res = 'PASS'; else, res = 'FAIL'; end
    fprintf('  %-38s %12.3e %10.1e  %s\n', tests{k,1}, tests{k,2}, tests{k,3}, res);
end
S5.anim = A;  S5.pass = pass;
if pass
    fprintf('  Stage 5: ALL CHECKS PASSED\n');
else
    warning('Stage 5: at least one check FAILED - see table above.');
end
end
