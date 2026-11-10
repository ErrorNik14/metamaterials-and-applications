function draw_lattice(ax, P, L, L0, lig_rgb)
%DRAW_LATTICE  Draw a deformed layout L over the as-printed outline L0.
%   Both are centred on the origin so the lattice grows/shrinks about
%   its middle. lig_rgb colours the ligaments (used for stress status).
cla(ax); hold(ax, 'on');
c  = mean(L.centre, 1);
c0 = mean(L0.centre, 1);
for k = 1:numel(L0.corners)                      % as-printed outline (ghost)
    C = L0.corners{k} - c0;
    patch(ax, C(:,1), C(:,2), 'w', 'FaceColor', 'none', ...
        'EdgeColor', P.plot.c_ghost, 'LineStyle', '--', 'LineWidth', 0.8);
end
for k = 1:numel(L.corners)                       % squares
    C = L.corners{k} - c;
    patch(ax, C(:,1), C(:,2), 'w', 'FaceColor', P.plot.c_square, ...
        'EdgeColor', P.plot.c_series1, 'LineWidth', 0.9);
end
for k = 1:numel(L.lig)                           % ligaments
    Q = L.lig{k} - c;
    plot(ax, Q(:,1), Q(:,2), '-', 'Color', lig_rgb, 'LineWidth', 2.5);
end
axis(ax, 'equal');  axis(ax, 'off');
end
