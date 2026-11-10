function h = draw_squares(ax, G, centre_shift, face_rgb, edge_rgb, line_style, show_hinges)
%DRAW_SQUARES  Draw every square of an array built by rs_build_array.
%   centre_shift [dx dy] moves the drawing (used to centre arrays on 0,0).
%   face_rgb = 'none' draws outlines only (used for the undeformed ghost).
if nargin < 7, show_hinges = false; end
hold(ax, 'on');
h = gobjects(0);
for k = 1:numel(G.corners)
    C = G.corners{k} + centre_shift;
    h(end+1) = patch(ax, C(:,1), C(:,2), 'w', ...
        'FaceColor', face_rgb, 'EdgeColor', edge_rgb, ...
        'LineStyle', line_style, 'LineWidth', 0.8); %#ok<AGROW>
end
if show_hinges
    Hp = G.hinges + centre_shift;
    plot(ax, Hp(:,1), Hp(:,2), 'o', 'MarkerSize', 3, ...
        'MarkerFaceColor', edge_rgb, 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
end
end
