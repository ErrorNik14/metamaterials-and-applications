%% Test: run this file to see an example
clear; clc; close all;
lat = build_lat(-20, 1.3*7, 7, 1.1, 7, 3, true);

%% Function
function lat = build_lat(theta_deg, h, l, t, nx, ny, show)
% BUILD_LATTICE  Re-entrant honeycomb as a 2D frame (nodes + wall list).
%
%   lat = build_lat(theta_deg, h, l, t, nx, ny)
%   lat = build_lat(-20, 1.3*7, 7, 1.1, 7, 3, true)   % plots it
%
%   theta_deg : wall angle from x1 (negative = re-entrant), degrees
%   h, l      : vertical and inclined wall lengths (mm)
%   t         : wall thickness (mm), only used for plotting / output
%   nx, ny    : number of cells along x1 and x2
%
%   lat.nodes : [N x 2] node coordinates (mm)
%   lat.elems : [E x 2] node indices of each wall
%   lat.type  : [E x 1] 1 = vertical (h) wall, 2 = inclined (l) wall
%   lat.size  : [width height] of the lattice (mm)
%
% Unit cell (pitch 2 l cos(theta) by 2(h + l sin(theta))):
%   vertical wall from (0,0) to (0,h), two inclined walls meeting at
%   (l cos, h + l sin), a vertical wall up from there, and two inclined
%   walls back down to the next row's joints.

if nargin < 7, show = false; end

c  = cosd(theta_deg);
s  = sind(theta_deg);
px = 2*l*c;            % cell pitch along x1
py = 2*(h + l*s);      % cell pitch along x2

% walls of one cell, as [x1 y1 x2 y2], with the cell origin at (0,0)
cellW = [ 0     0       0     h       1     % lower vertical wall
          0     h       l*c   h+l*s   2     % lower-left inclined
          l*c   h+l*s   px    h       2     % lower-right inclined
          l*c   h+l*s   l*c   2*h+l*s 1     % middle vertical wall
          l*c   2*h+l*s 0     py      2     % upper-left inclined
          l*c   2*h+l*s px    py      2 ];  % upper-right inclined

segs = [];
for i = 0:nx-1
    for j = 0:ny-1
        segs = [segs; cellW + [i*px j*py i*px j*py 0]]; %#ok<AGROW>
    end
end
% close the right edge with the last column of lower vertical walls
for j = 0:ny-1
    segs = [segs; nx*px j*py nx*px j*py+h 1]; %#ok<AGROW>
end

% merge coincident end points into unique nodes
pts = [segs(:,1:2); segs(:,3:4)];
[nodes, ~, id] = unique(round(pts, 6), 'rows');
E = size(segs, 1);
elems = [id(1:E) id(E+1:end)];

lat.nodes = nodes;
lat.elems = elems;
lat.type  = segs(:,5);
lat.size  = [max(nodes(:,1)) - min(nodes(:,1)), max(nodes(:,2)) - min(nodes(:,2))];
lat.par   = struct('theta', theta_deg, 'h', h, 'l', l, 't', t, 'nx', nx, 'ny', ny);

fprintf('%d nodes, %d walls, lattice %.1f x %.1f mm\n', ...
        size(nodes,1), E, lat.size(1), lat.size(2));

if show
    figure; hold on; axis equal; box on
    for e = 1:E
        xy = nodes(elems(e,:), :);
        plot(xy(:,1), xy(:,2), 'k-', 'LineWidth', 2*t, 'Color', 'green');
    end
    plot(nodes(:,1), nodes(:,2), 'r.', 'MarkerSize', 8, 'Color', 'blue');
    xlabel('x_1 (mm)'); ylabel('x_2 (mm)');
    title(sprintf('\\theta = %g°, h/l = %.2f, l = %g, t = %g mm', ...
          theta_deg, h/l, l, t));
end
end
