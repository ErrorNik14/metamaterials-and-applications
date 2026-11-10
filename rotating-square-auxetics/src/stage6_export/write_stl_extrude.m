function TR = write_stl_extrude(ps, t, file)
%WRITE_STL_EXTRUDE  Extrude a polyshape by thickness t and save as STL (mm).
%   Bottom face (z = 0, normals down), top face (z = t, normals up) from the
%   polyshape triangulation, plus one quad (2 triangles) per boundary edge.
%   Base MATLAB: triangulation + stlwrite. Returns the 3D triangulation.
T2 = triangulation(ps);
V2 = T2.Points;  F2 = T2.ConnectivityList;
nv = size(V2, 1);
V = [V2 zeros(nv,1); V2 t*ones(nv,1)];

% make every 2D triangle counter-clockwise (positive signed area)
a2 = (V2(F2(:,2),1) - V2(F2(:,1),1)) .* (V2(F2(:,3),2) - V2(F2(:,1),2)) - ...
     (V2(F2(:,3),1) - V2(F2(:,1),1)) .* (V2(F2(:,2),2) - V2(F2(:,1),2));
F2(a2 < 0, :) = F2(a2 < 0, [1 3 2]);
F = [F2(:, [1 3 2]); F2 + nv];                     % bottom faces down, top up

% side walls: for edge p1 -> p2 the triangle (p1, p2, p2_top) has normal
% along n = (dy, -dx). Keep it if n points OUT of the material, else flip.
h = 1e-4 * sqrt(area(ps));                         % small probe distance
for b = 1:numboundaries(ps)
    [x, y] = boundary(ps, b);
    for k = 1:numel(x) - 1
        i1 = find_vertex(V2, [x(k)   y(k)]);
        i2 = find_vertex(V2, [x(k+1) y(k+1)]);
        d  = [x(k+1) - x(k), y(k+1) - y(k)];
        n  = [d(2), -d(1)] / norm(d);
        mid = [x(k) + x(k+1), y(k) + y(k+1)] / 2;
        if isinterior(ps, mid(1) - h*n(1), mid(2) - h*n(2))   % material behind n
            F(end+1, :) = [i1, i2, i2 + nv]; %#ok<AGROW>
            F(end+1, :) = [i1, i2 + nv, i1 + nv]; %#ok<AGROW>
        else
            F(end+1, :) = [i1, i2 + nv, i2]; %#ok<AGROW>
            F(end+1, :) = [i1, i1 + nv, i2 + nv]; %#ok<AGROW>
        end
    end
end
TR = triangulation(F, V);
stlwrite(TR, file);
end

function i = find_vertex(V, p)
[d, i] = min(sum((V - p).^2, 2));
if d > 1e-12
    error('Boundary point not found in triangulation.');
end
end
