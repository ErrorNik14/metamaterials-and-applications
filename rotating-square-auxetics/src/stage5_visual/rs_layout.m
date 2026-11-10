function L = rs_layout(theta, a, Lh, theta0, N, n_pts)
%RS_LAYOUT  Drawable geometry of an N x N array at square rotation theta.
%
%   L = rs_layout(theta, a, Lh, theta0, N, n_pts)     angles in rad
%
%   Kinematics of Stage 2 (rigid squares, ligament arms):
%     pitch      p(theta) = a (cos th + sin th) + Lh cos(th - th0)
%     rotations  -theta (type A) / +theta (type B), checkerboard
%   Ligaments are drawn as cubic (Hermite) curves from the hinge corner of
%   one square to the facing corner of its neighbour. Each end leaves its
%   square at the angle it was printed, turned with that square - this is
%   the bent shape of a beam whose two ends are clamped into the squares.
%
%   Output: L.corners{j,i} (4x2), L.lig{k} (n_pts x 2), L.centre (N*N x 2),
%           L.p (pitch), L.eps (strain relative to theta0)
rot = @(al) [cos(al) -sin(al); sin(al) cos(al)];
CORNER = (a/2) * [1 1; -1 1; -1 -1; 1 -1];             % TR TL BL BR
HINGE  = struct('right', [1 4], 'left', [3 2], 'top', [2 1], 'bottom', [4 3]);
sgn = [-1 +1];

L.p   = a*(cos(theta) + sin(theta)) + Lh*cos(theta - theta0);
p0    = a*(cos(theta0) + sin(theta0)) + Lh;
L.eps = L.p / p0 - 1;

typ = zeros(N);  alpha = zeros(N);  cen = zeros(N, N, 2);
L.corners = cell(N);
for j = 1:N
    for i = 1:N
        typ(j,i)   = 1 + mod(i + j, 2);
        alpha(j,i) = sgn(typ(j,i)) * theta;
        c = [(i-1) (j-1)] * L.p;
        cen(j,i,:) = c;
        L.corners{j,i} = c + (rot(alpha(j,i)) * CORNER.').';
    end
end
L.centre = reshape(cen, [], 2);

% Hermite basis on n_pts points
s = linspace(0, 1, n_pts).';
h00 = 2*s.^3 - 3*s.^2 + 1;  h10 = s.^3 - 2*s.^2 + s;
h01 = -2*s.^3 + 3*s.^2;     h11 = s.^3 - s.^2;

L.lig = {};  L.lig_turn = [];  L.lig_chord = [];
for j = 1:N
    for i = 1:N
        nb = {};
        if i < N, nb(end+1,:) = {i+1, j, 'right', 'left', [1 0]}; end %#ok<AGROW>
        if j < N, nb(end+1,:) = {i, j+1, 'top', 'bottom', [0 1]}; end %#ok<AGROW>
        for q = 1:size(nb, 1)
            [i2, j2, s1, s2, e] = nb{q, :};
            P0 = L.corners{j,  i }(HINGE.(s1)(typ(j, i)),   :);
            P1 = L.corners{j2, i2}(HINGE.(s2)(typ(j2, i2)), :);
            dA = alpha(j, i)   - sgn(typ(j, i))   * theta0;   % turn since printing
            dB = alpha(j2, i2) - sgn(typ(j2, i2)) * theta0;
            t0 = Lh * (rot(dA) * e(:)).';
            t1 = Lh * (rot(dB) * e(:)).';
            L.lig{end+1} = h00*P0 + h10*t0 + h01*P1 + h11*t1;
            L.lig_turn(end+1) = dB - dA;                      % relative rotation
            L.lig_chord(end+1) = norm(P1 - P0);              % corner-to-corner
        end
    end
end
end
