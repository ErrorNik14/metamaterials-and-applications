function G = rs_build_array(theta, a, nx, ny)
%RS_BUILD_ARRAY  Build an nx-by-ny array of rotated squares by chaining hinges.
%
%   G = rs_build_array(theta, a, nx, ny)      theta [rad], a [mm]
%
%   This does NOT use the closed-form pitch formula. Each square is placed
%   so that its hinge corner lands exactly on the hinge corner of the square
%   already placed next to it (vector loop closure):
%       c_new = c_old + R(alpha_old) r_old - R(alpha_new) r_new
%   Row 1 is chained left -> right; every other square is chained from the
%   square below it. The horizontal hinges in rows 2..ny are therefore NOT
%   imposed - checking their gap afterwards tests that the mechanism is
%   compatible (that the loops close).
%
%   Outputs (i = column / x index, j = row / y index):
%     G.centre   ny x nx x 2   square centres [mm]
%     G.alpha    ny x nx       rotation of each square [rad]
%     G.corners  ny x nx cell, 4x2 corner coordinates (CCW: TR, TL, BL, BR)
%     G.hinges   H x 2         hinge positions [mm]
%     G.hinge_sq H x 2         linear indices of the two squares at each hinge
%     G.max_hinge_gap          largest distance between the two corners that
%                              are supposed to meet at a hinge [mm]

% Local corner coordinates, counter-clockwise from top-right.
CORNER = (a/2) * [ 1  1;    % 1 TR
                  -1  1;    % 2 TL
                  -1 -1;    % 3 BL
                   1 -1];   % 4 BR

% Pinwheel hinge topology. Type A squares (i+j even) rotate by -theta,
% type B squares rotate by +theta. Each corner of a square is pinned to
% exactly one neighbour; which corner depends on the type:
%                 type A   type B
HINGE.right  = [  1        4 ];    % A: TR,  B: BR
HINGE.left   = [  3        2 ];    % A: BL,  B: TL
HINGE.top    = [  2        1 ];    % A: TL,  B: TR
HINGE.bottom = [  4        3 ];    % A: BR,  B: BL
% (A's right corner meets B's left corner, etc. - the table is consistent:
%  this choice opens the cell for theta > 0; the mirror choice closes it.)

rot  = @(al) [cos(al) -sin(al); sin(al) cos(al)];
type = @(i, j) 1 + mod(i + j, 2);         % 1 = A (i+j even), 2 = B
sgn  = [-1, +1];                          % rotation sign by type

G.alpha  = zeros(ny, nx);
G.centre = zeros(ny, nx, 2);
for j = 1:ny
    for i = 1:nx
        G.alpha(j, i) = sgn(type(i, j)) * theta;
    end
end

% Corner of square (i,j) used for hinge 'side', in global orientation
hinge_vec = @(i, j, side) ...
    rot(G.alpha(j, i)) * CORNER(HINGE.(side)(type(i, j)), :).';

% --- Chain the squares -----------------------------------------------------
for j = 1:ny
    for i = 1:nx
        if i == 1 && j == 1
            c = [0; 0];
        elseif j == 1      % chain from the left neighbour
            c_old = squeeze(G.centre(j, i-1, :));
            c = c_old + hinge_vec(i-1, j, 'right') - hinge_vec(i, j, 'left');
        else               % chain from the neighbour below
            c_old = squeeze(G.centre(j-1, i, :));
            c = c_old + hinge_vec(i, j-1, 'top') - hinge_vec(i, j, 'bottom');
        end
        G.centre(j, i, :) = c;
    end
end

% --- Corners for plotting ----------------------------------------------------
G.corners = cell(ny, nx);
for j = 1:ny
    for i = 1:nx
        c = squeeze(G.centre(j, i, :)).';
        G.corners{j, i} = c + (rot(G.alpha(j, i)) * CORNER.').';
    end
end

% --- Hinge list and compatibility check ---------------------------------------
H = []; HS = []; gap = 0;
idx = @(i, j) sub2ind([ny nx], j, i);
for j = 1:ny
    for i = 1:nx
        cA = squeeze(G.centre(j, i, :));
        if i < nx   % horizontal hinge with the right neighbour
            cB = squeeze(G.centre(j, i+1, :));
            pA = cA + hinge_vec(i, j, 'right');
            pB = cB + hinge_vec(i+1, j, 'left');
            gap = max(gap, norm(pA - pB));
            H(end+1, :) = pA.'; HS(end+1, :) = [idx(i,j) idx(i+1,j)]; %#ok<AGROW>
        end
        if j < ny   % vertical hinge with the top neighbour
            cB = squeeze(G.centre(j+1, i, :));
            pA = cA + hinge_vec(i, j, 'top');
            pB = cB + hinge_vec(i, j+1, 'bottom');
            gap = max(gap, norm(pA - pB));
            H(end+1, :) = pA.'; HS(end+1, :) = [idx(i,j) idx(i,j+1)]; %#ok<AGROW>
        end
    end
end
G.hinges = H;
G.hinge_sq = HS;
G.max_hinge_gap = gap;
G.theta = theta; G.a = a; G.nx = nx; G.ny = ny;
end
