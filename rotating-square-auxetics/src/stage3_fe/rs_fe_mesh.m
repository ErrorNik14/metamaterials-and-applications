function M = rs_fe_mesh(g, N, periodic)
%RS_FE_MESH  Frame-element mesh of an N x N rotating-squares array.
%
%   M = rs_fe_mesh(g, N)            finite array (free edges)
%   M = rs_fe_mesh(g, N, true)      periodic cell (N even): the right column
%       is also joined to the left column, and the top row to the bottom row,
%       by "wrap-around" ligaments. M.wrap(e,:) = [sx sy] says that node 2 of
%       element e really sits at X + [sx*Lx, sy*Ly] (its periodic image).
%
%   g fields: a, Lh, w, t, theta0 [rad], E, n_el_lig, square_E_factor
%
%   Each SQUARE  = 1 centre node + 4 corner nodes, joined by 4 stiff
%                  beams (a rigid "star": rigid if its beams are stiff).
%   Each LIGAMENT = n_el_lig PLA beam elements (width w, depth t) running
%                  from a corner of one square to the facing corner of its
%                  neighbour. Its ends share the corner nodes (welded).
%
%   Same geometry and hinge table as Stages 1-2: pitch
%   p0 = a(cos th0 + sin th0) + Lh, squares rotated -/+ th0 (checkerboard).
%
%   Output: M.X, M.conn, M.E, M.A, M.I, M.is_lig (element flags),
%           M.centre (N x N node ids, row j, column i), M.corner (N x N x 4),
%           M.lig_len (actual ligament lengths, for checking)
if nargin < 3, periodic = false; end
if periodic && mod(N, 2) ~= 0
    error('A periodic cell needs an even N (the checkerboard repeats every 2).');
end
rot = @(al) [cos(al) -sin(al); sin(al) cos(al)];
CORNER = (g.a/2) * [1 1; -1 1; -1 -1; 1 -1];        % TR TL BL BR
HINGE  = struct('right', [1 4], 'left', [3 2], 'top', [2 1], 'bottom', [4 3]);
sgn = [-1 +1];
p0  = g.a * (cos(g.theta0) + sin(g.theta0)) + g.Lh;

X = zeros(0, 2);
M.centre = zeros(N, N);
M.corner = zeros(N, N, 4);
M.type   = zeros(N, N);
conn = zeros(0, 2);  is_lig = false(0, 1);

% --- squares -----------------------------------------------------------------
for j = 1:N
    for i = 1:N
        ty = 1 + mod(i + j, 2);
        c  = [(i-1)*p0, (j-1)*p0];
        Rm = rot(sgn(ty) * g.theta0);
        X(end+1, :) = c; %#ok<AGROW>
        M.centre(j, i) = size(X, 1);
        M.type(j, i) = ty;
        for k = 1:4
            X(end+1, :) = c + (Rm * CORNER(k, :).').'; %#ok<AGROW>
            M.corner(j, i, k) = size(X, 1);
            conn(end+1, :) = [M.centre(j, i), size(X, 1)]; %#ok<AGROW>
            is_lig(end+1, 1) = false; %#ok<AGROW>
        end
    end
end

% --- ligaments -----------------------------------------------------------------
M.lig_len = [];
M.lig_elems = {};                       % element ids of each ligament
Lper = N * p0;                          % period (used only if periodic)
wrap = zeros(size(conn, 1), 2);         % square arms never wrap
for j = 1:N
    for i = 1:N
        % {neighbour i, neighbour j, side, other side, periodic shift}
        nbrs = {};
        if i < N, nbrs(end+1, :) = {i+1, j, 'right', 'left', [0 0]};      %#ok<AGROW>
        elseif periodic, nbrs(end+1, :) = {1, j, 'right', 'left', [1 0]}; %#ok<AGROW>
        end
        if j < N, nbrs(end+1, :) = {i, j+1, 'top', 'bottom', [0 0]};      %#ok<AGROW>
        elseif periodic, nbrs(end+1, :) = {i, 1, 'top', 'bottom', [0 1]}; %#ok<AGROW>
        end
        for q = 1:size(nbrs, 1)
            [i2, j2, s1, s2, sh] = nbrs{q, :};
            nA = M.corner(j,  i,  HINGE.(s1)(M.type(j, i)));
            nB = M.corner(j2, i2, HINGE.(s2)(M.type(j2, i2)));
            xB = X(nB,:) + sh * Lper;            % image position of the far end
            M.lig_len(end+1, 1) = norm(xB - X(nA,:));
            % interior nodes along the ligament
            ids = nA;
            for m = 1:g.n_el_lig - 1
                X(end+1, :) = X(nA,:) + (m / g.n_el_lig) * (xB - X(nA,:)); %#ok<AGROW>
                ids(end+1) = size(X, 1); %#ok<AGROW>
            end
            ids(end+1) = nB; %#ok<AGROW>
            first = size(conn, 1) + 1;
            for m = 1:g.n_el_lig
                conn(end+1, :) = ids(m:m+1); %#ok<AGROW>
                is_lig(end+1, 1) = true; %#ok<AGROW>
                wrap(end+1, :) = sh * (m == g.n_el_lig); %#ok<AGROW>
            end
            M.lig_elems{end+1, 1} = first:size(conn, 1);
        end
    end
end
M.wrap = wrap;  M.Lper = Lper;  M.periodic = periodic;

% --- section properties ----------------------------------------------------------
ne = size(conn, 1);
M.X = X;  M.conn = conn;  M.is_lig = is_lig;
M.E = g.E * ones(ne, 1);
M.A = g.t * g.w * ones(ne, 1);
M.I = g.t * g.w^3 / 12 * ones(ne, 1);
% square arms: depth a (a solid square is far deeper than a ligament),
% and E multiplied by square_E_factor -> rigid compared with the hinges
M.E(~is_lig) = g.E * g.square_E_factor;
M.A(~is_lig) = g.t * g.a;
M.I(~is_lig) = g.t * g.a^3 / 12;
M.p0 = p0;  M.N = N;  M.g = g;
end
