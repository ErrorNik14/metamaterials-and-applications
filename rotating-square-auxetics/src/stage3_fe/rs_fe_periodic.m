function R = rs_fe_periodic(M, sigma_x)
%RS_FE_PERIODIC  Uniaxial stress on a PERIODIC FE cell (infinite lattice).
%
%   R = rs_fe_periodic(M, sigma_x)     M from rs_fe_mesh(g, N, true)
%
%   Unknowns: all nodal DOFs  u  plus the two average strains eps_x, eps_y.
%   A wrap-around element sees its far node at the periodic image, whose
%   displacement is  u_image = u_node + [eps_x * Lper, 0]   (x wrap)
%                                     + [0, eps_y * Lper]   (y wrap).
%   So each element's DOFs are a linear map  u_e = C * [u; eps_x; eps_y],
%   and its stiffness enters the global matrix as  C' * k_e * C.
%   Load: the force conjugate to eps_x is sigma_x * (cell area * t);
%         the one conjugate to eps_y is 0 (sigma_y = 0, uniaxial stress).
%   One node is pinned in u, v to remove rigid translation. Rigid rotation
%   is already impossible (it is not periodic without a shear strain DOF).
nn  = size(M.X, 1);
ne  = size(M.conn, 1);
nd  = 3*nn + 2;
iex = 3*nn + 1;  iey = 3*nn + 2;          % indices of eps_x, eps_y
L   = M.Lper;

ii = zeros(64*ne, 1); jj = ii; vv = ii; ptr = 0;
for e = 1:ne
    n1 = M.conn(e,1);  n2 = M.conn(e,2);
    x2 = M.X(n2,:) + M.wrap(e,:) * L;
    kg = fe_frame_k(M.E(e), M.A(e), M.I(e), M.X(n1,:), x2);
    dofs = [3*n1-2, 3*n1-1, 3*n1, 3*n2-2, 3*n2-1, 3*n2, iex, iey];
    C = [eye(6), zeros(6, 2)];
    C(4, 7) = M.wrap(e,1) * L;            % u2 += eps_x * L  (x wrap)
    C(5, 8) = M.wrap(e,2) * L;            % v2 += eps_y * L  (y wrap)
    ke = C.' * kg * C;
    [r, c] = ndgrid(dofs, dofs);
    idx = ptr + (1:64);  ptr = ptr + 64;
    ii(idx) = r(:); jj(idx) = c(:); vv(idx) = ke(:);
end
K = sparse(ii, jj, vv, nd, nd);
K = (K + K.') / 2;

V = L^2 * M.g.t;
f = zeros(nd, 1);  f(iex) = sigma_x * V;
pin   = M.centre(1, 1);
fixed = [3*pin-2; 3*pin-1];
free  = setdiff((1:nd).', fixed);
U = zeros(nd, 1);
U(free) = K(free, free) \ f(free);

R.eps_x = U(iex);  R.eps_y = U(iey);
R.nu    = -R.eps_y / R.eps_x;
R.E_eff = sigma_x / R.eps_x;
R.U = U(1:3*nn);
end
