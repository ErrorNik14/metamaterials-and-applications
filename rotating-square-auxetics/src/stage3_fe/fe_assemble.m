function K = fe_assemble(M)
%FE_ASSEMBLE  Global sparse stiffness matrix of a frame mesh.
%   M.X (nn x 2) nodes, M.conn (ne x 2), M.E, M.A, M.I (ne x 1).
nn = size(M.X, 1);
ne = size(M.conn, 1);
ii = zeros(36*ne, 1); jj = ii; vv = ii;
for e = 1:ne
    n1 = M.conn(e, 1);  n2 = M.conn(e, 2);
    kg = fe_frame_k(M.E(e), M.A(e), M.I(e), M.X(n1,:), M.X(n2,:));
    dofs = [3*n1-2, 3*n1-1, 3*n1, 3*n2-2, 3*n2-1, 3*n2];
    [r, c] = ndgrid(dofs, dofs);
    idx = (e-1)*36 + (1:36);
    ii(idx) = r(:); jj(idx) = c(:); vv(idx) = kg(:);
end
K = sparse(ii, jj, vv, 3*nn, 3*nn);
K = (K + K.') / 2;          % remove round-off asymmetry
end
