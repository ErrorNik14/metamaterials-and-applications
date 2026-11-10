function B = rs_fe_buckling(M, eps_ref, n_modes)
%RS_FE_BUCKLING  Linear (eigenvalue) buckling of the finite FE array, in-plane.
%
%   1) Linear solve at a compressive reference strain eps_ref (< 0), same
%      grips as rs_fe_uniaxial -> axial force N_e in every element.
%   2) Assemble the geometric stiffness KG from those forces.
%   3) Solve  (K + lambda * KG) phi = 0  on the free DOFs.
%      Critical strain = lambda_min * eps_ref.
%   This is the same method as an Ansys "Eigenvalue Buckling" analysis.
%   It finds IN-PLANE bifurcations only (2D model): out-of-plane buckling of
%   the sheet is not captured.
R  = rs_fe_uniaxial(M, eps_ref);
K  = R.K;  U = R.U;
nd = size(K, 1);
ne = size(M.conn, 1);

ii = zeros(36*ne,1); jj = ii; vv = ii;
B.N = zeros(ne, 1);
for e = 1:ne
    n1 = M.conn(e,1);  n2 = M.conn(e,2);
    dofs = [3*n1-2, 3*n1-1, 3*n1, 3*n2-2, 3*n2-1, 3*n2];
    [~, T, kl] = fe_frame_k(M.E(e), M.A(e), M.I(e), M.X(n1,:), M.X(n2,:));
    f = kl * (T * U(dofs));                   % local end forces
    B.N(e) = -f(1);                           % axial force, tension +
    kgG = fe_frame_kg(B.N(e), M.X(n1,:), M.X(n2,:));
    [r, c] = ndgrid(dofs, dofs);
    idx = (e-1)*36 + (1:36);
    ii(idx) = r(:); jj(idx) = c(:); vv(idx) = kgG(:);
end
KG = sparse(ii, jj, vv, nd, nd);  KG = (KG + KG.')/2;

% same supports as the static solve (grips hold u_x; one u_y pinned)
N = M.N;
left = M.centre(:,1); right = M.centre(:,N); mid = left(ceil(N/2));
fixed = [3*left-2; 3*right-2; 3*mid-1];
free  = setdiff((1:nd).', fixed);
Kf  = full(K(free,free));  KGf = full(KG(free,free));

% K phi = -lambda KG phi  <=>  (-KG) phi = mu K phi,  mu = 1/lambda
[V, D] = eig(-KGf, Kf, 'chol');
mu = diag(D);
pos = find(mu > 0);
[mu_s, ord] = sort(mu(pos), 'descend');
nm = min(n_modes, numel(mu_s));
B.lambda = 1 ./ mu_s(1:nm);
B.eps_cr = B.lambda * eps_ref;
B.modes  = zeros(nd, nm);
for k = 1:nm
    B.modes(free, k) = V(:, pos(ord(k)));
end
B.eps_ref = eps_ref;
end
