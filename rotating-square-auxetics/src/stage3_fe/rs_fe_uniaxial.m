function R = rs_fe_uniaxial(M, eps_applied)
%RS_FE_UNIAXIAL  Linear static uniaxial test of the FE lattice along x.
%
%   R = rs_fe_uniaxial(M, eps_applied)   (eps > 0 tension, < 0 compression)
%
%   Boundary conditions (like frictionless grips):
%     - centres of the left column:  u_x = 0
%     - centres of the right column: u_x = eps * (N-1) * p0
%     - one left-column centre (middle row): u_y = 0  (stops rigid sliding)
%     All other DOFs free: edge squares may rotate and slide in y, so
%     the lateral (y) strain is NOT constrained.
%
%   Measurement (avoids edge effects): fit u = eps*x + const on the
%   centres of the INTERIOR squares (outer ring removed when N >= 5):
%       eps_x = d u_x / d x,   eps_y = d u_y / d y,   nu = -eps_y/eps_x
%   Effective modulus from the grip force:
%       sigma = F_right / (N p0 t),   E_eff = sigma / eps_applied
N  = M.N;  p0 = M.p0;
K  = fe_assemble(M);
nd = size(K, 1);

left  = M.centre(:, 1);
right = M.centre(:, N);
mid   = left(ceil(N/2));
delta = eps_applied * (N - 1) * p0;

fixed = [3*left-2; 3*right-2; 3*mid-1];
uval  = [zeros(N,1); delta*ones(N,1); 0];
free  = setdiff((1:nd).', fixed);

U = zeros(nd, 1);
U(fixed) = uval;
U(free)  = K(free, free) \ (-K(free, fixed) * uval);
F = K * U;                                   % reactions at fixed DOFs

% --- strain from interior square centres -------------------------------------
if N >= 5, rng = 2:N-1; else, rng = 1:N; end
ids = M.centre(rng, rng);  ids = ids(:);
x = M.X(ids, 1);  y = M.X(ids, 2);
ux = U(3*ids-2);  uy = U(3*ids-1);
B  = [x y ones(size(x))];
cx = B \ ux;   cy = B \ uy;                  % u = c1 x + c2 y + c3
R.eps_x = cx(1);
R.eps_y = cy(2);
R.nu    = -R.eps_y / R.eps_x;

F_right = sum(F(3*right-2));
R.F_right = F_right;
R.F_left  = sum(F(3*left-2));
R.sigma   = F_right / (N * p0 * M.g.t);
R.E_eff   = R.sigma / eps_applied;
R.U = U;  R.K = K;
end
