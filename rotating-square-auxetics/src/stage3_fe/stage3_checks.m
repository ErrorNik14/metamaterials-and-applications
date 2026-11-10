function pass = stage3_checks(P, T, Tc, Pt, Pc, Wcmp, S2C)
%STAGE3_CHECKS  Verification of the FE solver, then validation against Stage 2.
%   T/Tc = finite array tension/compression, Pt/Pc = periodic cell
%   tension/compression, Wcmp = width comparison, S2C = Stage-2 design values.
tests = {};
E = hinge_modulus(P);

% 1. element test: cantilever ligament, tip load -> PL^3/(3EI) (exact for EB)
L = P.geom.Lh;  w = P.geom.w;  t = P.geom.t;  I = t*w^3/12;  F = 1;
n = P.fe.n_el_lig;
Mc.X = [linspace(0, L, n+1).' zeros(n+1, 1)];
Mc.conn = [(1:n).' (2:n+1).'];
Mc.E = E*ones(n,1); Mc.A = t*w*ones(n,1); Mc.I = I*ones(n,1);
K = fe_assemble(Mc);
free = 4:3*(n+1);
f = zeros(3*(n+1), 1);  f(end-1) = F;
u = zeros(size(f));  u(free) = K(free, free) \ f(free);
exact = F*L^3/(3*E*I);
tests(end+1,:) = {'cantilever tip vs PL^3/3EI (rel)', abs(u(end-1) - exact)/exact, 1e-10};

% 2. a free 2x2 lattice has exactly 3 rigid-body modes (no mechanism left)
M2 = rs_fe_mesh(fe_geometry_from_params(P), 2);
ev = sort(abs(eig(full(fe_assemble(M2)))));
n_zero = sum(ev < 1e-9 * ev(end));
tests(end+1,:) = {'rigid-body modes of free 2x2 (=3)', abs(n_zero - 3), 0};

% 3. ligaments in the mesh have the designed length
M  = rs_fe_mesh(fe_geometry_from_params(P), P.fe.N_cells);
tests(end+1,:) = {'mesh ligament length - Lh [mm]', max(abs(M.lig_len - P.geom.Lh)), P.tol.geom};

% 4. global equilibrium: left and right grip forces cancel
tests(end+1,:) = {'F_left + F_right (rel)', abs(T.F_left + T.F_right)/abs(T.F_right), 1e-8};

% 5. linear solver: compression is the exact mirror of tension
tests(end+1,:) = {'nu tension - nu compression', abs(T.nu - Tc.nu), 1e-10};

tests(end+1,:) = {'periodic: nu tension - compression', abs(Pt.nu - Pc.nu), 1e-10};

% 6. validation, like for like: periodic FE vs Stage 2 (both infinite lattices)
tests(end+1,:) = {'nu: FE periodic vs Stage 2', abs(Pt.nu - S2C.nu), P.fe.tol_nu_vs_S2};
tests(end+1,:) = {'E_eff: FE periodic vs Stage 2 (rel)', ...
    abs(Pt.E_eff - S2C.E_eff)/S2C.E_eff, P.fe.tol_nu_vs_S2};
tests(end+1,:) = {'nu vs w: max |FE periodic - S2|', ...
    max(abs(Wcmp.nu_FE_per - Wcmp.nu_S2_at_w)), P.fe.tol_nu_vs_S2};

% 7. finite array differs from the infinite lattice only by a small edge effect
tests(end+1,:) = {'finite array: |nu_NxN - nu_periodic|', abs(T.nu - Pt.nu), P.fe.tol_nu_vs_S2};

fprintf('\n  Stage 3 checkpoint tests\n');
fprintf('  %-36s %12s %10s  %s\n', 'test', 'error', 'tol', 'result');
pass = true;
for k = 1:size(tests, 1)
    ok = tests{k,2} <= tests{k,3};
    pass = pass && ok;
    if ok, res = 'PASS'; else, res = 'FAIL'; end
    fprintf('  %-36s %12.3e %10.1e  %s\n', tests{k,1}, tests{k,2}, tests{k,3}, res);
end
end
