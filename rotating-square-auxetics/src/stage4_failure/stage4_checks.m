function pass = stage4_checks(P, S4, S2)
%STAGE4_CHECKS  Checkpoint tests for Stage 4.
g = P.geom;  E = hinge_modulus(P);
tests = {};

% 1. geometric stiffness: pinned-pinned column -> pi^2 EI / L^2
L = 50; w = 2; t = 5; n = 16; I = t*w^3/12; A = t*w;
X = [linspace(0, L, n+1).' zeros(n+1,1)];
K = sparse(3*(n+1), 3*(n+1));  KG = K;
for e = 1:n
    d = 3*e-2 : 3*e+3;
    K(d,d)  = K(d,d)  + fe_frame_k(E, A, I, X(e,:), X(e+1,:));
    KG(d,d) = KG(d,d) + fe_frame_kg(-1, X(e,:), X(e+1,:));     % unit compression
end
fixed = [1 2 3*(n+1)-1];                       % pin left (u,v), roller right (v)
free = setdiff(1:3*(n+1), fixed);
lam = eig(full(K(free,free)), -full(KG(free,free)));
lam = min(lam(lam > 0));
exact = pi^2 * E * I / L^2;
tests(end+1,:) = {'column buckling vs pi^2EI/L^2 (rel)', abs(lam - exact)/exact, 1e-3};

% 2. hinge stress model vs FE (centre ligaments of a big array, 1 % strain)
Mb = rs_fe_mesh(fe_geometry_from_params(P), P.fe.N_convergence(end));
Rb = rs_fe_uniaxial(Mb, P.fe.applied_strain);
[~, s_c] = fe_ligament_stress(Mb, Rb.U, P.fail.Kt);
Hb = rs_hinge_stress(g.a, g.Lh, g.w, g.t, E, deg2rad(g.theta0_deg), ...
    P.model.root_length_factor, P.fail.Kt, Rb.eps_x);
tests(end+1,:) = {'hinge stress: theory vs FE centre', abs(Hb.sig - s_c)/s_c, 0.01};

% 3. ligament axial force from equilibrium matches the axial-spring value
C = S2.C0;  ea = 0.01;
dth = rs_hinge_stress(g.a, g.Lh, g.w, g.t, E, deg2rad(g.theta0_deg), ...
    P.model.root_length_factor, 1, ea).dtheta_per_eps * ea;
N_spring = C.k_a * (C.p0*ea - C.b*dth);
N_equil  = C.E_eff * ea * C.p0 * g.t;
tests(end+1,:) = {'ligament force: spring vs equilibrium', abs(N_spring - N_equil)/N_equil, 1e-6};

% 4. no contact in the as-printed state
tests(end+1,:) = {'as printed: squares do not touch', double(S4.contact.gap0 <= 0), 0};

% 5. safe range is inside every individual limit
inside = S4.safe(2) <= S4.ycross(2) && S4.safe(2) <= S4.lock_open && ...
         S4.safe(1) >= S4.ycross(1) && S4.safe(1) >= S4.contact.eps && ...
         S4.safe(1) >= S4.buckle.eps_cr(1);
tests(end+1,:) = {'safe range inside all limits', double(~inside), 0};

fprintf('\n  Stage 4 checkpoint tests\n');
fprintf('  %-38s %12s %10s  %s\n', 'test', 'error', 'tol', 'result');
pass = true;
for k = 1:size(tests, 1)
    ok = tests{k,2} <= tests{k,3};
    pass = pass && ok;
    if ok, res = 'PASS'; else, res = 'FAIL'; end
    fprintf('  %-38s %12.3e %10.1e  %s\n', tests{k,1}, tests{k,2}, tests{k,3}, res);
end
end
