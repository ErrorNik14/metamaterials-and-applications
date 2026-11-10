function pass = stage2_checks(P, C0, R0, SW, NL)
%STAGE2_CHECKS  Checkpoint tests for Stage 2. Prints PASS/FAIL per test.
E   = hinge_modulus(P);
g   = P.geom;
th0 = deg2rad(g.theta0_deg);
rf  = P.model.root_length_factor;
tests = {};

% 1. closed form vs independent numerical cell (nu and E)
tests(end+1,:) = {'nu: formula vs numerical cell', abs(C0.nu - R0.nu), P.tol.rel};
tests(end+1,:) = {'E_eff: formula vs numerical (rel)', ...
    abs(C0.E_eff - R0.E_eff) / C0.E_eff, P.tol.rel};
% 2. same comparison at a deliberately thick, different design
Ct = rs_compliant_props(8, 2, 3, g.t, E, deg2rad(25), rf);
Rt = rs_energy_check(8, 2, 3, g.t, E, deg2rad(25), rf);
tests(end+1,:) = {'nu: thick-hinge design', abs(Ct.nu - Rt.nu), P.tol.rel};
% 3. geometry of the numerical cell closes in the as-printed state
tests(end+1,:) = {'cell gap as printed [mm]', R0.max_gap_as_printed, P.tol.geom};
% 4. thin-hinge limit recovers Stage 1: nu -> -1
Cthin = rs_compliant_props(g.a, g.Lh, 1e-4, g.t, E, th0, rf);
tests(end+1,:) = {'w -> 0 gives nu = -1', abs(Cthin.nu + 1), P.tol.nu};
% 5. nu does not depend on E or t (pure geometry ratio)
C2 = rs_compliant_props(g.a, g.Lh, g.w, 3*g.t, 7*E, th0, rf);
tests(end+1,:) = {'nu independent of E and t', abs(C2.nu - C0.nu), 1e-12};
% 6. E_eff scales linearly with E and is independent of t
tests(end+1,:) = {'E_eff/E independent of E, t', abs(C2.E_ratio - C0.E_ratio)/C0.E_ratio, 1e-12};
% 7. nonlinear curve has slope E_eff at the as-printed state
tests(end+1,:) = {'NL slope at eps=0 vs E_eff (rel)', ...
    abs(NL.slope0 - C0.E_eff_rotation_only) / C0.E_eff_rotation_only, P.tol.rel};
% 8. sweeps are monotonic as theory says (stiffer with wider, softer with longer)
tests(end+1,:) = {'E_eff rises with w (# violations)', sum(diff(SW.vs_w.E_eff) <= 0), 0};
tests(end+1,:) = {'E_eff falls with Lh (# violations)', sum(diff(SW.vs_Lh.E_eff) >= 0), 0};

fprintf('\n  Stage 2 checkpoint tests\n');
fprintf('  %-36s %12s %10s  %s\n', 'test', 'error', 'tol', 'result');
pass = true;
for k = 1:size(tests, 1)
    ok = tests{k,2} <= tests{k,3};
    pass = pass && ok;
    if ok, res = 'PASS'; else, res = 'FAIL'; end
    fprintf('  %-36s %12.3e %10.1e  %s\n', tests{k,1}, tests{k,2}, tests{k,3}, res);
end
end
