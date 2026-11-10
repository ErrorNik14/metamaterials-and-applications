function pass = stage1_checks(P, K, Nm)
%STAGE1_CHECKS  Checkpoint tests for Stage 1. Prints PASS/FAIL per test.
a      = P.geom.a;
theta0 = deg2rad(P.geom.theta0_deg);
ok_range = K.theta > 0 & K.theta < K.theta_lock_open;   % open, not locked

tests = {};
% 1. ideal Poisson's ratio is -1 (closed form)
err = max(abs(K.nu_incremental(ok_range) + 1));
tests(end+1, :) = {'nu = -1 (closed form)', err, P.tol.nu};
% 2. ideal Poisson's ratio is -1 (numerical array, independent of formula)
okN = Nm.theta > 0 & Nm.theta < K.theta_lock_open;
err = max(abs(Nm.nu_incremental(okN) + 1), [], 'omitnan');
tests(end+1, :) = {'nu = -1 (numerical array)', err, P.tol.nu};
% 3. every hinge closes: the mechanism is compatible
tests(end+1, :) = {'hinge gap in array [mm]', Nm.max_gap, P.tol.geom};
% 4. numerical pitch equals closed form a(cos+sin)
err = max(abs(Nm.px - a*(cos(Nm.theta) + sin(Nm.theta))));
tests(end+1, :) = {'pitch numeric vs formula [mm]', err, P.tol.geom};
% 5. limiting cases: closed theta=0 -> pitch a; open 45 -> pitch a*sqrt(2)
Kl = rs_cell_kinematics([0; pi/4], a, theta0);
err = max(abs(Kl.Lx/2 - [a; a*sqrt(2)]));
tests(end+1, :) = {'pitch at 0 and 45 deg [mm]', err, P.tol.geom};
% 6. porosity at the open lock is exactly 50 %
err = abs(Kl.porosity(2) - 0.5);
tests(end+1, :) = {'porosity at 45 deg = 0.5', err, P.tol.geom};
% 7. inverse kinematics recovers theta
th_test = linspace(0.05, K.theta_lock_open - 0.05, 20).';
Kt = rs_cell_kinematics(th_test, a, theta0);
err = max(abs(rs_theta_from_strain(Kt.eps_x, theta0) - th_test));
tests(end+1, :) = {'theta(eps) inverse [rad]', err, 1e-10};

fprintf('\n  Stage 1 checkpoint tests\n');
fprintf('  %-34s %12s %10s  %s\n', 'test', 'error', 'tol', 'result');
pass = true;
for k = 1:size(tests, 1)
    ok = tests{k, 2} <= tests{k, 3};
    pass = pass && ok;
    fprintf('  %-34s %12.3e %10.1e  %s\n', tests{k,1}, tests{k,2}, tests{k,3}, ...
        ternary(ok, 'PASS', 'FAIL'));
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
