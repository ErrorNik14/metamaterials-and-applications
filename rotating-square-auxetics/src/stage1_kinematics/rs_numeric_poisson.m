function N = rs_numeric_poisson(theta, a, n, theta0)
%RS_NUMERIC_POISSON  Poisson's ratio measured from built arrays (no formula).
%
%   N = rs_numeric_poisson(theta, a, n, theta0)
%   For each angle, build an n x n array with rs_build_array, measure the
%   average pitch in x and in y from the square CENTRES, then
%     nu_incremental = -(dLy/dtheta)/(dLx/dtheta)  (central differences)
%     nu_secant      = -eps_y/eps_x  relative to the array at theta0
%   Centres are used (not the outer bounding box) because the outer corners
%   of the edge squares stick out by an amount that depends on theta - that
%   is an edge effect, not a property of the material.

theta = theta(:);
nt = numel(theta);
N.theta = theta;
N.px = zeros(nt, 1);  N.py = zeros(nt, 1);  N.max_gap = 0;

for k = 1:nt
    [N.px(k), N.py(k), gap] = pitch_of(theta(k));
    N.max_gap = max(N.max_gap, gap);
end
[px0, py0] = pitch_of(theta0);

N.eps_x = N.px / px0 - 1;
N.eps_y = N.py / py0 - 1;

dLx = gradient(N.px, theta);
dLy = gradient(N.py, theta);
N.nu_incremental = -dLy ./ dLx;
N.nu_incremental(abs(dLx) < 1e-6 * a) = NaN;   % undefined at the 45 deg lock

N.nu_secant = -N.eps_y ./ N.eps_x;
N.nu_secant(abs(N.eps_x) < 1e-9) = NaN;        % undefined at theta0

    function [px, py, gap] = pitch_of(th)
        G  = rs_build_array(th, a, n, n);
        cx = G.centre(:, :, 1);
        cy = G.centre(:, :, 2);
        % mean centre spacing along a row (x) and up a column (y)
        px = mean(cx(:, end) - cx(:, 1)) / (n - 1);
        py = mean(cy(end, :) - cy(1, :)) / (n - 1);
        gap = G.max_hinge_gap;
    end
end
