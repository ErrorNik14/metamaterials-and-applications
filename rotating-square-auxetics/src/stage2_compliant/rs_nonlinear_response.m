function NL = rs_nonlinear_response(a, Lh, w, t, E, theta0, root_factor, lock_margin)
%RS_NONLINEAR_RESPONSE  Large-rotation stress-strain curve, rotation mode only.
%
%   The lattice is driven through the mechanism path (ligaments do not
%   stretch, rho -> 0) with a LINEAR rotational spring at every hinge.
%   Geometry is exact (large rotation), the hinge spring is not
%   (beam theory is small-rotation) - so trust the shape, not the last digit.
%
%   Pitch with ligament arms (derived in docs/stage2_derivation.md):
%       p(theta)  = a (cos th + sin th) + Lh cos(th - th0)
%       p'(theta) = a (cos th - sin th) - Lh sin(th - th0)
%   Energy per 2x2 cell  U = 8 * 1/2 k_r (2 dtheta)^2 = 16 k_r dtheta^2
%   Work by sigma_x only (sigma_y = 0):  sigma * (Ly0 t) dLx = dU
%       sigma(theta) = 8 k_r (th - th0) / (p0 t p'(th))     [nominal, MPa]

L   = Lh + root_factor * w;
k_r = E * t * w^3 / 12 / L;

p   = @(th) a*(cos(th) + sin(th)) + Lh*cos(th - theta0);
dp  = @(th) a*(cos(th) - sin(th)) - Lh*sin(th - theta0);
p0  = p(theta0);
b0  = dp(theta0);

% Kinematic (open) lock: dp/dtheta = 0
NL.theta_lock = fzero(dp, [theta0, pi/2]);
NL.eps_lock   = p(NL.theta_lock) / p0 - 1;

% Plot path: from fully closed (theta = 0) to just before the lock
th_end = fzero(@(th) dp(th) - lock_margin*b0, [theta0, NL.theta_lock]);
NL.theta = linspace(0, th_end, 400).';
NL.eps   = p(NL.theta) / p0 - 1;
sig_of = @(th) 8 * k_r * (th - theta0) ./ (p0 * t * dp(th));
eps_of = @(th) p(th) / p0 - 1;
NL.sigma = sig_of(NL.theta);
NL.k_r = k_r;  NL.p0 = p0;  NL.b0 = b0;

% Tangent modulus exactly at the as-printed state (central difference)
dth = 1e-6;
NL.slope0 = (sig_of(theta0 + dth) - sig_of(theta0 - dth)) / ...
            (eps_of(theta0 + dth) - eps_of(theta0 - dth));
end
