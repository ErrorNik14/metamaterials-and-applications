function St = rs_state_at_theta(P, theta, edge_factor)
%RS_STATE_AT_THETA  Everything worth displaying at one square rotation theta.
%
%   St = rs_state_at_theta(P, theta)                theta in rad
%   St = rs_state_at_theta(P, theta, edge_factor)   (from Stage 4, default 1)
%
%   Uses the Stage 2 large-rotation path (rotation mode) and the Stage 4
%   hinge-stress formulas, so the slider app and the animation show the
%   same numbers as the analysis.
if nargin < 3, edge_factor = 1; end
g   = P.geom;  E = hinge_modulus(P);
th0 = deg2rad(g.theta0_deg);
rf  = P.model.root_length_factor;
L   = g.Lh + rf*g.w;
k_r = E * g.t * g.w^3 / 12 / L;
sy  = P.mat.sigma_y * P.mat.print_knockdown;

p   = @(th) g.a*(cos(th) + sin(th)) + g.Lh*cos(th - th0);
dp  = @(th) g.a*(cos(th) - sin(th)) - g.Lh*sin(th - th0);
p0  = p(th0);

St.theta_deg = rad2deg(theta);
St.eps       = p(theta) / p0 - 1;
St.porosity  = 1 - g.a^2 / p(theta)^2;          % squares / cell area (pitch^2)
C = rs_compliant_props(g.a, g.Lh, g.w, g.t, E, th0, rf);
St.nu_ideal     = -1;
St.nu_compliant = C.nu;                          % small-strain value
St.E_eff        = C.E_eff;
St.sigma = 8 * k_r * (theta - th0) / (p0 * g.t * dp(theta));     % lattice stress
St.sig_hinge = edge_factor * P.fail.Kt * ...
    (E * g.w * abs(theta - th0) / L + abs(St.sigma) * p0 / g.w);  % worst ligament

% traffic light, same thresholds as the Stage 4 safe range
if St.sig_hinge <= sy / P.fail.safety_factor
    St.status = 'SAFE';          St.rgb = P.plot.c_good;
elseif St.sig_hinge <= sy
    St.status = 'BEYOND SF';     St.rgb = P.plot.c_warn;
else
    St.status = 'YIELDED';       St.rgb = P.plot.c_lock;
end
% Above yield the elastic formula no longer gives a real stress
if St.sig_hinge <= sy
    St.sig_text = sprintf('%.0f MPa', St.sig_hinge);
else
    St.sig_text = sprintf('> yield (elastic est. %.0f MPa)', St.sig_hinge);
end
end
