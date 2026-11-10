function theta = rs_theta_from_strain(eps_x, theta0)
%RS_THETA_FROM_STRAIN  Inverse kinematics: square rotation for an applied strain.
%
%   theta = rs_theta_from_strain(eps_x, theta0)      [rad]
%
%   From  (1 + eps) = (cos th + sin th) / (cos th0 + sin th0)
%   and   cos th + sin th = sqrt(2) sin(th + pi/4):
%       theta = asin( (1+eps)(cos th0 + sin th0) / sqrt(2) ) - pi/4
%   Valid between the two locks (0 <= theta <= 45 deg). Strains outside
%   that range return NaN (the ideal mechanism cannot reach them).
s = (1 + eps_x) .* (cos(theta0) + sin(theta0)) / sqrt(2);
theta = asin(min(s, 1)) - pi/4;
theta(s > 1 | theta < -eps) = NaN;
end
