function K = rs_cell_kinematics(theta, a, theta0)
%RS_CELL_KINEMATICS  Closed-form kinematics of the ideal rotating-squares cell.
%
%   K = rs_cell_kinematics(theta, a, theta0)
%     theta  : rotation angle(s) of the squares [rad], vector allowed
%     a      : side of a square [mm]
%     theta0 : reference (as-printed) angle [rad]; strains are measured from it
%
%   MODEL: rigid squares, frictionless pins at the corners. Neighbouring
%   squares rotate by -theta and +theta (checkerboard). See
%   docs/stage1_derivation.md for the full derivation.
%
%   KEY RESULT: centre-to-centre pitch of neighbouring squares
%       p(theta) = a (cos theta + sin theta)      (same in x and in y)

theta = theta(:);
K.theta = theta;

% --- Pitch in each direction ---------------------------------------------
% x: hinge chain square -> right neighbour; y: square -> top neighbour.
% Both derivations (see docs) give the same expression, so the expansion is
% isotropic. They are written separately so the physics stays visible.
px  = a * (cos(theta)  + sin(theta));
py  = a * (cos(theta)  + sin(theta));
px0 = a * (cos(theta0) + sin(theta0));
py0 = a * (cos(theta0) + sin(theta0));

% Representative unit cell = 2 x 2 squares (one A-B-A-B checkerboard repeat)
K.Lx = 2 * px;
K.Ly = 2 * py;

% --- Stretch ratios and engineering strains relative to theta0 -----------
K.lambda_x = px / px0;
K.lambda_y = py / py0;
K.eps_x = K.lambda_x - 1;
K.eps_y = K.lambda_y - 1;

% --- Poisson's ratio -------------------------------------------------------
% Incremental (tangent) ratio: nu = -(d eps_y / d theta) / (d eps_x / d theta)
dpx = a * (cos(theta) - sin(theta));     % d px / d theta
dpy = a * (cos(theta) - sin(theta));     % d py / d theta
K.nu_incremental = -dpy ./ dpx;
K.nu_incremental(abs(dpx) < eps(a)) = NaN;   % 0/0 at the 45 deg lock

% Secant (total) ratio relative to theta0: nu = -eps_y / eps_x
K.nu_secant = -K.eps_y ./ K.eps_x;
K.nu_secant(abs(K.eps_x) < eps) = NaN;       % 0/0 at theta = theta0

% --- Areas -----------------------------------------------------------------
% A 2x2 cell holds 4 squares of area a^2 whatever theta is (rigid squares).
K.area_cell      = K.Lx .* K.Ly;             % = 4 a^2 (1 + sin 2theta)
K.area_solid     = 4 * a^2;
K.solid_fraction = K.area_solid ./ K.area_cell;   % = 1 / (1 + sin 2theta)
K.porosity       = 1 - K.solid_fraction;
K.area_ratio     = K.area_cell / (4 * px0 * py0);  % A / A0 = (1+eps)^2

% --- Locking limits ----------------------------------------------------------
% Opening: dp/dtheta = 0  ->  cos = sin  ->  theta = 45 deg (max pitch).
% Closing: squares touch edge-to-edge at theta = 0.
K.theta_lock_open   = pi/4;
K.theta_lock_closed = 0;
K.eps_lock_tension     = a*sqrt(2) / px0 - 1;   % strain from theta0 to 45 deg
K.eps_lock_compression = a          / px0 - 1;   % strain from theta0 to 0 deg
end
