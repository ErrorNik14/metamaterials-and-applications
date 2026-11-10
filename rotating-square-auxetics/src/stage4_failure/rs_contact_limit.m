function Ct = rs_contact_limit(a, Lh, theta0)
%RS_CONTACT_LIMIT  Compression lock: strain at which two squares first touch.
%
%   Squares (side a, rotation -/+ theta) sit on a square grid of pitch
%       p(theta) = a (cos th + sin th) + Lh cos(th - th0)     (Stage 2)
%   Under compression theta decreases from theta0 (it may go negative:
%   the pinwheel then turns the other way). For each theta we compute the
%   smallest gap between the centre square and its 8 neighbours with the
%   separating-axis test for convex polygons (gap < 0 means overlap).
%   Contact = first theta below theta0 where the gap reaches 0.
%   (The ligament itself is ignored - only square-to-square contact.)
p   = @(th) a*(cos(th) + sin(th)) + Lh*cos(th - theta0);
gap = @(th) min_gap(th, a, p(th));

Ct.gap0 = gap(theta0);
th_lo = -pi/4;                               % scan far enough
ths = linspace(theta0, th_lo, 2000);
g = arrayfun(gap, ths);
k = find(g <= 0, 1);
if isempty(k)
    Ct.theta = NaN;  Ct.eps = NaN;
else
    Ct.theta = fzero(gap, [ths(k), ths(k-1)]);
    Ct.eps   = p(Ct.theta) / p(theta0) - 1;
end
Ct.theta_scan = ths;  Ct.gap_scan = g;
end

function g = min_gap(th, a, pitch)
rot = @(al) [cos(al) -sin(al); sin(al) cos(al)];
S = (a/2) * [1 1; -1 1; -1 -1; 1 -1];
c0 = (rot(-th) * S.').';                     % centre square (type A)
g = inf;
for di = -1:1
    for dj = -1:1
        if di == 0 && dj == 0, continue; end
        s = (-1)^(di + dj);                   % neighbour type: B if odd step
        nb = (rot(s * -th) * S.').' + pitch * [di dj];
        g = min(g, sat_gap(c0, nb));
    end
end
end

function g = sat_gap(A, B)
% Separating-axis distance between convex polygons (lower bound of the
% true distance when separated; negative = penetration depth).
g = -inf;
for P = {A, B}
    V = P{1};
    for k = 1:size(V, 1)
        e = V(mod(k, size(V,1)) + 1, :) - V(k, :);
        n = [-e(2) e(1)] / norm(e);
        pa = A * n.';  pb = B * n.';
        g = max(g, max(min(pb) - max(pa), min(pa) - max(pb)));
    end
end
end
