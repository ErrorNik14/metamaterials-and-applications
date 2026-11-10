function [ps, parts] = rs_cad_polyshape(P)
%RS_CAD_POLYSHAPE  The printable 2D outline of the N x N part (as printed).
%
%   [ps, parts] = rs_cad_polyshape(P)
%   ps    : one polyshape = union of all squares and all ligament strips
%   parts : struct with the separate squares and strips (for checks)
%
%   Geometry (same as Stages 1-5, theta = theta0):
%   - squares of side a, centres on a grid of pitch p0, rotated -/+ theta0
%   - each ligament: a strip of width w along the line joining the two
%     hinge corners (length Lh corner-to-corner), extended lig_overlap*w
%     into each square so that it fuses with the squares into one solid.
%   The extension leaves small nubs where it pokes past a square's edge
%   near the corner; they are part of the real printed hinge root.
g   = P.geom;
th0 = deg2rad(g.theta0_deg);
N   = P.cad.N;
L0  = rs_layout(th0, g.a, g.Lh, th0, N, 2);       % straight ligaments at theta0

sq = polyshape.empty;
for k = 1:numel(L0.corners)
    sq(end+1) = polyshape(L0.corners{k}); %#ok<AGROW>
end
ext = P.cad.lig_overlap * g.w;
strips = polyshape.empty;
for k = 1:numel(L0.lig)
    A = L0.lig{k}(1, :);  B = L0.lig{k}(end, :);
    u = (B - A) / norm(B - A);  n = [-u(2) u(1)];
    A2 = A - ext*u;  B2 = B + ext*u;
    R = [A2 + n*g.w/2; B2 + n*g.w/2; B2 - n*g.w/2; A2 - n*g.w/2];
    strips(end+1) = polyshape(R); %#ok<AGROW>
end
ps = union([sq strips]);
c  = mean(L0.centre, 1);
ps = translate(ps, -c(1), -c(2));                  % part centred on (0,0)
parts.squares = sq;  parts.strips = strips;  parts.layout = L0;
end
