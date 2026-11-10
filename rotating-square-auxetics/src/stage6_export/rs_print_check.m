function PC = rs_print_check(P, ps)
%RS_PRINT_CHECK  Can an FDM printer make this part? Prints a PASS/WARN list.
%
%   Checks (FDM, part printed FLAT: lattice in the XY plane, t along Z):
%   1. ligament width w >= min_perims x nozzle  (needs >= 2 extrusion lines
%      to be a real wall; 1 line is weak and often missing)
%   2. w is close to a whole number of extrusion widths (else the slicer
%      leaves a gap or an over-filled line in the hinge)
%   3. narrowest slot between squares >= min_gap_factor x nozzle
%      (narrower slots tend to fuse shut / string)
%   4. thickness t is a whole number of layers
%   5. the part is one connected piece, with the expected (N-1)^2 pores
%   Why flat: the hinge bends IN the XY plane, so its tensile/compressive
%   fibres run along the extruded lines, not across weak layer bonds.
g = P.geom;  pr = P.print;
th0 = deg2rad(g.theta0_deg);
items = {};

w_min = pr.min_perims * pr.nozzle;
items(end+1,:) = {'ligament width w', sprintf('%.2f mm (min %.2f)', g.w, w_min), g.w >= w_min - 1e-9};

nlines = g.w / pr.nozzle;
items(end+1,:) = {'w / nozzle (whole lines?)', sprintf('%.2f lines', nlines), ...
    abs(nlines - round(nlines)) < 0.15};

Ct = rs_contact_limit(g.a, g.Lh, th0);
gap = Ct.gap0;
items(end+1,:) = {'narrowest slot between squares', ...
    sprintf('%.2f mm (min %.2f)', gap, P.cad.min_gap_factor*pr.nozzle), ...
    gap >= P.cad.min_gap_factor*pr.nozzle};

nl = g.t / pr.layer_height;
items(end+1,:) = {'thickness / layer height', sprintf('%.2f layers', nl), ...
    abs(nl - round(nl)) < 1e-6};

nreg  = ps.NumRegions;
nhole = ps.NumHoles;
items(end+1,:) = {'one piece, (N-1)^2 pores', sprintf('%d piece(s), %d pores (expect %d)', ...
    nreg, nhole, (P.cad.N-1)^2), nreg == 1 && nhole == (P.cad.N-1)^2};

fprintf('\n  Print-feasibility check (%.1f mm nozzle, %.2f mm layers, print flat)\n', ...
    pr.nozzle, pr.layer_height);
for k = 1:size(items, 1)
    if items{k,3}, r = 'OK  '; else, r = 'WARN'; end
    fprintf('  [%s] %-32s %s\n', r, items{k,1}, items{k,2});
end
PC.items = items;  PC.all_ok = all([items{:,3}]);  PC.gap = gap;
end
