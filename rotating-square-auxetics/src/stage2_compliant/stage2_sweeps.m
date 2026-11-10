function SW = stage2_sweeps(P)
%STAGE2_SWEEPS  Vary one parameter at a time around the defaults.
E   = hinge_modulus(P);
g   = P.geom;
th0 = deg2rad(g.theta0_deg);
rf  = P.model.root_length_factor;
n   = P.stage2.n_sweep;

SW.w  = linspace(P.stage2.w_range(1),  P.stage2.w_range(2),  n).';
SW.Lh = linspace(P.stage2.Lh_range(1), P.stage2.Lh_range(2), n).';
SW.a  = linspace(P.stage2.a_range(1),  P.stage2.a_range(2),  n).';

SW.vs_w  = rs_compliant_props(g.a,  g.Lh,  SW.w, g.t, E, th0, rf);
SW.vs_Lh = rs_compliant_props(g.a,  SW.Lh, g.w,  g.t, E, th0, rf);
SW.vs_a  = rs_compliant_props(SW.a, g.Lh,  g.w,  g.t, E, th0, rf);

% nu vs w for a family of square sizes
fam = P.stage2.a_family;
SW.nu_family = zeros(n, numel(fam));
for k = 1:numel(fam)
    C = rs_compliant_props(fam(k), g.Lh, SW.w, g.t, E, th0, rf);
    SW.nu_family(:, k) = C.nu;
end
end
