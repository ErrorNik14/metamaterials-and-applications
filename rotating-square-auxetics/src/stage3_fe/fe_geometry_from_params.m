function g = fe_geometry_from_params(P, w)
%FE_GEOMETRY_FROM_PARAMS  Pack the inputs rs_fe_mesh needs (w optional override).
if nargin < 2, w = P.geom.w; end
g.a  = P.geom.a;
g.Lh = P.geom.Lh;
g.w  = w;
g.t  = P.geom.t;
g.theta0 = deg2rad(P.geom.theta0_deg);
g.E  = hinge_modulus(P);
g.n_el_lig = P.fe.n_el_lig;
g.square_E_factor = P.fe.square_E_factor;
end
