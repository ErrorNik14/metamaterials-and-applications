function H = rs_hinge_stress(a, Lh, w, t, E, theta0, rf, Kt, eps, NL)
%RS_HINGE_STRESS  Peak stress in a horizontal ligament vs applied strain eps_x.
%
%   H = rs_hinge_stress(a, Lh, w, t, E, theta0, rf, Kt, eps)       small strain
%   H = rs_hinge_stress(..., eps, NL)   also along the large-rotation path NL
%                                        from rs_nonlinear_response
%
%   Two stress sources in each ligament (outer fibre, beam theory):
%   1) BENDING: constant moment M = EI (2 dtheta)/L   (Stage 2, pure bending)
%        sigma_b = M (w/2) / I = E w dtheta / L
%   2) AXIAL: by equilibrium of a vertical cut through the lattice, each
%      horizontal ligament carries N = sigma_x * p0 * t, so
%        sigma_a = N / (w t) = sigma_x p0 / w
%   Peak (outer fibre, worst side):  sigma_max = Kt * (sigma_b + |sigma_a|)
%
%   Small strain: square rotation per unit strain from Stage 2 equilibrium,
%     dtheta/eps_x = k_a b p0 (1 - nu) / (2 k_a b^2 + 8 k_r)

C = rs_compliant_props(a, Lh, w, t, E, theta0, rf);
L = Lh + rf * w;

H.dtheta_per_eps = C.k_a * C.b * C.p0 * (1 - C.nu) / (2*C.k_a*C.b^2 + 8*C.k_r);
H.eps   = eps(:);
dth     = H.dtheta_per_eps * H.eps;
H.sig_b = Kt * E * w * abs(dth) / L;
H.sig_a = Kt * abs(C.E_eff * H.eps) * C.p0 / w;
H.sig   = H.sig_b + H.sig_a;
H.slope = H.sig(end) / abs(H.eps(end));            % MPa per unit strain (linear)
H.C = C;

if nargin >= 10 && ~isempty(NL)
    H.NL.eps   = NL.eps;
    H.NL.sig_b = Kt * E * w * abs(NL.theta - theta0) / L;
    H.NL.sig_a = Kt * abs(NL.sigma) * NL.p0 / w;
    H.NL.sig   = H.NL.sig_b + H.NL.sig_a;
end
end
