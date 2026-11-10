function [s_max, s_centre] = fe_ligament_stress(M, U, Kt)
%FE_LIGAMENT_STRESS  Peak outer-fibre stress in the FE ligaments.
%   sigma = Kt * ( |N|/A + |M| (w/2) / I ) at both ends of every element.
%   s_max:    worst ligament anywhere in the array (often near an edge)
%   s_centre: worst ligament attached to the centre square (the closest
%             thing to the infinite lattice that the theory describes)
ne = size(M.conn, 1);
s_el = zeros(ne, 1);
for e = find(M.is_lig).'
    n1 = M.conn(e,1);  n2 = M.conn(e,2);
    dofs = [3*n1-2, 3*n1-1, 3*n1, 3*n2-2, 3*n2-1, 3*n2];
    [~, T, kl] = fe_frame_k(M.E(e), M.A(e), M.I(e), M.X(n1,:), M.X(n2,:));
    f = kl * (T * U(dofs));
    c = sqrt(M.I(e) * 12 / M.A(e)) / 2;          % half the width, w/2
    s_el(e) = Kt * (abs(f(1))/M.A(e) + max(abs(f([3 6]))) * c / M.I(e));
end
s_max = max(s_el);

% ligaments attached to the centre square
N = M.N;  jc = ceil(N/2);
centre_nodes = squeeze(M.corner(jc, jc, :));
s_centre = 0;
for k = 1:numel(M.lig_elems)
    el = M.lig_elems{k};
    if any(ismember(M.conn(el([1 end]), :), centre_nodes), 'all')
        s_centre = max(s_centre, max(s_el(el)));
    end
end
end
