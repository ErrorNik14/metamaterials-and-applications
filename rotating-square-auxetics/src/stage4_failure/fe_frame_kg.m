function kgG = fe_frame_kg(Nax, x1, x2)
%FE_FRAME_KG  Geometric (stress) stiffness of a 2D frame element, global axes.
%   Nax = axial force in the element (tension +) [N].
%   Consistent geometric stiffness of the cubic beam (local axes):
%     kG = N/L * [0   0     0      0   0     0
%                 0   6/5   L/10   0  -6/5   L/10
%                 0   L/10  2L^2/15 0 -L/10 -L^2/30
%                 0   0     0      0   0     0
%                 0  -6/5  -L/10   0   6/5  -L/10
%                 0   L/10 -L^2/30 0  -L/10  2L^2/15]
%   Compression (N < 0) makes the element softer to sideways motion:
%   buckling is when K + lambda*KG becomes singular.
d = x2(:) - x1(:);  L = norm(d);
c = d(1)/L;  s = d(2)/L;
k = zeros(6);
k([2 3 5 6], [2 3 5 6]) = [ 6/5    L/10    -6/5   L/10;
                            L/10   2*L^2/15 -L/10 -L^2/30;
                           -6/5   -L/10     6/5  -L/10;
                            L/10  -L^2/30  -L/10  2*L^2/15];
k = Nax / L * k;
R = [c s 0; -s c 0; 0 0 1];
T = blkdiag(R, R);
kgG = T.' * k * T;
end
