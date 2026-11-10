function E = hinge_modulus(P)
%HINGE_MODULUS  Modulus used for ligament bending/stretching [MPa].
%   Applies the print knock-down and, if chosen, the wide-beam correction.
E = P.mat.E * P.mat.print_knockdown;
switch P.model.hinge_E_mode
    case 'plane_stress'
        % 2D plane-stress lattice: consistent with a 2D Ansys model.
    case 'wide_beam'
        % Ligament depth t >> width w: anticlastic bending is suppressed.
        E = E / (1 - P.mat.nu^2);
    otherwise
        error('Unknown P.model.hinge_E_mode "%s"', P.model.hinge_E_mode);
end
end
