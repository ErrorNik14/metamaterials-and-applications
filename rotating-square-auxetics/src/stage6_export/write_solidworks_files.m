function files = write_solidworks_files(P, S2, folder)
%WRITE_SOLIDWORKS_FILES  Parameters for SolidWorks (CSV + equations file)
%   and the square centre/rotation table.
%
%   1) S6_solidworks_parameters.csv  : name, value, unit, description
%   2) S6_solidworks_equations.txt   : SolidWorks "Link to external file"
%      format, one  "name" = value  per line (Tools > Equations > Link to
%      external file). Lengths in mm, angles in degrees.
%   3) S6_square_positions.csv       : i, j, x, y, rotation of every square
g   = P.geom;  C = S2.C0;
th0 = deg2rad(g.theta0_deg);
names = {'a'; 'Lh'; 'w'; 't'; 'theta0'; 'pitch'; 'cell'; 'N'; 'part_size'; ...
         'lig_overlap'};
vals  = [g.a; g.Lh; g.w; g.t; g.theta0_deg; C.p0; 2*C.p0; P.cad.N; ...
         (P.cad.N - 1)*C.p0 + g.a*(cos(th0) + sin(th0)); P.cad.lig_overlap*g.w];
units = {'mm'; 'mm'; 'mm'; 'mm'; 'deg'; 'mm'; 'mm'; '-'; 'mm'; 'mm'};
desc  = {'square side'; 'ligament length (corner to corner)'; 'ligament width'; ...
         'extrusion thickness'; 'as-printed square rotation (A: -theta0, B: +theta0)'; ...
         'centre-to-centre pitch = a(cos th0 + sin th0) + Lh'; ...
         '2x2 repeat cell = 2*pitch (use for linear pattern)'; 'squares per side'; ...
         'overall part size, corner to corner of the rotated squares'; ...
         'ligament runs this far into each square'};
T = table(names, vals, units, desc, 'VariableNames', {'name','value','unit','description'});
files.csv = fullfile(folder, 'S6_solidworks_parameters.csv');
writetable(T, files.csv);

files.eq = fullfile(folder, 'S6_solidworks_equations.txt');
fid = fopen(files.eq, 'w');
for k = 1:numel(names)
    if strcmp(units{k}, '-')
        fprintf(fid, '"%s" = %.6g\n', names{k}, vals(k));
    else
        fprintf(fid, '"%s" = %.6g%s\n', names{k}, vals(k), units{k});
    end
end
fclose(fid);

% square table (as printed, part centred on the origin)
L0 = rs_layout(th0, g.a, g.Lh, th0, P.cad.N, 2);
c  = mean(L0.centre, 1);
[J, I] = ndgrid(1:P.cad.N, 1:P.cad.N);
typ = 1 + mod(I + J, 2);
rotd = g.theta0_deg * (2*typ - 3);                % A: -theta0, B: +theta0
Tq = table(I(:), J(:), L0.centre(:,1) - c(1), L0.centre(:,2) - c(2), rotd(:), ...
    'VariableNames', {'i_column', 'j_row', 'x_mm', 'y_mm', 'rotation_deg'});
files.squares = fullfile(folder, 'S6_square_positions.csv');
writetable(Tq, files.squares);
end
