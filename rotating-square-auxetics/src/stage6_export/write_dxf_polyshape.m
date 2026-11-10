function write_dxf_polyshape(ps, file)
%WRITE_DXF_POLYSHAPE  Write every boundary of a polyshape as DXF LINEs (mm).
%   Minimal ASCII DXF (R12 / AC1009): opens in SolidWorks, Fusion,
%   Inkscape, LibreCAD. In SolidWorks: File > Open > .dxf > "Import to a
%   new part as 2D sketch", units = millimetres. Then Extruded Boss by t.
%   Outer boundary and holes are all closed loops, so one extrude makes
%   the whole lattice.
fid = fopen(file, 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, '0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1009\n0\nENDSEC\n');
fprintf(fid, '0\nSECTION\n2\nENTITIES\n');
for b = 1:numboundaries(ps)
    [x, y] = boundary(ps, b);           % closed: last point = first point
    for k = 1:numel(x) - 1
        fprintf(fid, '0\nLINE\n8\nLATTICE\n10\n%.6f\n20\n%.6f\n30\n0.0\n11\n%.6f\n21\n%.6f\n31\n0.0\n', ...
            x(k), y(k), x(k+1), y(k+1));
    end
end
fprintf(fid, '0\nENDSEC\n0\nEOF\n');
end
