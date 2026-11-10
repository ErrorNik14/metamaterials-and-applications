function x0 = first_crossing(x, y, level)
%FIRST_CROSSING  First x (moving away from x = 0) where y reaches level.
%   x must contain 0 (or values on both sides). Linear interpolation.
%   Returns [x_negative_side, x_positive_side]; NaN if never reached.
x = x(:); y = y(:);
[x, k] = sort(x);  y = y(k);
x0 = [NaN NaN];
i0 = find(x >= 0, 1);
% positive side
for i = i0:numel(x)-1
    if y(i) < level && y(i+1) >= level
        x0(2) = x(i) + (level - y(i)) * (x(i+1) - x(i)) / (y(i+1) - y(i));
        break
    end
end
% negative side
for i = i0:-1:2
    if y(i) < level && y(i-1) >= level
        x0(1) = x(i) + (level - y(i)) * (x(i-1) - x(i)) / (y(i-1) - y(i));
        break
    end
end
end
