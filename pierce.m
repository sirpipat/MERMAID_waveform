function p = pierce(m, x0, x, method, precision)
% p = PIERCE(m, [x0 y0], [xi yi], method)
%
% Determine the location of piercing points of a line on an interface. The
% line is defined by the slope M and a point [X0 Y0] on the line. The
% interface is defined by a set of points [Xi Yi].
%
% INPUT:
% m             slope of the straight line
% [x0 y0]       xy-coordinate of a point on a line
% [xi yi]       xy-coordinate of points on an interface
% method        interpolation method    [Default: 'linear']
% precision     precision of the output [Default: 0.01]
%
% OUTPUT:
% p             xy-coordinate of the piercing point(s)
%
% Last modified by sirawich-at-princeton.edu: 04/04/2023

defval('method', 'linear')
defval('precision', 0.01)

% (v,w) coordinate of the interface minus the pierce line
v = x(:,1);
w = x(:,2) - y_line(m, x0, x(:,1));

% find the sections where the line intercepts on the horizontal axis
where = (diff(sign(w)) ~= 0);
where_begin = v([where; false]);
where_end = v([false; where]);

p = nan(length(where_begin), 2);
for ii = 1:length(where_begin)
    vmin = floor(where_begin(ii) / precision) * precision;
    vmax = ceil(where_end(ii) / precision) * precision;
    % use bracketing method to speed up when the percision is very small
    while (vmax - vmin) / precision > 100
        vq = linspace(vmin, vmax, 101)';
        wq = interp1(v, w, vq, method);
        where = (diff(sign(wq)) ~= 0);
        vmin = vq([where; false]);
        vmax = vq([false; where]);
        vmin = floor(vmin / precision) * precision;
        vmax = ceil(vmax / precision) * precision;
    end
    vq = (vmin:precision:vmax)';
    wq = interp1(v, w, vq, method);
    [wq0, jj] = min(abs(wq));
    p(ii,1) = vq(jj);
    p(ii,2) = wq0 + y_line(m, x0, vq(jj));
end
end

% compute the y-coordinates of points on the line given the x-coordinates
function y = y_line(m, x0, x)
y = m * (x - x0(1)) + x0(2);
end