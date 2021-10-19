function xq = lanczos(t, x, tq, a)
% xq = LANCZOS(t, x, tq, a)
% 
% Interpolate a signal using Lanczos kernel.
%
% INPUT:
% t         time of the original signal, must be equally spaced
% x         original signal
% tq        requested time for interpolated signal
% a         scale factor for Lanczos kernel
%           If a = 0, then the Whittaker-Shannon interpolation is used
%           instead.
%
% OUTPUT:
% xq        interpolated signal
%
% EXAMPLE:
% % run a demo
% lanczos('demo');
%
% Last modified by sirawich-at-princeton.edu, 10/19/2021

defval('a', 0)

% demo
if ischar(t)
    % original signal
    t = (-10:0.227:10)';
    A = [6 5 4]';
    f = [1 0.75 0.3]';
    p = [0.1 2/pi 0.4]';
    x = sin(2 * pi * (t * f' + p')) * A;
    % requested signal
    tq = (-15:0.01:15)';
    xq = lanczos(t, x, tq, 3);
    xq_expected = sin(2 * pi * (tq * f' + p')) * A;
    % plot the result
    figure
    subplot(3,1,1)
    plot(t, x, 'k-o')
    hold on
    grid on
    plot(tq, xq, 'r')
    legend('original', 'reconstruct')
    ylabel('x')
    title(sprintf(['lanzcos demo: dt-original = %.3f s, ' ...
        'dt-interpolated = %.3f s, a = %.2f'], t(2)-t(1), tq(2)-tq(1), 3));
    
    subplot(3,1,2)
    plot(tq, xq_expected, 'k');
    hold on
    grid on
    plot(tq, xq, 'r')
    legend('expected', 'reconstruct')
    ylabel('x')
    
    subplot(3,1,3)
    plot(tq, xq - xq_expected)
    grid on
    ylabel('difference')
    xlabel('t')
    return
end

if size(t, 2) > 1
    t = t';
end
if size(x, 2) > 1
    x = x';
end
if size(tq, 2) > 1
    tq = tq';
end

% sampling period
dt = t(2) - t(1);

if a == 0
    % Whittaker-Shannon interpolation formula
    xq = sinc((repmat(tq, 1, length(t)) - repmat(t', length(tq), 1)) / dt) * x;
else
    % apply Lanczos window
    xq = lanczoswindow((repmat(tq, 1, length(t)) - repmat(t', length(tq), 1)) / dt, a) * x;
end
end

function y = lanczoswindow(x, a)
y = sinc(x) .* sinc(x/a);
y(abs(x) > a) = 0;
end