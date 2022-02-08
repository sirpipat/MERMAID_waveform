function xs = removenoise(t, x, fs, arrival, n)
% xs = REMOVENOISE(t, x, fs, arrival, n)
%
% Determines the Fourier coefficients of the noise time-series before 
% ARRIVAL which separate the noise-only section from the other section and 
% subtracts the time-series generated from those coefficients from the
% signal.
%
% INPUT:
% t             time array or datetime array
% x             input signal
% fs            rate of sampling
% arrival       when the signal arrives, either time or datetime, must be
%               consistent with t
% n             number of frequencies to fit
%
% OUTPUT:
% xs            output signal with noise removed
%
% Last modified by sirawich-at-princeton.edu. 02/08/2022

if ischar(t) && strcmpi(t, 'demo')
    defval('x', 0.5)
    fs = 40;
    t = 0:(1/fs):10;
    arrival = 6;
    noiseamplitude = x;
    
    % original signal
    xnoise1 = 0.5*sin(2*pi*t);
    xnoise2 = 0.125*cos(2*pi*2*t+0.03);
    xnoise3 = 0.25*sin(2*pi*4*t+2);
    xnoise = xnoise1 + xnoise2 + xnoise3;
    xrandom = noiseamplitude*(2*rand(size(t)) - 1);
    xsignal = 8*exp(-(t-8).^2 *2).*cos(2*pi*4*t+1).*sin(2*pi/10*(t-8));
    x = xsignal + xnoise + xrandom;
    % bandpass [3 10] Hz for comparision
    xf = bandpass(x, fs, 3, 5, 2, 2, 'butter', 'linear');
    
    % remove the noise
    n = [1 2 3 4 6 10 20 40 60 100];
    xs = zeros(length(n), length(x));
    
    for ii = 1:length(n)
        xs(ii,:) = removenoise(t, x, fs, arrival, n(ii));
    end
    
    figure
    clf
    set(gcf, 'Unit', 'inches', 'Position', [12 6 12 9])
    ax = subplot('Position', [0.16 0.08 0.8 0.9]);
    
    % plot original signal
    plot(t, x, 'LineWidth', 1, 'Color', 'k')
    hold on
    % plot the denoised signals
    for ii = 1:length(n)
        plot(t, xs(ii,:) - 3*ii, 'Color', [0.8 0 0], 'LineWidth', 1)
    end
    axes(ax);
    % plot bandpass for comparison
    plot(t, xf - 3*(length(n)+1), 'LineWidth', 1, 'Color', [0.4 0.7 0.2]);
    % plot the interesting signal
    plot(t, xsignal - 3*(length(n)+2), 'LineWidth', 1, 'Color', 'b')
    plot(t, xnoise1 - 3*(length(n)+3), 'LineWidth', 1, 'Color', 'b')
    plot(t, xnoise2 - 3*(length(n)+4), 'LineWidth', 1, 'Color', 'b')
    plot(t, xnoise3 - 3*(length(n)+5), 'LineWidth', 1, 'Color', 'b')
    plot(t, xrandom - 3*(length(n)+6), 'LineWidth', 1, 'Color', 'b')
    grid on
    ylim(3 * ([-ii-6.5 0.8]))
    set(gca, 'FontSize', 12, 'TickDir', 'both')
    yticks(-48:3:0)
    yticklabels({'noise: random', 'noise: 4 Hz wave', ...
        'noise: 2 Hz wave', 'noise: 1 Hz wave', 'signal', ...
        'bandpass 3--5 Hz', ...
        'n = 100', 'n = 60', 'n = 40', 'n = 20', 'n = 10', 'n = 6', ...
        'n = 4', 'n = 3', 'n = 2', 'n = 1', 'original'})
    vline(gca, arrival, 'LineWidth', 2, 'Color', [0.6 0.6 0.6], ...
        'LineStyle', '-');
    xlabel('time (s)')
    figdisp('removenoise_demo.eps', [], [], 2, [], 'epstopdf');
    return
end

% convert all inputs to row vectors
if size(t, 1) > 1
    t = t';
end
if size(x, 1) > 1
    x = x';
    converted = true;
else
    converted = false;
end

% if datetimes are given, convert to doubles representing the number of
% seconds from the first sample
if isdatetime(t) && isdatetime(arrival)
    arrival = seconds(arrival - t(1));
    t = seconds(t - t(1));
elseif xor(isdatetime(t), isdatetime(arrival))
    error('T and ARRIVAL must be the same data type: either time or datetime.')
end

% slice for section to determine the noise
xb = x(t <= arrival);
tb = t(t <= arrival);

% determine target frequencies
max_f = fs/2;
min_f = 1/(tb(end)-tb(1));

f = min_f:min_f:(max_f-min_f);
max_n = length(f);
n = min(n, max_n);

% check whether sinefit throws a warning about singular matrix
lastwarn('', '');
[A, B, C, ~, F, P] = sinefit(tb, xb, [], f);
% now if a warning was raised, warnmsg and warnid will not be empty
[warnmsg, warnid] = lastwarn();
if ~isempty(warnid)
    error(warnmsg, warnid);
end

A = [0 sqrt(2)*A];
B = [C sqrt(2)*B];

% find k most significant frequency
[~, I] = maxk(P, n);

A = A(I);
B = B(I);
F = F(I);

% constrct the noise signal
xn = A * sin(2 * pi * F' * t) + B * cos(2 * pi * F' * t);

% remove the noise
xs = x - xn;

% do not forget to convert back
if converted
    xs = xs';
end
end