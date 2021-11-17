function [r, dd, n, val] = spectraldivision(d, u, w, reg, val)
% r = SPECTRALDIVISION(d, u, w, reg, val)
%
% compute r where d = conv(r, u) via spectral division
% See Pesce+2010 for method details
%
% INPUT:
% d         convolved function
% u         unconvolved function
% w         window function
% reg       regularization method (either 'damp' or 'water')
% val       regularization value
%
% OUTPUT:
% r         result time-series with the length of the longest input signal
% dd        convolution of 'r' and 'u' with the same length as 'd'
% n         residue norm between 'd' and 'dd'
% val       regularization value
%
% Last modified by sirawich-at-princeton.edu, 11/17/2021


% run the demo
if strcmpi(reg, 'demo')
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 6 6 9])
    
    ax1 = subplot('Position', [0.08 0.86 0.86 0.12]);
    [r, ~, n] = spectraldivision(d, u, w, 'damp', 0);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('damping factor = %0.5g, norm = %0.5g', 0, n))
    set(ax1, 'Box', 'on', 'TickDir', 'both');

    ax2 = subplot('Position', [0.08 0.66 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 0.01;
    [r, ~, n] = spectraldivision(d, u, w, 'damp', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('damping factor = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax2, 'Box', 'on', 'TickDir', 'both');

    ax3 = subplot('Position', [0.08 0.46 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 0.1;
    [r, ~, n] = spectraldivision(d, u, w, 'damp', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('damping factor = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax3, 'Box', 'on', 'TickDir', 'both');
    
    ax4 = subplot('Position', [0.08 0.26 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 1;
    [r, ~, n] = spectraldivision(d, u, w, 'damp', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('damping factor = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax4, 'Box', 'on', 'TickDir', 'both');

    ax5 = subplot('Position', [0.08 0.06 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 10;
    [r, dd, n] = spectraldivision(d, u, w, 'damp', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('damping factor = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax5, 'Box', 'on', 'TickDir', 'both');
    
    ax3.YLim = ax2.YLim;
    ax4.YLim = ax2.YLim;
    ax5.YLim = ax2.YLim;

    % save figure
    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_demo_damping.eps', mfilename);
    figdisp(savename, [], [], 2, [], 'epstopdf');
    
    figure(2)
    clf
    set(gcf, 'Units', 'inches', 'Position', [6 6 6 9])
    
    ax1 = subplot('Position', [0.08 0.86 0.86 0.12]);
    [r, ~, n] = spectraldivision(d, u, w, 'water', 0);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('water level = %0.5g, norm = %0.5g', 0, n))
    set(ax1, 'Box', 'on', 'TickDir', 'both');

    ax2 = subplot('Position', [0.08 0.66 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 1;
    [r, ~, n] = spectraldivision(d, u, w, 'water', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('water level = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax2, 'Box', 'on', 'TickDir', 'both');

    ax3 = subplot('Position', [0.08 0.46 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 10;
    [r, ~, n] = spectraldivision(d, u, w, 'water', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('water level = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax3, 'Box', 'on', 'TickDir', 'both');
    
    ax4 = subplot('Position', [0.08 0.26 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 100;
    [r, ~, n] = spectraldivision(d, u, w, 'water', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('water level = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax4, 'Box', 'on', 'TickDir', 'both');

    ax5 = subplot('Position', [0.08 0.06 0.86 0.12]);
    asq =  max(abs(u)) ^ 2;
    val = asq * 1000;
    [r, dd, n] = spectraldivision(d, u, w, 'water', val);
    plot(r, 'Color', 'k')
    xlim([1 length(r)])
    grid on
    title(sprintf('water level = %0.5g x input amplitude^2, norm = %0.5g', ...
        val/asq, n))
    set(ax5, 'Box', 'on', 'TickDir', 'both');
    
    ax3.YLim = ax2.YLim;
    ax4.YLim = ax2.YLim;
    ax5.YLim = ax2.YLim;
    
    % save figure
    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_demo_water.eps', mfilename);
    figdisp(savename, [], [], 2, [], 'epstopdf');
    return
end

% if the value is not specified, then figure out the value that gives the
% minimum norm.
if isempty(val)
    vals = 10.^(-4:0.1:4) * max(abs(u)) ^ 2;
    ns = zeros(size(vals));
    for ii = 1:length(vals)
        [~, ~, ns(ii)] = spectraldivision(d, u, w, reg, vals(ii));
    end
    [~, imin] = min(ns);
    val = vals(imin);
end

% convert input signals to a column vector
if size(d, 2) > 1
    d = d';
end
if size(u, 2) > 1
    u = u';
end

% pad zeros to obtain the same length
if size(d, 1) > size(u, 1)
    u = [u; zeros(size(d, 1) - size(u, 1), 1)];
elseif size(d, 1) < size(u, 1)
    d = [d; zeros(size(u, 1) - size(d, 1), 1)];
end

% setting a window to have the same dimension of the signals
defval('w', ones(size(u)))
if size(w, 2) > 1
    w = w';
end
if size(w, 1) > size(u, 1)
    w = w(1:size(u, 1), 1);
elseif size(w, 1) < size(u, 1)
    w = [w; zeros(size(u, 1) - size(w, 1), 1)];
end

% number of frequency is set to a power of 2 to speed up FFT and IFFT
nf = 2 ^ nextpow2(length(u));

% convert to frequency domain
D = fft(d .* w, nf);
U = fft(u .* w, nf);

if strcmpi(reg, 'damp')
    RR = (D .* conj(U)) ./ (U .* conj(U) + val);
elseif strcmpi(reg, 'water')
    RR = (D .* conj(U)) ./ max(U .* conj(U), val);
else
    error('invalid method: reg must be eigther ''damp'' or ''water''\n')
end

% convert back to time domain
r = ifft(RR, nf);
r = r(1:size(u, 1), 1);

% expected output signal from 
dd = conv(u, r);
dd = dd(1:size(u, 1), 1);

% norm of residue
n = norm(d - dd, 2);
end