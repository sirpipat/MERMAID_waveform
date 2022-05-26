function [fc, s] = freqselect(t, x, fs, plt, titlename, savename)
% [fc, s] = FREQSELECT(t, x, fs, plt, titlename, savename)
%
% Figures out the frequency band where the signal stands out the most from
% the background noise.
% 
% INPUT:
% t             time
% x             time-series data
% fs            sampling rate
% plt           whether to plot or not
% titlename     name to put on the title
% savename      filename for the saved figure
%
% OUTPUT:
% fc            best corner frequency
% s             best signal-to-noise ratio
%
% Last modified by sirawich-at-princeton.edu, 05/26/2022

% list of corner frequency candidates
fcs = 0:0.05:2;

% bandpass/lowpass SNR
A = NaN(length(fcs), length(fcs));
T = NaN(length(fcs), length(fcs));

% bandstop/highpass SNR
B = NaN(length(fcs), length(fcs));
U = NaN(length(fcs), length(fcs));

for ii = 1:length(fcs)
    for jj = (ii+1):length(fcs)
        % for zero lower frequency: lowpass
        if fcs(ii) == 0
            xf = lowpass(x, fs, fcs(jj), 2, 2, 'butter', 'linear');
        elseif fcs(jj) / fcs(ii) >= 1.2
            xf = bandpass(x, fs, fcs(ii), fcs(jj), 2, 2, 'butter', 'linear');
        else
            continue
        end
        xs = x - xf;
        halfwin = 2/sqrt((fcs(ii)+fcs(jj))/2);
        [A(jj, ii), T(jj, ii)] = snrvar(t, xf, [-1 1] * halfwin, ...
            -60, 80);
        [B(jj, ii), U(jj, ii)] = snrvar(t, xs, [-1 1] * halfwin, ...
            -60, 80);
    end
end

% pick the best corner frequencies
[M, I] = max(A./B, [], 1);
[MM, J] = max(M);
fc = [fcs(J) fcs(I(J))];

% best SNR
s = A(I(J), J);

% visualize the result
if plt
    % locate the best pixel on the plot
    [xx, yy] = boxcorner(fc(1) + 0.025 * [-1 1], fc(2) + 0.025 * [-1 1]);
    
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [16 12 8 9])
    
    % title
    ax0 = subplot('Position', [0.08 0.95 0.88 0.01]);
    title(titlename)
    set(ax0, 'FontSize', 12, 'Color', 'none')
    ax0.XAxis.Visible = 'off';
    ax0.YAxis.Visible = 'off';
    
    % bandpass SNR
    ax1 = subplot('Position', [0.08 0.64 0.37 0.26]);
    im1 = imagesc(fcs, fcs, log10(A));
    axis xy
    xlabel('lower corner frequency (Hz)')
    ylabel('upper corner frequency (Hz)')
    title('bandpass')
    setimagenan(ax1, im1, [1 1 1]);
    c1 = colorbar(ax1, 'EastOutside');
    [~,v1] = vline(ax1, ax1.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    [~,h1] = hline(ax1, ax1.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    plot(ax1, xx, yy, '-r', 'LineWidth', 1)
    c1.Label.String = 'log_{10} SNR';
    set(ax1, 'FontSize', 12, 'TickDir', 'out')
    set(c1, 'FontSize', 12, 'TickDir', 'out')
    
    ax2 = subplot('Position', [0.58 0.64 0.37 0.26]);
    im2 = imagesc(fcs, fcs, log10(B));
    axis xy
    xlabel('lower corner frequency (Hz)')
    ylabel('upper corner frequency (Hz)')
    title('bandstop')
    setimagenan(ax2, im2, [1 1 1]);
    c2 = colorbar(ax2, 'EastOutside');
    [~,v2] = vline(ax2, ax2.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    [~,h2] = hline(ax2, ax2.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    plot(ax2, xx, yy, '-r', 'LineWidth', 1)
    c2.Label.String = 'log_{10} SNR';
    set(ax2, 'FontSize', 12, 'TickDir', 'out')
    set(c2, 'FontSize', 12, 'TickDir', 'out')
    
    ax3 = subplot('Position', [0.08 0.28 0.37 0.26]);
    im3 = imagesc(fcs, fcs, log10(A./B));
    axis xy
    xlabel('lower corner frequency (Hz)')
    ylabel('upper corner frequency (Hz)')
    title('ratio')
    setimagenan(ax3, im3, [1 1 1]);
    c3 = colorbar(ax3, 'EastOutside');
    [~,v3] = vline(ax3, ax3.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    [~,h3] = hline(ax3, ax3.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    plot(ax3, xx, yy, '-r', 'LineWidth', 1)
    c3.Label.String = 'log_{10} SNR';
    set(ax3, 'FontSize', 12, 'TickDir', 'out')
    set(c3, 'FontSize', 12, 'TickDir', 'out')
    
    ax4 = subplot('Position', [0.58 0.28 0.37 0.26]);
    im4 = imagesc(fcs, fcs, T);
    axis xy
    xlabel('lower corner frequency (Hz)')
    ylabel('upper corner frequency (Hz)')
    title('optimal time for bandpass')
    c4 = colorbar(ax4, 'EastOutside');
    colormap(ax4, kelicol)
    setimagenan(ax4, im4, [1 1 1]);
    ax4.CLim = [-1 1] * max(abs(T), [], 'all');
    [~,v4] = vline(ax4, ax4.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    [~,h4] = hline(ax4, ax4.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    plot(ax4, xx, yy, '-r', 'LineWidth', 1)
    c4.Label.String = 'time (s)';
    set(ax4, 'FontSize', 12, 'TickDir', 'out')
    set(c4, 'FontSize', 12, 'TickDir', 'out')
    
    ax5 = subplot('Position', [0.08 0.08 0.88 0.13]);
    xf = bandpass(x, fs, fc(1), fc(2), 2, 2, 'butter', 'linear');
    plot(ax5, t, xf, 'Color', 'k', 'LineWidth', 1)
    xlim(10 * round([-30 30] / 10 * 1/fc(1)))
    grid on
    [~, v5] = vline(ax5, T(I(J), J), 'LineWidth', 1, 'Color', [0 0.6 1]);
    xlabel('time since picked first arrival (s)')
    ylabel('pressure (pa)')
    axb51 = addbox(ax5, [0 0 0.25 0.2]);
    text(0.15, 0.45, sprintf('%.2f--%.2f Hz', fc(1), fc(2)), 'FontSize', 12);
    axb52 = addbox(ax5, [0.75 0 0.25 0.2]);
    text(0.2, 0.45, sprintf('SNR = %.2f', A(I(J),J)), 'FontSize', 12);
    axes(axb51)
    set(ax5, 'FontSize', 12, 'TickDir', 'out')
    
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_%s', mfilename, savename), [], [], 2, [], ...
        'epstopdf');
end
end

% compute signal-to-noise ratio
%
% INPUT:
% t             time
% x             signal
% win_select    time window for time search
% t_begin       begining time for noise window
% t_end         ending time for signal window
%
% OUTPUT:
% s             signal-to-noise ratio
% t_max         time that give the maximum signal-to-noise ratio
function [s, t_max] = snrvar(t, x, win_select, t_begin, t_end)
tt = t(and(t >= win_select(1), t <= win_select(2)));
ss = zeros(size(tt));
for ii = 1:length(tt)
    ss(ii) = var(x(and(t >= tt(ii), t < t_end))) / ...
        var(x(and(t >= t_begin, t < tt(ii))));
end
s = max(ss);
t_max = tt(ss == s);
end