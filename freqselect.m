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
% Last modified by sirawich-at-princeton.edu, 05/31/2022

% Nyquist frequency
fNq = fs/2;

% list of corner frequency candidates
fcs = 0:0.05:2.05;

% bandpass/lowpass SNR
A = NaN(length(fcs), length(fcs));
T = NaN(length(fcs), length(fcs));

% bandstop/highpass SNR
B = NaN(length(fcs), length(fcs));
U = NaN(length(fcs), length(fcs));

for ii = 1:length(fcs)
    for jj = (ii+1):length(fcs)
        % for zero lower corner frequency: lowpass
        if fcs(ii) == 0 && jj < length(fcs)
            xf = lowpass(x, fs, fcs(jj), 2, 2, 'butter', 'linear');
        % for the highest upper corner frequency: high pass
        elseif fcs(ii) > 0 && jj == length(fcs)
            xf = hipass(x, fs, fcs(ii), 2, 2, 'butter', 'linear');
        % bandpass
        elseif fcs(ii) > 0 && fcs(jj) / fcs(ii) >= 1.5
            xf = bandpass(x, fs, fcs(ii), fcs(jj), 2, 2, 'butter', 'linear');
        % skip if the window is too narrow or [0 Inf]
        else
            continue
        end
        xs = x - xf;
        if jj < length(fcs)
            fmid = max(fcs(ii), 0.05); %(2 * fcs(ii) + fcs(jj)) / 3;
        else
            fmid = max(fcs(ii), 0.05); %(2 * fcs(ii) + fNq) / 3;
        end
        halfwin = 2/ fmid;
        [A(jj, ii), T(jj, ii)] = snrvar(t, xf, [-1 1] * halfwin/2, ...
            -60, 80, 1 * halfwin);
        [B(jj, ii), U(jj, ii)] = snrvar(t, xs, [-1 1] * halfwin/2, ...
            -60, 80, 1 * halfwin);
    end
end

% pick the best corner frequencies
[M, I] = max(A./B, [], 1);
[MM, J] = max(M);
fc = [fcs(J) fcs(I(J))];

% best SNR
s = A(I(J), J);

%% visualize the result
if plt
    % locate the best pixel on the plot
    [xx, yy] = boxcorner(fc(1) + 0.025 * [-1 1], fc(2) + 0.025 * [-1 1]);
    
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [16 12 8 12])
    
    % title
    ax0 = subplot('Position', [0.08 0.97 0.88 0.01]);
    title(titlename)
    set(ax0, 'FontSize', 12, 'Color', 'none')
    ax0.XAxis.Visible = 'off';
    ax0.YAxis.Visible = 'off';
    
    % bandpass SNR
    ax1 = subplot('Position', [0.08 0.77 0.37 0.19]);
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
    
    ax2 = subplot('Position', [0.58 0.77 0.37 0.19]);
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
    
    ax3 = subplot('Position', [0.08 0.52 0.37 0.19]);
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
    
    ax4 = subplot('Position', [0.58 0.52 0.37 0.19]);
    im4 = imagesc(fcs, fcs, T);
    axis xy
    xlabel('lower corner frequency (Hz)')
    ylabel('upper corner frequency (Hz)')
    title('optimal time for bandpass')
    c4 = colorbar(ax4, 'EastOutside');
    colormap(ax4, kelicol)
    clim = [-1 1] * max(abs(T), [], 'all');
    ax4.CLim = clim;
    setimagenan(ax4, im4, [1 1 1], clim(1), clim(2));
    [~,v4] = vline(ax4, ax4.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    [~,h4] = hline(ax4, ax4.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4]);
    plot(ax4, xx, yy, '-r', 'LineWidth', 1)
    c4.Label.String = 'time (s)';
    set(ax4, 'FontSize', 12, 'TickDir', 'out')
    set(c4, 'FontSize', 12, 'TickDir', 'out')
    
    % power spectral density
    ax5 = subplot('Position', [0.08 0.35 0.88 0.09]);
    % fix the upper corner frequency for hipass case
    if fc(2) == fcs(end)
        fc(2) = fNq;
    end
    p = specdensplot(x, round(fs * 40), fs, round(fs * 40), 70, 10, 's');
    p(1).Color = [0 0 0.5];
    p(1).LineWidth = 1;
    delete(p(2))
    delete(p(3))
    delete(p(4))
    hold on
    grid on
    % plot cosine filter when removing instrument response
    ff = (0:0.001:(fs/2))';
    fac = zeros(size(ff));
    fl = [0.01 0.02 10 20];
    for ii = 1:length(ff)
        freq = ff(ii);
        if freq < fl(1)
            fac(ii) = 0.0;
        elseif freq >= fl(1) && freq <= fl(2)
            fac(ii) = 0.5 * (1 - cos(pi*(freq - fl(1)) / (fl(2) - fl(1))));
        elseif freq >= fl(3) && freq <= fl(4)
            fac(ii) = 0.5 * (1 + cos(pi*(freq - fl(3)) / (fl(4) - fl(3))));
        elseif freq > fl(4)
            fac(ii) = 0.0;
        else
            fac(ii) = 1.0;
        end
    end
    fac = ax5.YLim(1) + fac * 0.5 * (ax5.YLim(2) - ax5.YLim(1));
    semilogx(ax5, ff, fac, 'LineWidth', 1, 'Color', 'k');
    xlim([fs/round(fs * 40) fs/2])
    ylabel('spectral density (Pa^2/Hz)')
    [~, v5] = vline(ax5, fc, 'LineWidth', 1, 'Color', [0 0.6 1]);
    set(ax5, 'FontSize', 12, 'TickDir', 'out', 'Color', 'none')
    % add x-tick label at the axis limits
    if ax5.XTick(1) / ax5.XLim(1) >= 2
        ax5.XTick = [ax5.XLim(1) ax5.XTick];
    end
    if ax5.XLim(2) / ax5.XTick(end) >= 2
        ax5.XTick = [ax5.XTick ax5.XLim(2)];
    end
    % add period axis
    ax5s = doubleaxes(ax5);
    inverseaxis(ax5s.XAxis, 'period (s)');
    % highlight the window
    ax5h = doubleaxes(ax5);
    [xbox, ybox] = boxcorner(fc, ax5.YLim);
    pgon = polyshape(xbox, ybox);
    bx = plot(ax5h, pgon, 'FaceColor', [1 0.9 0.4], 'FaceAlpha', 0.4, ...
        'EdgeAlpha', 0);
    ax5h.XAxis.Visible = 'off';
    ax5h.YAxis.Visible = 'off';
    ax5h.XScale = 'log';
    set(ax5h, 'Box', 'on', 'TickDir', 'both', 'XLim', ax5.XLim, 'YLim', ...
        ax5.YLim, 'Position', ax5.Position);
    
    % filtered seismogram
    ax6 = subplot('Position', [0.08 0.2 0.88 0.09]);
    if fc(1) > 0 && fc(2) < fNq
        xf = bandpass(x, fs, fc(1), fc(2), 2, 2, 'butter', 'linear');
    elseif fc(1) == 0 && fc(2) < fNq
        xf = lowpass(x, fs, fc(2), 2, 2, 'butter', 'linear');
    elseif fc(1) > 0 && fc(2) == fNq
        xf = hipass(x, fs, fc(1), 2, 2, 'butter', 'linear');
    else
        keyboard;
    end
    plot(ax6, t, xf, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
    hold on
    grid on
    
    % highlight noise window
    t0 = T(I(J), J);
    fmid = max(fc(1), 0.05);%(2*fc(1) + fc(2)) / 3;
    halfwin = 2 / fmid;
    t_start = max(-60, t0 - 1 * halfwin);
    t_end = min(80, t0 + 1 * halfwin);
    xn = xf(and(t >= t_start, t < t0));
    tn = t(and(t >= t_start, t < t0));
    plot(ax6, tn, xn, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
    
    % highlight signal window
    xs = xf(and(t >= t0, t < t_end));
    ts = t(and(t >= t0, t < t_end));
    plot(ax6, ts, xs, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
    
    % fix y-limits
    ylimits = ax6.YLim;
    ylimits = [-1 1] * abs(max(ylimits));
    ylim(ylimits)
    
    % fix x_limits
    xlimits = [max(10 * round(-20 / 10 * 1/fc(1)), min(t)), ...
        min(10 * round(20 / 10 * 1/fc(1)), max(t))];
    xlim(xlimits)
    
    [~, v6] = vline(ax6, T(I(J), J), 'LineWidth', 1, 'Color', [0 0.6 1]);
    xlabel('time since picked first arrival (s)')
    ylabel('pressure (pa)')
    axb61 = addbox(ax6, [0 0 0.25 0.2]);
    text(0.15, 0.45, sprintf('%.2f--%.2f Hz', fc(1), fc(2)), 'FontSize', 12);
    axb62 = addbox(ax6, [0.75 0 0.25 0.2]);
    if A(I(J), J) > 5
        text(0.2, 0.45, sprintf('SNR = %d', round(A(I(J),J))), 'FontSize', 12);
    else
        text(0.2, 0.45, sprintf('SNR = %.2f', A(I(J),J)), 'FontSize', 12);
    end
    axes(axb61)
    set(ax6, 'FontSize', 12, 'TickDir', 'out')
    
    % unfiltered seismogram
    ax7 = subplot('Position', [0.08 0.06 0.88 0.09]);
    plot(ax7, t, x, 'Color', 'k', 'LineWidth', 1)
    grid on
    set(ax7, 'FontSize', 12, 'TickDir', 'out')
    xlim(ax6.XLim)
    ylim([-1 1] * abs(max(x)))
    [~, v7] = vline(ax7, T(I(J), J), 'LineWidth', 1, 'Color', [0 0.6 1]);
    xlabel('time since picked first arrival (s)')
    ylabel('pressure (pa)')
    axb71 = addbox(ax7, [0 0 0.25 0.2]);
    text(0.45, 0.45, 'all', 'FontSize', 12);
    axes(axb61)
    
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_%s', mfilename, savename), [], [], 2, [], ...
        'epstopdf');
else
    if fc(2) == fcs(end)
        fc(2) = fNq;
    end
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
% t_length      length of noise and signal windows
%
% OUTPUT:
% s             signal-to-noise ratio
% t_max         time that give the maximum signal-to-noise ratio
function [s, t_max] = snrvar(t, x, win_select, t_begin, t_end, t_length)
tt = t(and(t >= win_select(1), t <= win_select(2)));
ss = zeros(size(tt));
for ii = 1:length(tt)
    ss(ii) = var(x(and(t >= tt(ii), t < min(t_end, tt(ii) + t_length)))) / ...
        var(x(and(t >= max(t_begin, tt(ii) - t_length), t < tt(ii))));
end
s = max(ss);
t_max = tt(ss == s);
end