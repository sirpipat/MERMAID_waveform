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
% Last modified by sirawich-at-princeton.edu, 10/24/2022

% Nyquist frequency
fNq = fs/2;

% list of corner frequency candidates
fcs = 0.4:0.05:2.05;

% pass SNR and optimal time
A = NaN(length(fcs), length(fcs));
T = NaN(length(fcs), length(fcs));

% stop SNR and optimal time
B = NaN(length(fcs), length(fcs));
U = NaN(length(fcs), length(fcs));

% detrend the original signal
x = detrend(x .* shanning(length(x), 0.05, 0), 1);

% remove frequency content below the lowest lower corner frequency
x = hipass(x, fs, fcs(1), 2, 2, 'butter', 'linear');

for ii = 1:length(fcs)
    for jj = (ii+1):length(fcs)
        % for zero lower corner frequency: lowpass
        if fcs(ii) == 0 && jj < length(fcs)
            xf = lowpass(x, fs, fcs(jj), 2, 2, 'butter', 'linear');
        % for the highest upper corner frequency: high pass
        elseif fcs(ii) > 0 && jj == length(fcs)
            continue
            xf = hipass(x, fs, fcs(ii), 2, 2, 'butter', 'linear');
        % bandpass
        elseif fcs(ii) > 0 && fcs(jj) - fcs(ii) >= 0.4995
            xf = bandpass(x, fs, fcs(ii), fcs(jj), 2, 2, 'butter', 'linear');
        % skip if the window is too narrow or [0 Inf]
        else
            continue
        end
        xs = bandstop(x, fs, fcs(ii), fcs(jj), 2, 2, 'butter', 'linear');
        if jj < length(fcs)
            fmid = max(fcs(ii), 0.05); %(2 * fcs(ii) + fcs(jj)) / 3;
        else
            fmid = max(fcs(ii), 0.05); %(2 * fcs(ii) + fNq) / 3;
        end
        halfwin = 2/ fmid;
        [A(jj, ii), T(jj, ii)] = snrvar(t, xf, [-1 1] * halfwin/2, ...
            -60, 80, 1 * halfwin);
        % for bandstop window length = 2 / lowest freq exist in bandstop
        if jj < length(fcs)
            [B(jj, ii), U(jj, ii)] = snrvar(t, xs, [-1 1] * 5/2, ...
                -60, 80, 1 * 5);
        else
            [B(jj, ii), U(jj, ii)] = snrvar(t, xs, [-1 1] / fcs(ii), ...
                -60, 80, 1 * 2 / fcs(ii));
        end
    end
end

% pick the best corner frequencies
[M, I] = max(A./1, [], 1);
[MM, J] = max(M);
fc = [fcs(J) fcs(I(J))];

% best SNR
s = A(I(J), J);

%% visualize the result
if plt
    % locate the best pixel on the plot
    [xx, yy] = boxcorner(fc(1) + 0.025 * [-1 1], fc(2) + 0.025 * [-1 1]);
    
%     figure(1)
%     clf
%     set(gcf, 'Units', 'inches', 'Position', [16 12 8 12])
%     
%     % title
%     ax0 = subplot('Position', [0.08 0.97 0.88 0.01]);
%     title(titlename)
%     set(ax0, 'FontSize', 12, 'Color', 'none')
%     ax0.XAxis.Visible = 'off';
%     ax0.YAxis.Visible = 'off';
%     
%     % bandpass SNR
%     ax1 = subplot('Position', [0.08 0.77 0.37 0.19]);
%     im1 = imagesc(fcs, fcs, log10(A));
%     axis xy
%     xlabel('lower corner frequency (Hz)')
%     ylabel('upper corner frequency (Hz)')
%     title('pass')
%     setimagenan(ax1, im1, [1 1 1]);
%     c1 = colorbar(ax1, 'EastOutside');
%     [~,v1] = vline(ax1, ax1.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     [~,h1] = hline(ax1, ax1.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     plot(ax1, xx, yy, '-r', 'LineWidth', 1)
%     c1.Label.String = 'log_{10} SNR';
%     set(ax1, 'FontSize', 12, 'TickDir', 'out')
%     set(c1, 'FontSize', 12, 'TickDir', 'out')
%     
%     ax2 = subplot('Position', [0.58 0.77 0.37 0.19]);
%     im2 = imagesc(fcs, fcs, log10(B));
%     axis xy
%     xlabel('lower corner frequency (Hz)')
%     ylabel('upper corner frequency (Hz)')
%     title('stop')
%     setimagenan(ax2, im2, [1 1 1]);
%     c2 = colorbar(ax2, 'EastOutside');
%     [~,v2] = vline(ax2, ax2.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     [~,h2] = hline(ax2, ax2.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     plot(ax2, xx, yy, '-r', 'LineWidth', 1)
%     c2.Label.String = 'log_{10} SNR';
%     set(ax2, 'FontSize', 12, 'TickDir', 'out')
%     set(c2, 'FontSize', 12, 'TickDir', 'out')
%     
%     ax3 = subplot('Position', [0.08 0.52 0.37 0.19]);
%     im3 = imagesc(fcs, fcs, log10(A./B));
%     axis xy
%     xlabel('lower corner frequency (Hz)')
%     ylabel('upper corner frequency (Hz)')
%     title('ratio')
%     setimagenan(ax3, im3, [1 1 1]);
%     c3 = colorbar(ax3, 'EastOutside');
%     [~,v3] = vline(ax3, ax3.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     [~,h3] = hline(ax3, ax3.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     plot(ax3, xx, yy, '-r', 'LineWidth', 1)
%     c3.Label.String = 'log_{10} SNR';
%     set(ax3, 'FontSize', 12, 'TickDir', 'out')
%     set(c3, 'FontSize', 12, 'TickDir', 'out')
%     
%     ax4 = subplot('Position', [0.58 0.52 0.37 0.19]);
%     im4 = imagesc(fcs, fcs, T);
%     axis xy
%     xlabel('lower corner frequency (Hz)')
%     ylabel('upper corner frequency (Hz)')
%     title('optimal time for pass')
%     c4 = colorbar(ax4, 'EastOutside');
%     colormap(ax4, kelicol)
%     clim = [-1 1] * max(max(abs(T)));
%     ax4.CLim = clim;
%     setimagenan(ax4, im4, [1 1 1], clim(1), clim(2));
%     [~,v4] = vline(ax4, ax4.XTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     [~,h4] = hline(ax4, ax4.YTick, 'LineWidth', 1, 'LineStyle', ':', ...
%         'Color', [0.4 0.4 0.4]);
%     plot(ax4, xx, yy, '-r', 'LineWidth', 1)
%     c4.Label.String = 'time (s)';
%     set(ax4, 'FontSize', 12, 'TickDir', 'out')
%     set(c4, 'FontSize', 12, 'TickDir', 'out')
%     
%     % power spectral density
%     ax5 = subplot('Position', [0.08 0.35 0.88 0.09]);
%     if true
%         % fix the upper corner frequency for hipass case
%         if fc(2) == fcs(end)
%             fc(2) = fNq;
%         end
%         p = specdensplot(x, round(fs * 40), fs, round(fs * 40), 70, 10, 's');
%         p(1).Color = [0 0 0.5];
%         p(1).LineWidth = 1;
%         delete(p(2))
%         delete(p(3))
%         delete(p(4))
%         hold on
%         grid on
%         % plot cosine filter when removing instrument response
%         ff = (0:0.001:(fs/2))';
%         fac = zeros(size(ff));
%         fl = [0.01 0.02 10 20];
%         for ii = 1:length(ff)
%             freq = ff(ii);
%             if freq < fl(1)
%                 fac(ii) = 0.0;
%             elseif freq >= fl(1) && freq <= fl(2)
%                 fac(ii) = 0.5 * (1 - cos(pi*(freq - fl(1)) / (fl(2) - fl(1))));
%             elseif freq >= fl(3) && freq <= fl(4)
%                 fac(ii) = 0.5 * (1 + cos(pi*(freq - fl(3)) / (fl(4) - fl(3))));
%             elseif freq > fl(4)
%                 fac(ii) = 0.0;
%             else
%                 fac(ii) = 1.0;
%             end
%         end
%         fac = ax5.YLim(1) + fac * 0.5 * (ax5.YLim(2) - ax5.YLim(1));
%         semilogx(ax5, ff, fac, 'LineWidth', 1, 'Color', 'k');
%         xlim([fs/round(fs * 40) fs/2])
%         ylabel('spectral density (Pa^2/Hz)')
%         [~, v5] = vline(ax5, fc, 'LineWidth', 1, 'Color', [0 0.6 1]);
%         set(ax5, 'FontSize', 12, 'TickDir', 'out', 'Color', 'none')
%         % add x-tick label at the axis limits
%         if ax5.XTick(1) / ax5.XLim(1) >= 2
%             ax5.XTick = [ax5.XLim(1) ax5.XTick];
%         end
%         if ax5.XLim(2) / ax5.XTick(end) >= 2
%             ax5.XTick = [ax5.XTick ax5.XLim(2)];
%         end
%         % add period axis
%         ax5s = doubleaxes(ax5);
%         inverseaxis(ax5s.XAxis, 'period (s)');
%         % highlight the window
%         ax5h = doubleaxes(ax5);
%         [xbox, ybox] = boxcorner(fc, ax5.YLim);
%         pgon = polyshape(xbox, ybox);
%         bx = plot(ax5h, pgon, 'FaceColor', [1 0.9 0.4], 'FaceAlpha', 0.4, ...
%             'EdgeAlpha', 0);
%         ax5h.XAxis.Visible = 'off';
%         ax5h.YAxis.Visible = 'off';
%         ax5h.XScale = 'log';
%         set(ax5h, 'Box', 'on', 'TickDir', 'both', 'XLim', ax5.XLim, 'YLim', ...
%             ax5.YLim, 'Position', ax5.Position);
%     else
%         plot(ax5, t, x, 'Color', 'k', 'LineWidth', 1)
%         grid on
%         set(ax5, 'FontSize', 12, 'TickDir', 'out')
%         %xlim(ax6.XLim)
%         ylim([-1 1] * abs(max(x)))
%         [~, v5] = vline(ax5, T(I(J), J), 'LineWidth', 1, 'Color', [0 0.6 1]);
%         xlabel('time since picked first arrival (s)')
%         ylabel('pressure (pa)')
%         title('original: highpass 0.4 Hz')
%         axb51 = addbox(ax5, [0 0 0.25 0.2]);
%         text(0.15, 0.45, sprintf('0.4--%.2f Hz', fNq), 'FontSize', 12);
%         axes(axb51)
%     end
%     
%     % filtered seismogram
%     ax6 = subplot('Position', [0.08 0.2 0.88 0.09]);
%     if fc(1) > 0 && fc(2) < fNq
%         xf = bandpass(x, fs, fc(1), fc(2), 2, 2, 'butter', 'linear');
%     elseif fc(1) == 0 && fc(2) < fNq
%         xf = lowpass(x, fs, fc(2), 2, 2, 'butter', 'linear');
%     elseif fc(1) > 0 && fc(2) == fNq
%         xf = hipass(x, fs, fc(1), 2, 2, 'butter', 'linear');
%     else
%         keyboard;
%     end
%     plot(ax6, t, xf, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
%     hold on
%     grid on
%     
%     % highlight noise window
%     t0 = T(I(J), J);
%     fmid = max(fc(1), 0.05);%(2*fc(1) + fc(2)) / 3;
%     halfwin = 2 / fmid;
%     t_start = max(-60, t0 - 1 * halfwin);
%     t_end = min(80, t0 + 1 * halfwin);
%     xn = xf(and(t >= t_start, t < t0));
%     tn = t(and(t >= t_start, t < t0));
%     plot(ax6, tn, xn, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
%     
%     % highlight signal window
%     xs = xf(and(t >= t0, t < t_end));
%     ts = t(and(t >= t0, t < t_end));
%     plot(ax6, ts, xs, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
%     
%     % fix y-limits
%     ylimits = ax6.YLim;
%     ylimits = [-1 1] * abs(max(ylimits));
%     ylim(ylimits)
%     
%     % fix x_limits
%     xlimits = [max(10 * round(-20 / 10 * 1/fc(1)), min(t)), ...
%         min(10 * round(20 / 10 * 1/fc(1)), max(t))];
%     xlim(xlimits)
%     
%     [~, v6] = vline(ax6, T(I(J), J), 'LineWidth', 1, 'Color', [0 0.6 1]);
%     xlabel('time since picked first arrival (s)')
%     ylabel('pressure (pa)')
%     axb61 = addbox(ax6, [0 0 0.25 0.2]);
%     text(0.15, 0.45, sprintf('%.2f--%.2f Hz', fc(1), fc(2)), 'FontSize', 12);
%     axb62 = addbox(ax6, [0.75 0 0.25 0.2]);
%     if A(I(J), J) > 5
%         text(0.2, 0.45, sprintf('SNR = %d', round(A(I(J),J))), 'FontSize', 12);
%     else
%         text(0.2, 0.45, sprintf('SNR = %.2f', A(I(J),J)), 'FontSize', 12);
%     end
%     axes(axb61)
%     set(ax6, 'FontSize', 12, 'TickDir', 'out')
%     
%     % bandstop seismogram
%     xt = x - xf;
%     ax7 = subplot('Position', [0.08 0.06 0.88 0.09]);
%     plot(ax7, t, xt, 'Color', [0.4 0.4 0.4], 'LineWidth', 1)
%     grid on
%     hold on
%     set(ax7, 'FontSize', 12, 'TickDir', 'out')
%     
%     % highlight noise window
%     t0 = U(I(J), J);
%     fmid = fcs(1);
%     halfwin = 2 / fmid;
%     t_start = max(-60, t0 - 1 * halfwin);
%     t_end = min(80, t0 + 1 * halfwin);
%     xn = xt(and(t >= t_start, t < t0));
%     tn = t(and(t >= t_start, t < t0));
%     plot(ax7, tn, xn, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
%     
%     % highlight signal window
%     xs = xt(and(t >= t0, t < t_end));
%     ts = t(and(t >= t0, t < t_end));
%     plot(ax7, ts, xs, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
%     
%     ax7.XLim = ax6.XLim;
%     %xlim([-20 20])
%     ylim([-1 1] * abs(max(x-xf)))
%     [~, v7] = vline(ax7, U(I(J), J), 'LineWidth', 1, 'Color', [0 0.6 1]);
%     xlabel('time since picked first arrival (s)')
%     ylabel('pressure (pa)')
%     axb71 = addbox(ax7, [0 0 0.25 0.2]);
%     text(0.40, 0.45, 'stop', 'FontSize', 12);
%     axb72 = addbox(ax7, [0.75 0 0.25 0.2]);
%     if B(I(J), J) > 5
%         text(0.2, 0.45, sprintf('SNR = %d', round(B(I(J),J))), 'FontSize', 12);
%     else
%         text(0.2, 0.45, sprintf('SNR = %.2f', B(I(J),J)), 'FontSize', 12);
%     end
%     axes(axb71)
%     
%     ax5.XLim = ax6.XLim;
%     
%     set(gcf, 'Renderer', 'painters')
%     figdisp(sprintf('%s_%s', mfilename, savename), [], [], 2, [], ...
%         'epstopdf');
    
    fig = figure(2);
    clf;
    set(fig, 'Units', 'inches', 'Position', [0 1 10 6]);
    sp0 = subplot('Position', [0.06 0.94 0.9 0.02]);
    title(titlename)
    set(sp0, 'FontSize', 12, 'Color', 'none')
    sp0.XAxis.Visible = 'off';
    sp0.YAxis.Visible = 'off';
    
    sp1 = subplot('Position', [0.06 0.5 0.5 0.36]);
    [~,~,~,~,F] = timspecplot_ns(x, 400, fs, 400,0.7, t(1), 's', 'log');
    hold on
    % grid
    [~,vvv] = vline(sp1, sp1.XTick, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.5, 'LineStyle', ':');
    [~,hhh] = hline(sp1, sp1.YTick, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.5, 'LineStyle', ':');
    % mark where the corner frequencies are
    hline(sp1, lin2logpos(fc, F(2), F(end)), 'Color', 'k', 'LineWidth', 1);
    % fix the precision of the time on XAxis label
    sp1.XAxis.Label.String = sprintf('time since picked arrival (s): %d s window', round(400/fs));
    sp1.YAxis.Label.String = 'frequency (Hz)';
    sp1.Title.String = '';
    set(sp1, 'FontSize', 12, 'TickDir', 'out', 'XAxisLocation', 'top')
    % insert colorbar
    cc1 = colorbar;
    cc1.Label.String = 'spectral density 10log_{10}(Pa^2/Hz)';
    cPosition = cc1.Position;
    
	sp2 = subplot('Position', [0.67 0.5 0.25 0.36]);
    p = specdensplot(x, 400, fs, 400, 70, 10, 's');
    grid on
    set(p(1), 'LineWidth', 1, 'Color', 'r')
    set(p(2), 'LineWidth', 0.75, 'Color', [0.4 0.4 0.4])
    set(p(3), 'LineWidth', 0.75, 'Color', [0.4 0.4 0.4])
    [~, vv2] = vline(sp2, fc, 'LineWidth', 1, 'LineStyle', '-', 'Color', 'k');
    % fix the precision of the time on XAxis label
    sp2.XAxis.Label.String = sprintf('frequency (Hz): %d s window', round(400/fs));
    sp2.YAxis.Label.String = 'spectral energy 10log_{10}(Pa^2/Hz)';
    set(sp2, 'FontSize', 12, 'TickDir', 'out', 'XAxisLocation', 'top', ...
        'YAxisLocation', 'right')
    
	sp3 = subplot('Position', [0.06 0.09 0.4389 0.36]);
    plot(sp3, t, x / max(abs(x(and(t >= -20, t < 20)))), ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1)
    hold on
    % highlight noise window
    [s_raw,t_max_raw] = snrvar(t, x, [-1 1] * 5/2, -60, 80, 1 * 5);
    t_start = max(-60, t_max_raw - 1 * 5);
    t_end = min(80, t_max_raw + 1 * 5);
    scale_all = max(abs(x(and(t >= -20, t < 20))));
    xn = x(and(t >= t_start, t < t_max_raw)) / scale_all;
    tn = t(and(t >= t_start, t < t_max_raw));
    plot(sp3, tn, xn, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
    % highlight signal window
    xs = x(and(t >= t_max_raw, t < t_end)) / scale_all;
    ts = t(and(t >= t_max_raw, t < t_end));
    plot(sp3, ts, xs, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
    plot(sp3, [1 1] * t_max_raw, [-0.8 0.8], 'Color', [0 0.6 1], ...
        'LineWidth', 1)
    text(sp3, -19, 0.7, 'all', 'FontSize', 11);
    
    if fc(1) > 0 && fc(2) < fNq
        xf = bandpass(x, fs, fc(1), fc(2), 2, 2, 'butter', 'linear');
    elseif fc(1) == 0 && fc(2) < fNq
        xf = lowpass(x, fs, fc(2), 2, 2, 'butter', 'linear');
    elseif fc(1) > 0 && fc(2) == fNq
        xf = hipass(x, fs, fc(1), 2, 2, 'butter', 'linear');
    else
        keyboard;
    end
    scale_pass = max(abs(xf(and(t >= -20, t < 20))));
    
    % bandstop signal
    xt = bandstop(x, fs, fc(1), fc(2), 2, 2, 'butter', 'linear');
    scale_stop = max(abs(xt(and(t >= -20, t < 20))));
    
    % determine which scale to use
    scale_used = max(scale_pass, scale_stop);
    
    plot(sp3, t, xf / scale_used - 2, ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1)
    % highlight noise window
    t0 = T(I(J), J);
    fmid = max(fc(1), 0.05);
    halfwin = 2 / fmid;
    t_start = max(-60, t0 - 1 * halfwin);
    t_end = min(80, t0 + 1 * halfwin);
    xn = xf(and(t >= t_start, t < t0)) / scale_used;
    tn = t(and(t >= t_start, t < t0));
    plot(sp3, tn, xn - 2, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
    % highlight signal window
    xs = xf(and(t >= t0, t < t_end)) / scale_used;
    ts = t(and(t >= t0, t < t_end));
    plot(sp3, ts, xs - 2, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
    plot(sp3, [1 1] * t0, [-0.8 0.8] - 2, 'Color', [0 0.6 1], ...
        'LineWidth', 1)
    text(sp3, -19, -1.3, sprintf('%.2f-%.2f Hz (x %.2f)', fc(1), fc(2), ...
        scale_all / scale_used), 'FontSize', 11);
    
    
    
    % bandstop signal
    plot(sp3, t, xt / scale_used - 4, ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1)
    % highlight noise window
    t0 = U(I(J), J);
    fmid = fcs(1);
    halfwin = 2 / fmid;
    t_start = max(-60, t0 - 1 * halfwin);
    t_end = min(80, t0 + 1 * halfwin);
    
    xn = xt(and(t >= t_start, t < t0)) / scale_used;
    tn = t(and(t >= t_start, t < t0));
    plot(sp3, tn, xn - 4, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
    % highlight signal window
    xs = xt(and(t >= t0, t < t_end)) / scale_used;
    ts = t(and(t >= t0, t < t_end));
    plot(sp3, ts, xs - 4, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
    plot(sp3, [1 1] * t0, [-0.8 0.8] - 4, 'Color', [0 0.6 1], ...
        'LineWidth', 1)
    text(sp3, -19, -3.3, sprintf('bandstop (x %.2f)', scale_all / scale_used), ...
        'FontSize', 11);
    
    grid on
    sp3.YTick = [-4 -2 0];
    sp3.YTickLabel = [round(B(I(J),J)) round(s) round(s_raw)];
    sp3.XLabel.String = 'time since picked arrival (s)';
    sp3.YLabel.String = 'SNR';
    sp3.YAxis.TickLabelRotation = 90;
    set(sp3, 'FontSize', 12, 'TickDir', 'out', 'XLim', [-20 20], ...
        'YLim', [-5.2 1.2])
    
	sp4 = subplot('Position', [0.67 0.09 0.3 0.36]);
    iim = imagesc(fcs(1:(end-10)), fcs(11:end), log10(A(11:end, 1:(end-10))));
    axis xy
    xlabel('lower corner frequency (Hz)')
    ylabel('upper corner frequency (Hz)')
    colormap(sp4, 'gray')
    setimagenan(sp4, iim, [1 1 1]);
    cc4 = colorbar(sp4, 'EastOutSide');
    cc4.Label.String = 'log_{10} SNR';
    cc4.Label.FontSize = 12;
    hold on
    plot(sp4, [0 2], [0 2], 'Color', 'k', 'LineWidth', 2)
    add_fc_grid(sp4)
    plot(sp4, xx, yy, '-r', 'LineWidth', 2)
    sp4.XLim = [0.375 1.525];
    sp4.YLim = [0.875 2.025];
    set(sp4, 'FontSize', 12, 'TickDir', 'out')
    set(gcf, 'Renderer', 'painters')
    
    % align the axes
    sp1.Position = [0.06 0.50 0.48 0.36];
    sp2.Position = [0.66 0.50 0.25 0.36];
    sp3.Position = [0.06 0.09 0.48 0.36];
    sp4.Position = [0.66 0.09 0.25 0.36];
    
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

% add grid to the frequency grid plot to an axes ax
%
% INPUT:
% ax        target axes
function add_fc_grid(ax)
min_lower_fc = 0.375;
max_lower_fc = 1.525;
min_upper_fc = 0.875;
max_upper_fc = 2.025;
df = 0.05;

x = [];
y = [];

% vertical grid
for freq = min_lower_fc:df:max_lower_fc
    x = [x freq freq NaN];
    y = [y max_upper_fc max(min_upper_fc, freq) NaN];
end

% horizontal grid
for freq = min_upper_fc:df:max_upper_fc
    x = [x min_lower_fc min(max_lower_fc, freq) NaN];
    y = [y freq freq NaN];
end

plot(ax, x, y, 'LineWidth', 0.5, 'Color', 'k')
end