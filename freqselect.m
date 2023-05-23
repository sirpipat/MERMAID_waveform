function [fc, s] = freqselect(t, x, fs, plt, titlename, savename)
% [fc, s] = FREQSELECT(t, x, fs, plt, titlename, savename)
%
% Figures out the frequency band where the signal stands out the most from
% the background noise.
% 
% INPUT:
% t             time [for an example, go between -100 and 100]
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
% Last modified by sirawich-at-princeton.edu, 05/18/2023

% Check the sampling rate consistency
diferm(unique(diff(t)),1/fs)
    
% Nyquist frequency
fNq = fs/2;

% list of corner frequency candidates
delf=0.05;
fcs = 0.4:delf:2.05;

% spread i.e. minimum bandwidth
fspread = 0.4995;

% number of poles and passes
npoles = 4;
npasss = 2;

% pass SNR and optimal time
A = NaN(length(fcs), length(fcs));
T = NaN(length(fcs), length(fcs));

% stop SNR and optimal time
B = NaN(length(fcs), length(fcs));
U = NaN(length(fcs), length(fcs));

% detrend the original signal
x = detrend(x .* shanning(length(x), 0.05, 0), 1);

% remove frequency content below the lowest lower corner frequency
x = hipass(x, fs, fcs(1), npoles, npasss, 'butter', 'linear');

% try all lower corner frequency candidates
for ii = 1:length(fcs)
    % try all upper corner frequency candidates
    for jj = (ii+1):length(fcs)
        % for zero lower corner frequency: lowpass
        if fcs(ii) == 0 && jj < length(fcs)
            xf = lowpass(x, fs, fcs(jj), npoles, npasss, 'butter', 'linear');
        % for the highest upper corner frequency: high pass
        elseif fcs(ii) > 0 && jj == length(fcs)
            continue
            xf = hipass(x, fs, fcs(ii), npoles, npasss, 'butter', 'linear');
        % bandpass
        elseif fcs(ii) > 0 && fcs(jj) - fcs(ii) >= fspread
            xf = bandpass(x, fs, fcs(ii), fcs(jj), npoles, npasss, 'butter', 'linear');
        % skip if the window is too narrow or [0 Inf]
        else
            continue
        end
        xs = bandstop(x, fs, fcs(ii), fcs(jj), npoles, npasss, 'butter', 'linear');
        if jj < length(fcs)
            fmid = max(fcs(ii), delf); %(2 * fcs(ii) + fcs(jj)) / 3;
        else
            fmid = max(fcs(ii), delf); %(2 * fcs(ii) + fNq) / 3;
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
    [~,vvv] = vline(sp1, sp1.XTick, 'Color', [0.25 0.25 0.25], ...
        'LineWidth', 1.5, 'LineStyle', ':');
    [~,hhh] = hline(sp1, sp1.YTick, 'Color', [0.25 0.25 0.25], ...
        'LineWidth', 1.5, 'LineStyle', ':');
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
        xf = bandpass(x, fs, fc(1), fc(2), npoles, npasss, 'butter', 'linear');
    elseif fc(1) == 0 && fc(2) < fNq
        xf = lowpass(x, fs, fc(2), npoles, npasss, 'butter', 'linear');
    elseif fc(1) > 0 && fc(2) == fNq
        xf = hipass(x, fs, fc(1), npoles, npasss, 'butter', 'linear');
    else
        keyboard;
    end
    scale_pass = max(abs(xf(and(t >= -20, t < 20))));
    
    % bandstop signal
    xt = bandstop(x, fs, fc(1), fc(2), npoles, npasss, 'butter', 'linear');
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
    
    % added plot labels
    sp1b = boxedlabel(sp1, 'northwest', 0.27, [], 'a', 'FontSize', 12);
    axes(sp1b)
    sp2b = boxedlabel(sp2, 'northeast', 0.27, [], 'b', 'FontSize', 12);
    axes(sp2b)
    sp3b = boxedlabel(sp3, 'southwest', 0.27, [], 'c', 'FontSize', 12);
    axes(sp3b)
    sp4b = boxedlabel(sp4, 'southeast', 0.27, [], 'd', 'FontSize', 12);
    axes(sp4b)
    
    figdisp(sprintf('%s_%s', mfilename, savename), [], [], 2, [], ...
        'epstopdf');
else
    if fc(2) == fcs(end)
        fc(2) = fNq;
    end
end
end

% compute signal-to-noise ratio for same-length fixed-length segments
% whose breakpoint can vary within a certain time window
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
