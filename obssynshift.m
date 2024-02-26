function [t_shift, CCmax, lag, CC, x_o, x_s, dts] = obssynshift(obs, syn)
% [t_shift, CCmax, lag, CC, x_o, x_s, dts] = OBSSYNSHIFT(obs, syn)
%
% Find the optimal time shift between the observed seismogram and
% synthetic seismogram. The seismograms are read from SAC files. The
% synthetic seismogram is interpolated to the same time stamps as the
% observed seismograms using Whittaker-Shannon interpolation.
%
% INPUT:
% obs           full path to a SAC file containging the observed seismogram
% syn           full path to a SAC file containging the synthetic seismogram
%
% OUTPUT:
% t_shift       best time shift
% CCmax         correlation coeficient at the best time shift
% lag           lag times
% CC            correlation coeficient at any lag times
% x_o           observed seismograms used
% x_s           synthetic seismograms used
% dts           datetimes of these two seimograms
%
% SEE ALSO:
% CCSHIFT
%
% Last modified by sirawich-at-princeton.edu, 02/26/2024

badval = -12345;

% read the seismogram
[SeisData_o, HdrData_o, ~, ~, tims_o] = readsac(obs);
[SeisData_s, HdrData_s, ~, ~, tims_s] = readsac(syn);
[dt_ref_o, dt_begin_o, dt_end_o, fs_o, npts_o, dts_o, ~] = gethdrinfo(HdrData_o);
[dt_ref_s, dt_begin_s, dt_end_s, fs_s, npts_s, dts_s, ~] = gethdrinfo(HdrData_s);

if isnan(HdrData_o.T0) || HdrData_o.T0 == badval
    return
end

% remove instrument response + filter the observed seismogram
x_o = real(counts2pa(SeisData_o));
x_o = bandpass(x_o, fs_o, 1, 2, 2, 2, 'butter', 'linear');

% remove the first 5 seconds (containing artefact from filtering)
wh = (tims_o - HdrData_o.B > 5);
x_o = x_o(wh);
dts_o_original = dts_o;
dts = dts_o(wh);

% interpolate the synthetic seismograms to have the same timestamps as the
% observed seismograms
x_s = shannon(dts_s, SeisData_s, dts_o_original);

x_s = bandpass(x_s, fs_o, 1, 2, 2, 2, 'butter', 'linear');
% remove the first 5 seconds (containing artefact from filtering)
x_s = x_s(wh);

% find the envelope
x_o_e = envelope(x_o);
x_s_e = envelope(x_s);

% find best timeshift using envelopes
[CC_e, lag_e] = xcorr(detrend(x_o_e, 1), detrend(x_s_e, 1), 'coeff');
lag_e = lag_e * HdrData_o.DELTA;
[CCmax_e, index] = max(CC_e);
t_shift_e = lag_e(index);
% [t_shift_e, CCmax_e, lag_e, CC_e] = ccshift(x_o_e, x_s_e, dts(1), dts(1), ...
%     fs_o, seconds(50)); % half window size is 100 and I remove the first 10 seconds due to the artifact from filtering

% find best timeshift using original waveforms
[CC, lag] = xcorr(detrend(x_o, 1), detrend(x_s, 1), 'coeff');
lag = lag * HdrData_o.DELTA;
[CCmax, index] = max(CC);
t_shift = lag(index);
% [t_shift, CCmax, lag, CC] = ccshift(x_o, x_s, dts(1), dts(1), ...
%     fs_o, seconds(50)); % half window size is 100 and I remove the first 10 seconds due to the artifact from filtering

%% plot the result
figure(7)
set(gcf, 'Renderer', 'painters', 'Units', 'inches', 'Position', ...
    [2 6 12 8]);
clf
% plot the lag
ax1 = subplot('Position', [0.07 0.72 0.88 0.23], 'FontSize', 12);
plot(lag, CC, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
hold on
grid on
box on
xlabel('lag time (s)')
ylabel('correlation coefficient')
plot(lag_e, CC_e, 'k', 'LineWidth', 1)
scatter(t_shift_e, CCmax_e, 80, 'Marker', 'v', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b');
scatter(t_shift, CCmax, 80, 'Marker', 'v', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'r', 'MarkerEdgeAlpha', 0.5, ...
    'MarkerFaceAlpha', 0.5);
legend('original CC', 'envelope CC', 'Location', 'best')
ax1.XLim = [-100 100];
title(ax1, sprintf('best shift: original = %.2f s, envelope = %.2f s', ...
    t_shift, t_shift_e))

% plot the shifted seimograms using best timeshift from original waveforms
ax2 = subplot('Position', [0.07 0.38 0.88 0.23], 'FontSize', 12);
plot(dts, x_o / max(abs(x_o)), 'k', 'LineWidth', 1)
hold on
grid on
box on
plot(dts + seconds(t_shift), x_s / max(abs(x_s)), 'b', 'LineWidth', 2)
plot(dts,  x_o_e / max(abs(x_o)), 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
plot(dts + seconds(t_shift),  x_s_e / max(abs(x_s)), 'Color', [0.3 0.7 1], 'LineWidth', 1)
vline(ax2, dt_ref_o + seconds(HdrData_o.T0), 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
ax2.XLim = dt_ref_o + seconds(HdrData_o.T0) + seconds([-10 30]);
legend('observed', 'synthetic -- shifted', 'observed envelope', ...
    'synthetic envelope -- shifted', 'Location', 'northeast');
title('applied best timeshift for original waveforms')

% plot the shifted seismogram using best timeshift from envelopes
ax3 = subplot('Position', [0.07 0.05 0.88 0.23], 'FontSize', 12);
plot(dts, x_o / max(abs(x_o)), 'k', 'LineWidth', 1)
hold on
grid on
box on
plot(dts + seconds(t_shift_e), x_s / max(abs(x_s)), 'b', 'LineWidth', 2)
plot(dts,  x_o_e / max(abs(x_o)), 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
plot(dts + seconds(t_shift_e),  x_s_e / max(abs(x_s)), 'Color', [0.3 0.7 1], 'LineWidth', 1)
vline(ax3, dt_ref_o + seconds(HdrData_o.T0), 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
ax3.XLim = dt_ref_o + seconds(HdrData_o.T0) + seconds([-10 30]);
legend('observed', 'synthetic -- shifted', 'observed envelope', ...
    'synthetic envelope -- shifted', 'Location', 'northeast');
title('applied best timeshift for envelope functions')

% save the figure to $EPS
savename = sprintf('%s_%d_%s', mfilename, HdrData_o.USER7, ...
    HdrData_o.KSTNM(ismember(HdrData_o.KSTNM, 33:126)));
% check whether the filename is already taken
if exist([getenv('EPS') savename '.pdf'], 'file')
    % append 01, 02, ... if the filename is already taken
    ii = 1;
    savename_plain = savename;
    savename = sprintf('%s_%02d', savename_plain, ii);
    while and(exist([getenv('EPS') savename '.pdf'], 'file'), ii < 99)
        ii = ii + 1;
        savename = sprintf('%s_%02d', savename_plain, ii);
    end
end
figdisp([savename '.eps'], [], [], 2, [], 'epstopdf');
end