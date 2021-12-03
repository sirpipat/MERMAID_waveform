function comparepressure(seis_s, hdr_s, seis_o, hdr_o, seis_r, t_r)
% COMPAREPRESSURE(seis_s, hdr_s, seis_o, hdr_o, seis_r, t_r)
%
% interpolates synthetic and response. Then, convolves the two and compares
% with the observed with instrument response removed.
%
% INPUT:
% seis_s        synthetic seismogram
% hdr_s         SAC header of the synthetic seismogram
% seis_o        observed seismogram (pressure record at hydrophone)
% hdr_o         SAC header of the observed seismogram
% seis_r        receiver function
% t_r           time of the receiver function
%
% OUTPUT:
%
% Last modified by sirawich-at-princeton.edu, 12/03/2021

% sampling rate
[~, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);
[dt_ref_o, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);
fs_r = 1 / (t_r(2) - t_r(1));

% downsample the response to about 100 Hz
d_factor = ceil(fs_r / 100);
if d_factor > 1
    seis_r = lowpass(seis_r, fs_r, fs_r / d_factor, 2, 2, 'butter', 'linear');
    seis_r = downsample(seis_r, d_factor);
    t_r = downsample(t_r, d_factor);
    fs_r = fs_r / d_factor;
end

% resample to MERMAID datetimes
seis_s = shannon(dts_s, seis_s, dts_o);
tr_r_interp = 0:(1/fs_o):t_r(end);
seis_r = shannon(t_r, seis_r, tr_r_interp);

% convolve
pres_s = conv(seis_s, seis_r) ;
pres_s = pres_s(1:length(seis_o), 1);
pres_s = bandpass(pres_s, fs_o, 0.5, 2, 4, 2, 'butter', 'linear');

% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);
pres_o = bandpass(pres_o, fs_o, 0.5, 2, 4, 2, 'butter', 'linear');

% compute the cross correlation
% cut the short windows to do cross correlation
best_lags_time_e = 0;
ep = (1 / fs_o) / 100;
dt_start = dt_ref_o + seconds(hdr_o.T0 - 10);
dt_end = dt_ref_o + seconds(hdr_o.T0 + 20);
pres_o1 = pres_o(and(geq(dts_o, dt_start, ep), leq(dts_o, dt_end, ep)));
pres_s1 = pres_s(and(geq(dts_o + seconds(best_lags_time_e), dt_start, ep), ...
    leq(dts_o + seconds(best_lags_time_e), dt_end, ep)));

% first correlate by the envelope
[best_lags_time_e, ~, ~, ~, ~] = ccscale(pres_o1, pres_s1, ...
    dt_begin_o, dt_begin_o, fs_o, seconds(15), true);

% then correlate by the waveform
[best_lags_time, ~, lags_time, cc, s] = ccscale(pres_o1, pres_s1, ...
    dt_start, dt_start + seconds(best_lags_time_e), fs_o, seconds(5), false);

best_lags_time = best_lags_time + best_lags_time_e;
lags_time = lags_time + best_lags_time_e;

best_lags = round(best_lags_time * fs_o);
% [cc, lags] = xcorr(pres_o, pres_s, 'coeff');
% lags_time = lags / fs_o;
% best_lags = lags(cc == max(cc));
% best_lags_time = best_lags / fs_o;

% plot the result
figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [0 6 6 6])

% plot two pressure records: observed vs synthetic
ax1 = subplot('Position', [0.08 0.68 0.86 0.24]);
cla
plot(dts_o', pres_o, 'k')
hold on
plot(dts_o'+seconds(best_lags_time), s * pres_s, 'b', 'LineWidth', 1)
grid on
xlim(dt_ref_o + seconds(hdr_o.T0 + [-10 40]))
ylim([-1.1 1.1] * max(max(abs(pres_o)), max(abs(s * pres_s))))
vline(ax1, dt_ref_o + seconds(hdr_o.T0), 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
legend('observed', 'synthetic')
ylabel('acoustic pressure (Pa)')
title('pressure record')
set(ax1, 'Box', 'on')

% residue: observed - shifted synthetic
ax2 = subplot('Position', [0.08 0.34 0.86 0.24]);
cla
plot(dts_o', pres_o, 'k')
hold on
plot(dts_o'+seconds(best_lags_time), s * pres_s, 'b', 'LineWidth', 1)
grid on
xlim([dts_o(1) dts_o(end)])
ylim([-1.1 1.1] * max(max(abs(pres_o)), max(abs(s * pres_s))))
vline(ax2, dt_ref_o + seconds(hdr_o.T0), 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
legend('observed', 'synthetic')
ylabel('acoustic pressure (Pa)')
title('pressure record')
set(ax2, 'Box', 'on')
% plot(dts_o', pres_o - circshift(s * pres_s, best_lags), 'k')
% grid on
% xlim(dt_ref_o + seconds(hdr_o.T0 + [-10 40]))
% ylim([-1.1 1.1] * max(abs(pres_o - circshift(pres_s, best_lags))));
% vline(ax2, dt_ref_o + seconds(hdr_o.T0), 'LineWidth', 2, ...
%     'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
% title('residue')
% ylabel('residue (Pa)')
% set(ax2, 'Box', 'on')

% correlation coefficient
ax3 = subplot('Position', [0.08 0.08 0.86 0.17]);
cla
plot(lags_time, cc, 'k')
hold on
scatter(best_lags_time, max(cc), 80, 'Marker', 'v', 'MarkerEdgeColor', ...
    'k', 'MarkerFaceColor', 'b');
grid on
xlim([-20 20])
title('correlation coefficient')
xlabel('timeshift (s)')
ylabel('correlation coefficient')
set(ax3, 'Box', 'on')

% save figure
set(gcf, 'Renderer', 'painters')
savename = sprintf('%s_%d_%s.eps', mfilename, hdr_o.USER7, ...
    replace(hdr_o.KSTNM, ' ', ''));
figdisp(savename, [], [], 2, [], 'epstopdf');
end

function r = leq(a, b, ep)
r = or(a < b, abs(a - b) < ep);
end

function r = geq(a, b, ep)
r = or(a > b, abs(a - b) < ep);
end