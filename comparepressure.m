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
% Last modified by sirawich-at-princeton.edu, 11/16/2021

% sampling rate
[~, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);
[~, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);
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
pres_s = conv(seis_s, seis_r) * 100;
pres_s = pres_s(1:length(seis_o), 1);
pres_s = bandpass(pres_s, fs_o, 0.5, 2, 4, 2, 'butter', 'linear');

% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);
pres_o = bandpass(pres_o, fs_o, 0.45, 0.775, 4, 2, 'butter', 'linear');

% compute the cross correlation
[cc, lags] = xcorr(pres_o, pres_s, 'coeff');
lags_time = lags / fs_o;
best_lags = lags(cc == max(cc));
best_lags_time = best_lags / fs_o;

% plot the result
figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [0 6 6 6])

ax1 = subplot('Position', [0.08 0.75 0.86 0.19]);
cla
plot(dts_o', pres_o, 'k')
hold on
plot(dts_o'+seconds(best_lags_time), pres_s, 'r')
grid on
xlim([dts_o(1) dts_o(end)])
legend('observed', 'synthetic')
ylabel('acoustic pressure (Pa)')
title('pressure record')
set(ax1, 'Box', 'on')

ax2 = subplot('Position', [0.08 0.41 0.86 0.19]);
cla
plot(dts_o', pres_o - circshift(pres_s, best_lags))
grid on
xlim([dts_o(1) dts_o(end)])
title('residue')
ylabel('residue (Pa)')
set(ax2, 'Box', 'on')

ax3 = subplot('Position', [0.08 0.08 0.86 0.19]);
cla
plot(lags_time, cc)
grid on
xlim([-10 10])
title('correlation coefficient')
xlabel('timeshift (s)')
ylabel('correlation coefficient')
set(ax3, 'Box', 'on')
end