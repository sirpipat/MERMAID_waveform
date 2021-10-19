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
% Last modified by sirawich-at-princeton.edu, 10/19/2021

% read the seismogram
[SeisData_o, HdrData_o, ~, ~, tims_o] = readsac(obs);
[SeisData_s, HdrData_s, ~, ~, tims_s] = readsac(syn);
[dt_ref_o, dt_begin_o, dt_end_o, fs_o, npts_o, dts_o, ~] = gethdrinfo(HdrData_o);
[dt_ref_s, dt_begin_s, dt_end_s, fs_s, npts_s, dts_s, ~] = gethdrinfo(HdrData_s);

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

% find best timeshift
[t_shift, CCmax, lag, CC] = ccshift(x_o_e, x_s_e, dts(1), dts(1), ...
    fs_o, seconds(50)); % half window size is 100 and I remove the first 10 seconds due to the artifact from filtering

figure(7)
clf
% plot the lag
ax1 = subplot(2, 1, 1);
plot(lag, CC, 'k', 'LineWidth', 1)
hold on
grid on
xlabel('lag time (s)')
ylabel('correlation coefficient')
scatter(t_shift, CCmax, 80, 'Marker', 'v', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b');

% plot the seismogram
ax2 = subplot(2, 1, 2);
plot(dts, x_o / max(abs(x_o)), 'k', 'LineWidth', 1)
hold on
grid on
plot(dts + seconds(t_shift), x_s / max(abs(x_s)), 'b', 'LineWidth', 2)
plot(dts,  x_o_e / max(abs(x_o)), 'Color', [0.7 0.7 0.7], 'LineWidth', 1)
plot(dts + seconds(t_shift),  x_s_e / max(abs(x_s)), 'Color', [0.3 0.7 1], 'LineWidth', 1)
vline(ax2, dt_ref_o + seconds(HdrData_o.T0), 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
end