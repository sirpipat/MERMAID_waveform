function timedefinitionfigure
% TIMEDEFINITIONFIGURE
%
% Makes a figure that visually defines terms for travel-time measurements
% for the paper
%
% Last modified by sirawich-at-princeton.edu, 04/01/2024

% TODO: Move to $MERMAID2/FIGURE??
obsfile = fullfile(getenv('MERMAID2'), 'DATA', 'Figure9', ...
    '20180817T154351.09_5B77394A.MER.DET.WLT5.sac');
synfile = fullfile(getenv('MERMAID2'), 'DATA', 'Figure9', ...
    '20180817T153501_09_0_SYNTHETIC.sac');
bathdir = fullfile(getenv('MERMAID2'), 'DATA', 'Figure8', 'bathymetry/');

% corner frequency
fcorners = [0.4 2];

% travel time anomaly
t_shift = 4.8501;

% dT = Instaseis - TauP
presidue = -1.6004;

% define color here
BLUE = [0 0.4 1];
GREEN = [0 0.6 0];
RED = [0.8 0 0];
PALEGREEN = [0.7 0.9 0.7];   % for filtered Instasei

[seis_o, hdr_o] = readsac(obsfile);
[seis_s, hdr_s] = readsac(synfile);

[dt_ref_o, ~, ~ , fs_o, ~, dts_o, ~] = gethdrinfo(hdr_o);
[dt_ref_s, ~, ~ , fs_s, ~, dts_s, ~] = gethdrinfo(hdr_s);

[~, ~, t_r, seis_r] = cctransplot(bathdir, bathdir, [], {'bottom', 'displacement'}, ...
    {'hydrophone', 'pressure'}, [], fs_o, false);

% resample to MERMAID datetimes
seis_s = shannon(dts_s, seis_s, dts_o);
tr_r_interp = 0:(1/fs_o):t_r(end);
seis_r = shannon(t_r, seis_r, tr_r_interp);

% convolve
pres_s = conv(seis_s, seis_r) ;
pres_s = pres_s(1:length(seis_o), 1);
pres_s = bandpass(pres_s, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);
pres_o = bandpass(pres_o, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

% filter instaseis for plotting
seis_sf = bandpass(seis_s, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

dt_arrival = dt_ref_o + seconds(hdr_o.T0);
dt_predict = dt_ref_s + seconds(hdr_s.T0);

% time
t = seconds(dts_o - dt_arrival);

%% Make a plot
figure(8)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 10 8])
ax = subplot('Position', [0.09 0.52 0.82 0.42]);

% observed pressure
plot(t, pres_o / max(abs(pres_o)), 'Color', BLUE, 'LineWidth', 2)
hold on
grid on
% mark Joel's pick
scatter(ax, 0, 0, 100, 'k', 'filled');

% synthetic pressure
plot(t + 0 * t_shift, pres_s / max(abs(pres_s)) - 3, 'Color', RED, 'LineWidth', 2)

% instaseis seismogram
plot(t, seis_sf / max(abs(seis_sf)) - 6, 'Color', PALEGREEN, 'LineWidth', 1.5)
plot(t, seis_s / max(abs(seis_s)) - 6, 'Color', GREEN, 'LineWidth', 2)
% mark TauP prediction
yseis_s = interp1(t, seis_s, seconds(dt_predict - dt_arrival), 'linear');
scatter(seconds(dt_predict - dt_arrival), ...
    yseis_s / max(abs(seis_s)) - 6, 100, GREEN, '^', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1)
% mark my picked arrival on Instaseis seismogram
yseis_s = interp1(t, seis_s, ...
    seconds(dt_predict - dt_arrival) + presidue, 'linear');
scatter(seconds(dt_predict - dt_arrival) + presidue, ...
    yseis_s / max(abs(seis_s)) - 6, 100, GREEN, 'v', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1)


xlabel('time since first picked arrival (s)')
nolabels(gca, 2);

xlim([-15 20])
ylim([-8.5 2.5])
yticks([-6 -3 0])

plot([-1.5 1] * t_shift + [0 5 nan nan]', [-4.5 3], 'LineWidth', 0.5, 'Color', [0.8 0.4 0.4], 'LineStyle', '-');

% note arrival time
[~,vpick] = vline(gca, 0, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-');
scatter(-hdr_o.USER4, 0, 140, BLUE, 's', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1);
[~,vtaup] = vline(gca, -hdr_o.USER4, 'LineWidth', 1.5, 'Color', BLUE, 'LineStyle', ':');
vpete = plot(gca, -t_shift * [1 1], [-4.5 -1], 'LineWidth', 1.5, 'Color', [0.6 0 0], 'LineStyle', ':');
vtobs = plot(gca, seconds(dt_predict - dt_arrival) * [1 1], [-9 -4.5], 'LineWidth', 1.5, 'Color', GREEN, 'LineStyle', ':');
viobs = plot(gca, (seconds(dt_predict - dt_arrival) + presidue) * [1 1], [-9 -4.5], 'LineWidth', 1, 'Color', GREEN, 'LineStyle', '--');

% horizontal arrow labels
hpick = arrow(-hdr_o.USER4, 1.5, hdr_o.USER4, 0);
set(hpick, 'LineWidth', 2, 'Color', BLUE, 'MaxHeadSize', 2);
hpete = arrow(-t_shift, -1.5, t_shift, 0, 1.05);
set(hpete, 'LineWidth', 2, 'Color', RED, 'MaxHeadSize', 0.4);
hinst = arrow(seconds(dt_predict - dt_arrival), -7.5, presidue, 0, 1, 1);
set(hinst, 'LineWidth', 2, 'Color', GREEN, 'MaxHeadSize', 4);

% text label
text(-13, 0.7, 'MERMAID', 'Color', BLUE, 'FontSize', 14)
text(-13, -2.3, 'synthetic', 'Color', RED, 'FontSize', 14)
text(-13, -5.3, 'Instaseis', 'Color', GREEN, 'FontSize', 14)
%text(0.5, -4.3, 'picked arrival from Simon et al. (2022)', 'Color', 'k', 'FontSize', 14)
%text(-19.8, 0.6, 'prediction with adjusted ak135 model', 'Color', BLUE, 'FontSize', 14)

%text(0.5, 1.5, '\Delta\tau - Simon et al. 2022', 'Color', BLUE, 'FontSize', 14)
%text(0.5, -1.5, '\Delta\tau - this study', 'Color', [0.6 0 0], 'FontSize', 14)

set(gca, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

% legend axes
axl = subplot('Position', [0.25 0.02 0.60 0.38]);
hold on
% symbol and text labels
x0 = 0.75;
xleft = 0.3;
xright = 1.2;
xtext = 2;

y0 = 1.9;
scatter(x0, y0, 100, 'k', 'filled');
% plot([0 0] + x0, [-0.08 0.08] + y0, 'LineStyle', '-', 'Color', 'k', ...
%     'LineWidth', 1.5)
text(xtext, y0, 'picked arrival from Simon et al. (2022)', ...
    'Color', 'k', 'FontSize', 14)

y0 = 1.7;
scatter(x0, y0, 140, BLUE, 's', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot([0 0] + x0, [-0.08 0.08] + y0, 'LineStyle', '-', 'Color', BLUE, ...
%     'LineWidth', 1.5)
text(xtext, y0, 'TauP prediction adjusted for the water column', ...
    'Color', BLUE, 'FontSize', 14)

y0 = 1.5;
scatter(x0, y0, 100, GREEN, 'v', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot([0 0] + x0, [-0.08 0.08] + y0, 'LineStyle', '-', ...
%     'Color', GREEN, 'LineWidth', 1.5)
text(xtext, y0, 'picked arrival on Instaseis seismogram', ...
    'Color', GREEN, 'FontSize', 14)

y0 = 1.3;
scatter(x0, y0, 100, GREEN, '^', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot([0 0] + x0, [-0.08 0.08] + y0, 'LineStyle', '-', ...
%     'Color', GREEN, 'LineWidth', 1.5)
text(xtext, y0, 'TauP prediction at the ocean bottom', ...
    'Color', GREEN, 'FontSize', 14)

y0 = 1.1;
plot([0 0] + xleft, [-0.08 0.08] + y0, 'LineStyle', ':', 'Color', BLUE, ...
    'LineWidth', 1.5)
plot([0 0] + xright, [-0.08 0.08] + y0, 'LineStyle', '-', 'Color', 'k', ...
    'LineWidth', 1.5)
hpick = arrow(xleft, y0, xright-xleft, 0, 1, 1);
set(hpick, 'LineWidth', 2, 'Color', BLUE, 'MaxHeadSize', 0.5);
text(xtext, y0, 'travel-time anomaly (t_{res}^*) from Simon et al. (2022)', 'Color', BLUE, 'FontSize', 14)

y0 = 0.9;
plot([0 0] + xleft, [-0.08 0.08] + y0, 'LineStyle', ':', ...
    'Color', [0.6 0 0], 'LineWidth', 1.5)
plot([0 0] + xright, [-0.08 0.08] + y0, 'LineStyle', '-', 'Color', 'k', ...
    'LineWidth', 1.5)
hpete = arrow(xleft, y0, xright - xleft, 0, 1, 1);
set(hpete, 'LineWidth', 2, 'Color', RED, 'MaxHeadSize', 0.5);
text(xtext, y0, 'travel-time anomaly (\Delta\tau) in this study', 'Color', RED, 'FontSize', 14)

y0 = 0.7;
plot([0 0] + xleft, [-0.08 0.08] + y0, 'LineStyle', '--', ...
    'Color', GREEN, 'LineWidth', 1)
plot([0 0] + xright, [-0.1 0.1] + y0, 'LineStyle', ':', ...
    'Color', GREEN, 'LineWidth', 1.5)
hpete = arrow(xright, y0, xleft - xright, 0, 1, 1);
set(hpete, 'LineWidth', 2, 'Color', GREEN, 'MaxHeadSize', 0.5);
text(xtext, y0, 'travel-time anomaly adjustment', 'Color', GREEN, 'FontSize', 14)

% sign indication
y0 = 0.5;
text(xleft, y0, 'sign indication', 'Color', 'k', 'FontSize', 14, ...
    'FontWeight', 'bold')

y0 = 0.3;
hplus = arrow(xleft, y0, xright - xleft, 0, 1, 1);
set(hplus, 'LineWidth', 2, 'Color', 'k', 'MaxHeadSize', 0.5);
text(xtext, y0, 'positive', 'Color', 'k', 'FontSize', 14)

y0 = 0.1;
hminus = arrow(xright, y0, xleft - xright, 0, 1, 1);
set(hminus, 'LineWidth', 2, 'Color', 'k', 'MaxHeadSize', 0.5);
text(xtext, y0, 'negative', 'Color', 'k', 'FontSize', 14)

xlim([0 10])
ylim([-0.2 2])

axl.XAxis.Visible = 'off';
axl.YAxis.Visible = 'off';
set(axl, 'Box', 'on')

set(gcf, 'Renderer', 'painters')
figdisp('timedefinitionfigure', [], [], 2, [], 'epstopdf')
end