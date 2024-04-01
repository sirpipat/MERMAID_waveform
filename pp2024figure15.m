function varargout = pp2024figure15(obs_struct)
% fig = PP2024FIGURE15
%
% Makes figure 15 of Pipatprathanporn+2024
%
% INPUT:
% obs_struct        A struct containing
%   - snr               signal-to-noise ratio
%   - fcorners          [lower upper] corner frequencies
%   - CCmaxs            maximum correlation coefficients for
%                       [flat bath] cases
%   - metadata          SAC Headers associated to the obsfile
%   - presiduals        InstaSeis arrival - TauP prediction for first P
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 04/01/2024

%% Define plot parameters here
axmargin = [0.10 0.20 0.08 0.10];
figmargin = [0.03 0.03 0.03 0.01];

COLOR_NAVY  = [45 114 183] / 255;
COLOR_GRAY  = [0.7 0.7 0.7];
COLOR_GREEN = [58 118 119]/ 255;

%% Calculations for derived parameters
% travel time
ttravel =  obs_struct.metadata.USER5 + obs_struct.metadata.USER6;
% relative travel time residual
percent_dlnt = obs_struct.t_shifts(:,2) ./ ttravel * 100;

% occupied bandwidth
bandwidth = obs_struct.fcorners(:,2) - obs_struct.fcorners(:,1);

% largest acceptable error for value comparisons
epsilon = 1e-6;
% corner frequencies bins
fc_bin_lower = (0.4:0.05:1.95)';
fc_bin_upper = fc_bin_lower + 0.05;
fc_bin_mid = (fc_bin_lower + fc_bin_upper) / 2;
% number of traces that occupy the frequency bands
OB = zeros(size(fc_bin_mid, 1), 1);

for ii = 1:length(obs_struct.fcorners(:,1))
    wh = and(fc_bin_lower - obs_struct.fcorners(ii, 1) >= -epsilon, ...
        obs_struct.fcorners(ii, 2) - fc_bin_upper >= -epsilon);
    OB(wh) = OB(wh) + 1;
end

%% Making figure
fig = figure(15);
cla
set(gcf, 'Units', 'inches', 'Position', [0 1 9 7])

ax1 = subplot('Position', subplotposition(3, 3, 1, axmargin, figmargin));
histogram(obs_struct.t_shifts(:,2), 'BinWidth', 2, ...
    'FaceColor', COLOR_NAVY, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xlabel('travel-time residual (s)')
set(ax1, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
boxedlabel(ax1, 'northwest', 0.28, 'inches', 'a', 'FontSize', 15);

ax2 = subplot('Position', subplotposition(3, 3, 2, axmargin, figmargin));
histogram(log10(obs_struct.snr(:,2)), 'BinWidth', 0.2, ...
    'FaceColor', COLOR_GRAY, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
ylim([0 210])
yticks([0 70 140 210])
xlabel('log(signal-to-noise ratio)')
set(ax2, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
axb2 = boxedlabel(ax2, 'northwest', 0.28, 'inches', 'b', 'FontSize', 15);
set(axb2.Children, 'Position', [0.5 0.5 0])

ax3 = subplot('Position', subplotposition(3, 3, 3, axmargin, figmargin));
histogram(obs_struct.CCmaxs(:,2), 'BinWidth', 0.05, ...
    'FaceColor', COLOR_GRAY, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xlabel('cross-correlation')
set(ax3, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
boxedlabel(ax3, 'northwest', 0.28, 'inches', 'c', 'FontSize', 15);

ax4 = subplot('Position', subplotposition(3, 3, 4, axmargin, figmargin));
histogram(percent_dlnt, 'BinWidth', 0.5, ...
    'FaceColor', COLOR_NAVY, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xlim([-5 5])
xticks([-4 -2 0 2 4])
xlabel('travel-time residual (%)')
set(ax4, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
axb4 = boxedlabel(ax4, 'northwest', 0.28, 'inches', 'd', 'FontSize', 15);
set(axb4.Children, 'Position', [0.5 0.5 0])

ax5 = subplot('Position', subplotposition(3, 3, 5, axmargin, figmargin));
histogram(bandwidth, (0.475:0.05:1.625), ...
    'FaceColor', COLOR_GRAY, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xticks([0.5 1 1.6])
xlabel('bandwidth (Hz)')
set(ax5, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
boxedlabel(ax5, 'northwest', 0.28, 'inches', 'e', 'FontSize', 15);

ax6 = subplot('Position', subplotposition(3, 3, 6, axmargin, figmargin));
bar(fc_bin_mid, OB, 1, 'FaceColor', COLOR_GRAY, 'LineWidth', 1)
grid on
xlim([0.3 2.1])
xticks([0.4 0.8 1.2 1.6 2])
ylim([0 1200])
yticks([0 400 800 1200])
xlabel('frequency (Hz)')
set(ax6, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
axb6 = boxedlabel(ax6, 'northwest', 0.28, 'inches', 'f', 'FontSize', 15);
set(axb6.Children, 'Position', [0.5 0.5 0])

ax7 = subplot('Position', subplotposition(3, 3, 7, axmargin, figmargin));
histogram(obs_struct.metadata.GCARC, 'BinWidth', 7.5, ...
    'FaceColor', COLOR_GREEN, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xlim([-10 190])
xticks([0 45 90 135 180])
ylim([0 210])
yticks([0 70 140 210])
xlabel('distance (�)')
set(ax7, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
boxedlabel(ax7, 'northeast', 0.28, 'inches', 'g', 'FontSize', 15);

ax8 = subplot('Position', subplotposition(3, 3, 8, axmargin, figmargin));
histogram(obs_struct.metadata.EVDP, 'BinWidth', 33, ...
    'FaceColor', COLOR_GREEN, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xlim([-30 700])
xticks([0 100 250 410 660])
ylim([0 450])
yticks([0 150 300 450])
xlabel('event depth (km)')
set(ax8, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
axb8 = boxedlabel(ax8, 'northeast', 0.28, 'inches', 'h', 'FontSize', 15);
set(axb8.Children, 'Position', [0.5 0.5 0])

ax9 = subplot('Position', subplotposition(3, 3, 9, axmargin, figmargin));
histogram(obs_struct.metadata.BAZ, 'BinWidth', 18, ...
    'FaceColor', COLOR_GREEN, 'FaceAlpha', 1, 'LineWidth', 1)
grid on
xlim([-20 380])
xticks((0:90:360)')
ylim([0 240])
yticks([0 80 160 240])
xlabel('backazimuth (�)')
set(ax9, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1)
axb9 = boxedlabel(ax9, 'northwest', 0.28, 'inches', 'i', 'FontSize', 15);
set(axb9.Children, 'Position', [0.5 0.5 0])

set(gcf, 'Renderer', 'painters')
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');

if nargout > 0
    varargout{1} = fig;
end
end