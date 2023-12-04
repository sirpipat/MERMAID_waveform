function fig = pp2023figureA2(obs_struct)
% fig = PP2023FIGUREA2(obs_struct)
%
% Makes figure A2 of Pipatprathanporn+2023 containing the histogram of the
% travel time corrections (pick arrival time in Instaseis seismogram minus
% the ray-thoery based arrival time by TauP for first P-wave arrival) and
% scatter plots of travel time correction against event depth, epicentral
% distance, and event magnitude.
%
% INPUT:
% obs_struct        A struct containing
%   - snr               signal-to-noise ratio
%   - fcorners          [lower upper] corner frequencies
%   - CCmaxs            maximum correlation coefficients for
%                       [flat bath] cases
%   - metadata          SAC Headers associated to the obsfile
%   - presiduals        InstaSeis arrival - TauP prediction for first P
%                       arrival
%
% OUTPUT:
% fig               figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 11/26/2023

fig = figure(2);
set(gcf, 'Units', 'inches', 'Position', [2 2 10 6])
ax1 = subplot('Position', [0.06 0.62 0.42 0.36]);
histogram(obs_struct.presiduals, 'BinWidth', 0.1, 'FaceColor', [1 1 1] * 0.75)
xlim([-4.5 1.5])
ylim([0 160])
xticks(-4:1)
yticks(0:40:160)
xlabel('correction (s)')
ylabel('counts')
grid on
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 0.75, 'Box', 'on')
boxedlabel(ax1, 'northwest', 0.25, [], 'a', 'FontSize', 14)

ax2 = subplot('Position', [0.56 0.62 0.42 0.36]);
scatter(obs_struct.metadata.EVDP, obs_struct.presiduals, 10, 'k', 'filled')
xlim([0 660])
ylim([-5 1])
xticks([0 100 250 410 660])
yticks(-5:1)
xlabel('event depth (km)')
ylabel('correction (s)')
grid on
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 0.75, 'Box', 'on')
boxedlabel(ax2, 'northwest', 0.25, [], 'b', 'FontSize', 14)

ax3 = subplot('Position', [0.06 0.12 0.42 0.36]);
scatter(obs_struct.metadata.GCARC, obs_struct.presiduals, 10, 'k', 'filled')
xlim([0 180])
ylim([-5 1])
xticks(0:30:180)
yticks(-5:1)
xlabel('distance (\circ)')
ylabel('correction (s)')
grid on
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 0.75, 'Box', 'on')
boxedlabel(ax3, 'northwest', 0.25, [], 'c', 'FontSize', 14)

ax4 = subplot('Position', [0.56 0.12 0.42 0.36]);
scatter(obs_struct.cmt.Mw, obs_struct.presiduals, 10, 'k', 'filled')
ylim([-5 1])
yticks(-5:1)
xlabel('magnitude')
ylabel('correction (s)')
grid on
set(gca, 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 0.75, 'Box', 'on')
boxedlabel(ax4, 'northwest', 0.25, [], 'd', 'FontSize', 14)

set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s.eps', mfilename), [], [], 2, [], 'epstopdf');
end