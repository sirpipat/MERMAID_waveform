function plotpaths(sacfiles)
% PLOTPATHS(sacfiles)
%
% Plot source-receivers pairs from all SACFILES. The sources and receivers
% are marked as yellow stars and orange upside-down triangles,
% respectively. Each source-receiver pair is connected by a light-blue
% line.
%
% INPUT
% sacfiles          cell array of sac file names
%
% Last modified by sirawich-at-princeton.edu, 05/16/2024

BADVAL = -12345;

figure
set(gcf, 'Unit', 'inches', 'Position', [0 0 12 8]);
ax = subplot(1,1,1);

% plot coastlines
[axlim,handl,XYZ] = plotcont([0 90], [360 -90], 1, 0);
% plot plate boundaries
[handlp, XYp] = plotplates([0 90], [360 -90], 1);
handlp.Color = 'r';

stlos = zeros(size(sacfiles));
stlas = zeros(size(sacfiles));
evlos = zeros(size(sacfiles));
evlas = zeros(size(sacfiles));

for ii = 1:length(sacfiles)
    [~, hdr] = readsac(sacfiles{ii});
    stlos(ii) = mod(hdr.STLO, 360);
    stlas(ii) = hdr.STLA;
    evlos(ii) = mod(hdr.EVLO, 360);
    evlas(ii) = hdr.EVLA;
end

% remove invalid entries
wh = and(and(stlos ~= BADVAL, stlas ~= BADVAL), ...
    and(evlos ~= BADVAL, evlas ~= BADVAL));

stlos = stlos(wh);
stlas = stlas(wh);
evlos = evlos(wh);
evlas = evlas(wh);

n = length(stlos);

% plot paths
for ii = 1:n
    plottrack(ax, [evlos(ii) evlas(ii)], [stlos(ii) stlas(ii)], 0, ...
          100, 'LineWidth', 0.25, 'Color', [0.7 0.8 0.9]);
end


% plot sources 
scatter(evlos, evlas, 80, 'Marker', 'p', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y')

% plot stations
scatter(stlos, stlas, 40, 'Marker', 'v', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.4], ...
    'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 1)

grid on
xlim([0 360]);
ylim([-90 90]);

% add box
ax.Box = 'on';

% ticks label
ax.XTick = 0:30:360;
ax.XTickLabel = {'0', '30', '60', '90', '120', '150', '180', '-150', ...
                  '-120', '-90', '-60', '-30', '0'};
ax.YTick = -90:30:90;
ax.YTickLabel = {'-90', '-60', '-30', '0', '30', '60', '90'};
ax.TickDir = 'both';
axs = doubleaxes(ax);
axs.Position = ax.Position;
axs.DataAspectRatio = ax.DataAspectRatio;
axs.TickDir = 'both';

% title
n_events = length(unique(1000*evlas+evlos));
title(axs, sprintf('Number of events = %d', n_events))

% save figure
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s.eps', mfilename), [], [], 2, [], 'epstopdf');
end