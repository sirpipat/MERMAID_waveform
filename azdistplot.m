function [ax, axs] = azdistplot(obs, syn)
% [ax, axs] = AZDISTPLOT(obs, syn)
%
% Plot a source-receiver map with the source at a center of the polar plot.
%
% INPUT:
% obs           observed seismogram sac files
% syn           synthetic seismogram sac files
%
% OUTPUT:
% ax            axes handle to the plot
% axs           axes handle of the beachball
%
% Last modified by sirawich-at-princeton.edu, 10/07/2021

% convert degrees to radians
deg2rad = pi / 180;

figure(1)
clf
% create a polar axes
ax = polaraxes('ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
hold on

% get metadata for map plotting
[~, HdrData, ~, ~, ~] = readsac(obs{1});
eventloc = [HdrData.EVLO HdrData.EVLA];

% plot coastlines
figure(10)      % placeholder for the plotcont (unused)
[~, ~, xy, ~] = plotcont([0 90], [360 -90], 1, 0);
delete(figure(10))
figure(1)
distDeg = zeros(size(xy, 1), 1);
azDeg = zeros(size(xy, 1), 1);
for ii = 1:length(distDeg)          
    [~,distDeg(ii)] = grcdist(xy(ii,:), eventloc);
    [azDeg(ii),~] = azim(eventloc, xy(ii,:));  
end
polarplot(azDeg * deg2rad, distDeg, 'k')

% plot plates
figure(10)
[~, xy] = plotplates([0 90], [360 -90], 1);
delete(figure(10))
figure(1)
distDeg = zeros(size(xy, 1), 1);
azDeg = zeros(size(xy, 1), 1);
for ii = 1:length(distDeg)          
    [~,distDeg(ii)] = grcdist(xy(ii,:), eventloc);
    [azDeg(ii),~] = azim(eventloc, xy(ii,:));  
end
polarplot(azDeg * deg2rad, distDeg, 'r')

% plot stations
max_dist = 0;
for ii = 1:length(obs)
    [SeisData, HdrData, ~, ~, tims] = readsac(obs{ii});
    polarscatter(HdrData.AZ * deg2rad, HdrData.GCARC, ...
        'Marker', 'v', ...
        'MarkerEdgeColor', rgbcolor('k'), ...
        'MarkerFaceColor', rgbcolor('r'));
    if max_dist < HdrData.GCARC
        max_dist = HdrData.GCARC;
    end
end


% plot the source
axs = axes('Position', ax.Position);
axs.XLim = [-180 180];
axs.YLim = [-180 180];
addfocalmech(axs, [0 0], 'PublicID', sprintf('%d', HdrData.USER7));
axs.Visible = 'off';
axes(ax)

% plot the outer circle of the plot (equivalent to 'box on')
ax.RLim = [0 max_dist*1.1];
polarplot(linspace(0,2*pi,500), linspace(0,2*pi,500).^0 * max_dist * 1.1, 'k')
axes(axs)
end