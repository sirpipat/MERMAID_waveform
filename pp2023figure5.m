function fig = pp2023figure5
% fig = PP2023FIGURE5
%
% Makes figure 5 of Pipatprathanporn+2023 containing the source-receiver
% map, bathymetry profiles, and the P-wave ray paths.
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 03/20/2023

%% Load data
obsfiles = allfile(sprintf('%sDATA/Figure5/observed/10996154/', ...
    getenv('MERMAID2')))';
synfiles = allfile(sprintf('%sDATA/Figure5/synthetic/10996154/', ...
    getenv('MERMAID2')))';

metadata_o = getheaderarray(obsfiles);
metadata_s = getheaderarray(synfiles);

% ddir = '/Users/sirawich/research/remote_specfem2d/flat_11104502_P0024/';
% name = 'flat_11104502_P0024';
ddir = '/Users/sirawich/research/remote_specfem2d/flat_10996154_P0009/';
name = 'flat_10996154_P0009';

% number of MERMAID floats
n = length(metadata_o.STLO);
i_example = 2;

% MERMAID number
stnm = nan(size(metadata_o.STLO));
for ii = 1:n
    stnm(ii) = str2double(indeks(metadata_o.KSTNM{ii}, '2:end'));
end

%% Plot
figure(22)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 8])
clf

% source-receivers map
ax1 = subplot('Position', [0.04 0.47 0.42 0.51]);

stlo = mod(metadata_o.STLO, 360);
stla = metadata_o.STLA;
evlo = mod(metadata_o.EVLO(1), 360);
evla = metadata_o.EVLA(1);

% map extent
latmin = min(min(stla), evla);
latmax = max(max(stla), evla);
lonmin = min(min(stlo), evlo);
lonmax = max(max(stlo), evlo);

% extend the max extent by 20%
latmid = (latmin + latmax) / 2;
halfheight = (latmax - latmin) / 2;
lonmid = (lonmin + lonmax) / 2;
halfwidth = (lonmax - lonmin) / 2;

latmin = latmid - 1.2 * halfheight;
latmax = latmid + 1.2 * halfheight;
lonmin = lonmid - 1.2 * halfwidth;
lonmax = lonmid + 1.2 * halfwidth;

% zoom in the map
original_x2y_ratio = (ax1.XLim(2)-ax1.XLim(1))/(ax1.YLim(2)-ax1.YLim(1));
new_x2y_ratio = (lonmax-lonmin)/(latmax-latmin);
if new_x2y_ratio > original_x2y_ratio
    latmid = (latmin + latmax) / 2;
    latmin = latmid - (lonmax - lonmin) / original_x2y_ratio / 2;
    latmax = latmid + (lonmax - lonmin) / original_x2y_ratio / 2;
else
    lonmid = (lonmin + lonmax) / 2;
    lonmin = lonmid - (latmax - latmin) * original_x2y_ratio / 2;
    lonmax = lonmid + (latmax - latmin) * original_x2y_ratio / 2;
end

[lons,lats,elev,~,~] = bathymetry([], [lonmin lonmax], [latmin latmax]);
imagesc(lons, lats, elev', [-11000 9000]);
axis xy;

grid on
hold on

[~, cont] = plotcont();
plate = plotplates();
plate.Color = 'r';

% add colorbar
[cb,cm] = cax2dem([-7000 3500], 'hor');
cb.Label.String = 'elevation (m)';
cb.Label.FontSize = 13;
cb.TickDirection = 'both';

% plot trajectory
plottrack(ax1, [evlo evla], [stlo(i_example) stla(i_example)], 0, 101, ...
    'LineWidth', 2, 'Color', csscolor('salmon'));

xlim([lonmin lonmax])
ylim([latmin latmax])

ax1s = doubleaxes(ax1);
axes(ax1s);
scatter(ax1s, mod(stlo, 360), stla, 80, (1:n)', 'filled', 'Marker', 'v', ...
    'MarkerEdgeColor', 'k');

cmap = jet(n);
colormap(ax1s, cmap);
% c = colorbar(ax1s, 'Ticks', 1:n, 'TickLabels', stnm);
% c.Label.String = 'MERMAID number';
% c.Label.FontSize = 12;
% caxis([0.5 n+0.5])
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')
set(ax1, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)

ax1s.XLim = ax1.XLim;
ax1s.YLim = ax1.YLim;
ax1s.Position = [0.04 0.5915 0.42 0.381];
ax1s.XAxis.Visible = 'off';
ax1s.YAxis.Visible = 'off';
set(ax1s, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'Color', 'none')

% beachball diagram
ax1b = doubleaxes(ax1s);
axes(ax1b)
addfocalmech(ax1b, [evlo evla], 'PublicID', ...
    string(metadata_o.USER7(1)), 40, csscolor('salmon'));
ax1b.XAxis.Visible = 'off';
ax1b.YAxis.Visible = 'off';
ax1b.Color = 'none';

% text label
ax1b08 = addbox(ax1s, [0.51 0.4 1/15 1/16]);
text(ax1b08, 0.17, 0.47, '08');
ax1b09 = addbox(ax1s, [0.59 0.56 1/15 1/16]);
text(ax1b09, 0.17, 0.47, '09');
ax1b10 = addbox(ax1s, [0.68 0.38 1/15 1/16]);
text(ax1b10, 0.17, 0.47, '10');
ax1b11 = addbox(ax1s, [0.78 0.49 1/15 1/16]);
text(ax1b11, 0.17, 0.47, '11');
ax1b12 = addbox(ax1s, [0.88 0.37 1/15 1/16]);
text(ax1b12, 0.17, 0.47, '12');

%% bathymetry of all trajectory
ax2 = subplot('Position', [0.54 0.51 0.38 0.47]);
hold on
y_elev = nan(n, 1);
for ii = 1:n
    [x, z] = bathymetryprofile(20000, 401, [stlo(ii) stla(ii)], ...
        mod(metadata_o.BAZ(ii) + 180, 360));
    x = (x - 10000) / 1000;
    z = z - metadata_o.STEL(ii);
    z = z / 2000;
    y_elev(ii) = interp1(x, z, -9) + n + 1 - ii;
    if ii == i_example
        plot(x, z + n + 1 - ii, 'LineWidth', 2, 'Color', 'k');
    else
        plot(x, z + n + 1 - ii, 'LineWidth', 1, 'Color', 'k');
    end
end
grid on
ax2.YTick = 1:n;
% add label
x_label = ones(n,1) * (-9.25);
y_label = y_elev + 0.2;
scatter(ax2, x_label, y_label, 60, (1:n)', 'filled', 'Marker', 'v', ...
    'MarkerEdgeColor', 'k');
colormap(ax2, cmap);

% depth label
scatter(ax2, zeros(8,1), (1:8)', 120, 'k', 'Marker', 'x');
text(ax2, ones(n,1) * (-1.5), n + 1 - (1:n)' - 0.2, ...
    string(-metadata_o.STEL) + ' m');

% distance label
text(ax2, x_label - 0.4, y_label - 0.35, ...
    '\Delta = ' + string(round(metadata_o.GCARC, 2)) + '^{\circ}');

% elevation scale
plot(ax2, [5 5], [2 2.5], 'k', 'LineWidth', 2)
text(ax2, 5.5, 2.3, '2000 m')

ax2.YLim = [0.5 n+0.5];
ax2.YTickLabel = flip(stnm);
ax2.XLabel.String = 'distance along the profile (km)';
ax2.YLabel.String = 'MERMAID number';
set(ax2, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)

ax2s = doubleaxes(ax2);
ax2s.XAxis.Visible = 'off';
ax2s.YTickLabel = round(flip(mod(metadata_o.BAZ + 180, 360)), 0);
ax2s.YLabel.String = 'azimuth (degrees)';
set(ax2s, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, 'Color', 'none')

%% SPECFEM2D setting of the selected path
ax3 = subplot('Position', [0.08 0.07 0.88 0.36]);

% create a figure using drawsetting and then copy to here
drawsetting(ddir, name, [], [], false);
fig = gcf;
copyobj(fig.Children.Children, ax3);
ax3.Children(end-1).FaceColor = [0.3 0.6 0.8];
set(ax3, 'XLim', [0 20000], 'YLim', [0 9600]);
delete(fig);

set(ax3, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, ...
    'DataAspectRatio', [1 1 1])

hold on
b = -metadata_o.STEL(i_example);
b_relfect = 5000;
d = metadata_o.STDP(i_example);
theta_f = asin(metadata_s.USER9(i_example) * 1.5 / (6371 - b_relfect/1000));
theta_s = asin(metadata_s.USER9(i_example) * 3.4 / (6371 - b_relfect/1000));

% plot paths in the ocean layer
xu = [0 d+(-b_relfect)*(1:11)] * tan(theta_f) + 10000;
zu = 9600 - [d 0.5*b_relfect*(1-(-1).^(1:11))];
plot(ax3, xu, zu, 'LineWidth', 2, 'Color', csscolor('lightgreen'));
xd = [0, -d-b_relfect*(0:11)] * tan(theta_f) + 10000;
zd = 9600 - [d 0.5*b_relfect*(1+(-1).^(1:12))];
plot(ax3, xd, zd, 'LineWidth', 2, 'Color', csscolor('salmon'));

% plot paths in the solid earth
xub = xu(2:2:end);
xu0 = xub - 2 * (9600 - b_relfect) * tan(theta_s);
xus = [xub; xu0; nan(size(xub))];
xus = reshape(xus, [1 numel(xus)]);
zus = [zu(2:2:end); zeros(size(zu(2:2:end))); nan(size(zu(2:2:end)))];
zus = reshape(zus, [1 numel(zus)]);
plot(ax3, xus, zus, 'LineWidth', 2, 'Color', csscolor('lightgreen'));
xdb = xd(3:2:end);
xd0 = xdb - 2 * (9600 - b_relfect) * tan(theta_s);
xds = [xdb; xd0; nan(size(xub))];
xds = reshape(xds, [1 numel(xds)]);
zds = zus;
plot(ax3, xds, zds, 'LineWidth', 2, 'Color', csscolor('salmon'));

% plot mermaid icon
mermaid_pole_shade = plot(ax3, [10000 10000], [10120 9080] - d, ...
    'LineWidth', 4.5, 'Color', 'k');
mermaid_pole = plot(ax3, [10000 10000], [10100 9100] - d, ...
    'LineWidth', 3.5, 'Color', csscolor('orange'));
mermaid_ball = scatter(ax3, 10000, 9600 - d, 160, ...
    csscolor('orange'), 'filled', 'Marker', 'o', ...
    'MarkerEdgeColor', 'k');

% OBS icon
obs = scatter(ax3, 10000, 9600 - b, 160, x11color('salmon2'), 'filled', ...
    'Marker', '^', 'LineWidth', 1, 'MarkerFaceColor', [0.75 0.2 0.9], ...
    'MarkerEdgeColor', [0.9 0.9 0.1]);

ax3.XTickLabel = -10:2:10;
ax3.XLabel.String = 'distance along profile (km)';
ax3.YTick = [0 1600:2000:9600];
ax3.YTickLabel = {'9.6', '8.0', '6.0', '4.0', '2.0', '0.0'}';
ax3.YLabel.String = 'depth (km)';

set(gcf, 'Renderer', 'painters')

if nargout > 0
    fig = gcf;
end

%% Save
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');
end