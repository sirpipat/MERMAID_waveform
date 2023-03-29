function varargout = pp2023figure5
% fig = PP2023FIGURE5
%
% Makes figure 5 of Pipatprathanporn+2023 containing the source-receiver
% map, bathymetry profiles, and the P-wave ray paths.
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 03/29/2023

%% Load data
obsfiles = allfile(sprintf('%sDATA/Figure5/observed/10996154/', ...
    getenv('MERMAID2')))';
synfiles = allfile(sprintf('%sDATA/Figure5/synthetic/10996154/', ...
    getenv('MERMAID2')))';

metadata_o = getheaderarray(obsfiles);
metadata_s = getheaderarray(synfiles);

ddir = '/Users/sirawich/research/remote_specfem2d/flat_10996154_P0009/';
name = 'flat_10996154_P0009';

% number of MERMAID floats
n = length(metadata_o.STLO);
i_example = 2;

% constants
vp = 3400;   % P-wave speed in the crust in m/s
vw = 1500;   % acoustic wave speed in water in m/s
dt = 0.5; % wave front spacing in seconds

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
ax1 = subplot('Position', [0.04 0.53 0.42 0.45]);

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
ax1s.DataAspectRatio = [1 1 1];
ax1s.Position = ax1.Position;
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
ax1b.DataAspectRatio = [1 1 1];

% text label
ax1b08 = addbox(ax1s, [0.51 0.4 1/15 1/16]);
text(ax1b08, 0.17, 0.47, '08');
ax1b09 = addbox(ax1s, [0.57 0.56 1/15 1/16]);
text(ax1b09, 0.17, 0.47, '09');
ax1b10 = addbox(ax1s, [0.64 0.38 1/15 1/16]);
text(ax1b10, 0.17, 0.47, '10');
ax1b11 = addbox(ax1s, [0.72 0.49 1/15 1/16]);
text(ax1b11, 0.17, 0.47, '11');
ax1b12 = addbox(ax1s, [0.81 0.37 1/15 1/16]);
text(ax1b12, 0.17, 0.47, '12');

%% bathymetry of all trajectory
ax2 = subplot('Position', [0.54 0.57 0.38 0.41]);
hold on
y_elev_left = nan(n, 1);
y_elev_right = nan(n, 1);
for ii = 1:n
    [x, z] = bathymetryprofile(20000, 401, [stlo(ii) stla(ii)], ...
        mod(metadata_o.BAZ(ii) + 180, 360));
    x = (x - 10000) / 1000;
    z = z - metadata_o.STEL(ii);
    z = z / 2000;
    y_elev_left(ii) = interp1(x, z, -9) + n + 1 - ii;
    y_elev_right(ii) = interp1(x, z, 9) + n + 1 - ii;
    if ii == i_example
        plot(x, z + n + 1 - ii, 'LineWidth', 2, 'Color', 'k');
    else
        plot(x, z + n + 1 - ii, 'LineWidth', 1, 'Color', 'k');
    end
end
grid on
ax2.YTick = 1:n;
% add label
x_label_left = ones(n,1) * (-9.25);
y_label_left = y_elev_left + 0.2;
x_label_right = ones(n,1) * 9.25;
y_label_right = y_elev_right + 0.2;
scatter(ax2, x_label_right, y_label_right, 60, (1:n)', 'filled', ...
    'Marker', 'v', 'MarkerEdgeColor', 'k');
colormap(ax2, cmap);

% depth label
scatter(ax2, zeros(8,1), (1:8)', 120, 'k', 'Marker', 'x');
text(ax2, ones(n,1) * (-1.5), n + 1 - (1:n)' - 0.2, ...
    string(-metadata_o.STEL) + ' m');

% distance label
text(ax2, x_label_left - 0.4, y_label_left - 0.35, ...
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
ax3 = subplot('Position', [0.08 0.07 0.88 0.42]);

% draw the seafloor setting
interfacefile = [ddir 'DATA/interfaces_' name '.dat'];
ax3 = drawbackground(interfacefile, ax3);
ax3.Children(end-1).FaceColor = [0.3 0.6 0.8];
set(ax3, 'XLim', [0 20000], 'YLim', [0 9600]);

set(ax3, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12, ...
    'DataAspectRatio', [1 1 1])

hold on
b = -metadata_o.STEL(i_example);
b_relfect = 5000;
d = metadata_o.STDP(i_example);
theta_f = asin(metadata_s.USER9(i_example) * 1.5 / (6371 - b_relfect/1000));
theta_s = asin(metadata_s.USER9(i_example) * 3.4 / (6371 - b_relfect/1000));
[x, z] = bathymetryprofile(20000, 401, [stlo(i_example) stla(i_example)], ...
        mod(metadata_o.BAZ(i_example) + 180, 360));

% plot paths in the ocean layer
xu = [0 d+(-b_relfect)*(1:11)] * tan(theta_f) + 10000;
bu_profile = interp1(x, -z, xu(2:end));
zu = 9600 - [d 0.5*bu_profile.*(1-(-1).^(1:11))];
plot(ax3, xu, zu, 'LineWidth', 2, 'Color', csscolor('lightgreen'));
xd = [0, -d-b_relfect*(0:11)] * tan(theta_f) + 10000;
bd_profile = interp1(x, -z, xd(2:end));
zd = 9600 - [d 0.5*bd_profile.*(1+(-1).^(1:12))];
plot(ax3, xd, zd, 'LineWidth', 2, 'Color', csscolor('salmon'));

% plot paths in the solid earth
xub = xu(2:2:end);
xu0 = xub - zu(2:2:end) * tan(theta_s);
xus = [xub; xu0; nan(size(xub))];
xus = reshape(xus, [1 numel(xus)]);
zus = [zu(2:2:end); zeros(size(zu(2:2:end))); nan(size(zu(2:2:end)))];
zus = reshape(zus, [1 numel(zus)]);
plot(ax3, xus, zus, 'LineWidth', 2, 'Color', csscolor('lightgreen'));
xdb = xd(3:2:end);
xd0 = xdb - zd(3:2:end) * tan(theta_s);
xds = [xdb; xd0; nan(size(xdb))];
xds = reshape(xds, [1 numel(xds)]);
zds = [zd(3:2:end); zeros(size(zd(3:2:end))); nan(size(zd(3:2:end)))];
zds = reshape(zds, [1 numel(zds)]);
plot(ax3, xds, zds, 'LineWidth', 2, 'Color', csscolor('salmon'));

% wavefronts in the solid earth
xw0 = 1000:(vp*dt/sin(theta_s)):((9600+z(end)) / tan(theta_s) + 20000);
slope = -tan(theta_s);
xws = [];
zws = [];
xwb = [];
zwb = [];
for ii = 1:length(xw0)
    p = pierce(slope, [xw0(ii) 0], [x 9600+z], 'linear', 0.01);
    if isempty(p)
        xws = [xws; xw0(ii); 0; nan];
        zws = [zws; 0; xw0(ii) * tan(theta_s); nan];
    else
        xws = [xws; xw0(ii); p(end,1); nan];
        zws = [zws; 0; p(end,2); nan];
        xwb = [xwb; p(1,1)];
        zwb = [zwb; p(1,2)];
    end
end
plot(ax3, xws, zws, 'LineWidth', 1, 'Color', csscolor('white'));

% wavefront in the water
xw = [];
zw = [];
for ii = 1:length(xwb)
    p = pierce(-tan(theta_f), [xwb(ii) zwb(ii)], [x 9600+z], 'linear', 0.01);
    if or(size(p, 1) <= 1, p(1, 1) >= xwb(ii))
        xw = [xw; xwb(ii); 0; nan];
        zw = [zw; zwb(ii); zwb(ii) + tan(theta_f) * xwb(ii); nan];
    else
        % find where the wave interesct the interface to the left
        [xin, iin] = max(p((p(:,1) < xwb(ii)),1));
        xw = [xw; xwb(ii); xin; nan];
        zw = [zw; zwb(ii); p(iin,2); nan];
    end
end
zw_more_left = zw(end-1):(vw*dt/cos(theta_f)):(20000 * tan(theta_f) + 9600);
zw_more_left = zw_more_left(2:end)';
xw_more_left = zeros(size(zw_more_left));
zw_more_right = zw_more_left - 20000 * tan(theta_f);
xw_more_right = ones(size(zw_more_right)) * 20000;
zw_more = [zw_more_left zw_more_right nan(size(zw_more_left))];
zw_more = reshape(zw_more', [numel(zw_more) 1]);
xw_more = [xw_more_left xw_more_right nan(size(xw_more_left))];
xw_more = reshape(xw_more', [numel(xw_more) 1]);
xw = [xw; xw_more];
zw = [zw; zw_more];
plot(ax3, xw, zw, 'LineWidth', 1, 'Color', csscolor('white'));

% plot mermaid icon
mermaid = scatter(ax3, 10000, 9600 - d, 160, ...
    cmap(i_example), 'filled', 'Marker', 'v', ...
    'MarkerEdgeColor', 'k');

% OBS icon
obs = scatter(ax3, 10000, 9600 - b, 160, x11color('salmon2'), 'filled', ...
    'Marker', '^', 'LineWidth', 1, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k');

ax3.XTickLabel = -10:2:10;
ax3.XLabel.String = 'distance along profile (km)';
ax3.YTick = [0 1600:2000:9600];
ax3.YTickLabel = {'9.6', '8.0', '6.0', '4.0', '2.0', '0.0'}';
ax3.YLabel.String = 'depth (km)';

set(gcf, 'Renderer', 'painters')

if nargout > 0
    varargout{1} = gcf;
end

%% Save
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');
end