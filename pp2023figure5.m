function pp2023figure5
% Make figure 5 of Pipatprathanporn+2023

%% TODO: Handle optional inputs

%% Load data
% TODO: Copy these files to $MERMAID2/IFILES
obsfiles = ...
    {'/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155349.08_5D627726.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155458.09_5D6AA390.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155515.10_5D6AAB9E.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155557.11_5D61B0DE.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155628.12_5D61AD96.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155905.16_5D61ACD6.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155758.24_5D62D93C.MER.DET.WLT5.sac', ...
     '/Users/sirawich/research/processed_data/MERMAID_reports_updated/11104502/20190824T155757.25_5D61AC08.MER.DET.WLT5.sac'};
synfiles = ...
    {'/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_08_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_09_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_10_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_11_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_12_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_16_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_24_0_SYNTHETIC.sac', ...
     '/Users/sirawich/research/processed_data/SYNTHETICS/11104502/20190824T155127_25_0_SYNTHETIC.sac'};

[seis_o, hdr_o, ~, ~, tims_o] = readsac(obsfiles{end-1});
[dt_ref_o, dt_begin_o, ~, fs_o] = gethdrinfo(hdr_o);

metadata_o = getheaderarray(obsfiles);
metadata_s = getheaderarray(synfiles);

% number of MERMAID floats
n = length(metadata_o.STLO);

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
ax1 = subplot('Position', [0.08 0.43 0.38 0.51]);

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
ax1s.Position = [0.08 0.5515 0.38 0.381];
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

% bathymetry of all trajectory
ax2 = subplot('Position', [0.58 0.43 0.38 0.51]);
set(ax2, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)

% SPECFEM2D setting of the selected path
ax3 = subplot('Position', [0.08 0.08 0.88 0.23]);
set(ax3, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)

set(gcf, 'Renderer', 'painters')
%% Save
figname = sprintf('%s.eps', mfilename);
%figdisp(figname, [], [], 2, [], 'epstopdf');
end