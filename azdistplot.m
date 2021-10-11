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
% SEE ALSO:
% AZIMPROJ
%
% Last modified by sirawich-at-princeton.edu, 10/11/2021

% convert degrees to radians
deg2rad = pi / 180;

figure(1)
clf
% create a polar axes
ax = polaraxes('ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
hold on

% read all sac files and stored in a stream
st_o = cell(length(obs), 1);
for ii = 1:length(obs)
    [SeisData, HdrData, ~, ~, ~] = readsac(obs{ii});
    tr.data = SeisData;
    tr.hdr = HdrData;
    st_o{ii} = tr;
end

% read all synthetic seismograms
st_s = cell(length(syn), 1);
for ii = 1:length(syn)
    [SeisData, HdrData, ~, ~, ~] = readsac(syn{ii});
    tr.data = SeisData;
    tr.hdr = HdrData;
    st_s{ii} = tr;
end

% get metadata for map plotting
eventloc = [HdrData.EVLO HdrData.EVLA];

% plot coastlines
figure(10)      % placeholder for the plotcont (unused)
[~, ~, xy, ~] = plotcont([0 90], [360 -90], 1, 0);
delete(figure(10))
figure(1)
[azDeg, distDeg] = azimproj(eventloc, xy);
polarplot(azDeg * deg2rad, distDeg, 'k')

% plot plates
figure(10)
[~, xy] = plotplates([0 90], [360 -90], 1);
delete(figure(10))
figure(1)
[azDeg, distDeg] = azimproj(eventloc, xy);
polarplot(azDeg * deg2rad, distDeg, 'r')

% plot stations
max_dist = 0;
for ii = 1:length(st_o)
    az = st_o{ii}.hdr.AZ;
    dist = st_o{ii}.hdr.GCARC;
    polarscatter(az * deg2rad, dist, ...
        'Marker', 'v', ...
        'MarkerEdgeColor', rgbcolor('k'), ...
        'MarkerFaceColor', rgbcolor('r'));
    if max_dist < dist
        max_dist = dist;
    end
end

% plot the outer circle of the plot (equivalent to 'box on')
ax.RLim = [0 max_dist*1.1];
polarplot(linspace(0,2*pi,500), linspace(0,2*pi,500).^0 * max_dist * 1.1, 'k')

% plot traces
window_left = -30;
window_right = 30;
axbs = cell(length(st_o), 1);
for ii = 1:length(st_o)
    az = st_o{ii}.hdr.AZ;
    dist = st_o{ii}.hdr.GCARC;
    [dt_ref_o, ~, ~, fs_o, ~, dts_o] = gethdrinfo(st_o{ii}.hdr);
    % remove instrument response
    data_o = real(counts2pa(st_o{ii}.data, fs_o));
    % filter between 1-2 Hz
    data_o = bandpass(data_o, fs_o, 1, 2, 2, 2, 'butter', 'linear');
    % cut and scale
    wh_o = and(dts_o >= dt_ref_o + seconds(st_o{ii}.hdr.T0 + window_left), ...
        dts_o <= dt_ref_o + seconds(st_o{ii}.hdr.T0 + window_right));
    data_o = data_o(wh_o);
    dts_o = dts_o(wh_o);
    data_o = data_o / max(abs(data_o));
    
    [dt_ref_s, ~, ~, fs_s, ~, dts_s] = gethdrinfo(st_s{ii}.hdr);
    data_s = st_s{ii}.data;
    % filter between 1-2 Hz
    data_s = bandpass(data_s, fs_s, 1, 2, 2, 2, 'butter', 'linear');
    % scale and cut
    data_s = data_s / max(abs(data_s));
    wh_s = and(dts_s >= dt_ref_o + seconds(st_o{ii}.hdr.T0 + window_left), ...
        dts_s <= dt_ref_o + seconds(st_o{ii}.hdr.T0 + window_right));
    data_s = data_s(wh_s);
    dts_s = dts_s(wh_s);
    
    % find normalized positions for the trace plot
    x = dist * sin(az * deg2rad);
    y = dist * cos(az * deg2rad);
    
    x_norm = (x + max_dist) / (2 * max_dist);
    y_norm = (y + max_dist) / (2 * max_dist);
    
    axb = addbox(ax, [x_norm-0.05, y_norm + 0.03, 0.1, 0.03]);
    plot(dts_o, data_o, 'Color', 'k')
    hold on
    plot(dts_s, data_s, 'Color', 'r')
    axb.XLim = dt_ref_o + seconds([window_left window_right] + st_o{ii}.hdr.T0);
    axb.Title.String = '';
    axb.XLabel.Visible = 'off';
    axb.XTickLabel = '';
    axb.YTickLabel = '';
    axbs{ii} = axb;
end

for ii = 1:length(axbs)
    axes(axbs{ii});
end


% plot the source
axs = axes('Position', ax.Position);
axs.XLim = [-180 180];
axs.YLim = [-180 180];
addfocalmech(axs, [0 0], 'PublicID', sprintf('%d', HdrData.USER7));
axs.Visible = 'off';

eqid = st_o{1}.hdr.USER7;
ax.Title.String = sprintf('Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', eqid, st_o{1}.hdr.MAG, st_o{1}.hdr.EVDP);

% save the figure
eqid = st_o{1}.hdr.USER7;
figdisp(sprintf('%s_%d.eps', mfilename, eqid), [], [], 2, [], 'epstopdf')
end