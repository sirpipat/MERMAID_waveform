function plotoccupiedbandwidth(obs_struct)
% PLOTOCCUPIEDBANDWIDTH(obs_struct)
%
% Plots occupied bandwidths of all traces given the chosen corner
% frequencies.
%
% INPUT:
% obs_struct        A struct containing
%   - snr               signal-to-noise ratio
%   - fcorners          [lower upper] corner frequencies
%   - CCmaxs            maximum correlation coefficients for
%                       [flat bath] cases
%   - metadata          SAC Headers associated to the obsfile
%
% Last modified by sirawich-at-princeton.edu: 10/30/2023

% largest acceptable error for value comparisons
epsilon = 1e-6;

% corner frequencies bins
fc_bin_lower = 0.4:0.05:1.95;
fc_bin_upper = fc_bin_lower + 0.05;
fc_bin_mid = (fc_bin_lower + fc_bin_upper) / 2;

% construct the containers for occupied bandwidth info
OB_snr = nan(size(obs_struct.fcorners,1), size(fc_bin_mid,2));
OB_cc  = nan(size(obs_struct.fcorners,1), size(fc_bin_mid,2));

for ii = 1:size(obs_struct.fcorners,1)
    wh = and(fc_bin_lower - obs_struct.fcorners(ii, 1) >= -epsilon, ...
        obs_struct.fcorners(ii, 2) - fc_bin_upper >= -epsilon);
    OB_snr(ii, wh) = obs_struct.snr(ii);
    OB_cc(ii, wh)  = obs_struct.CCmaxs(ii, 2);
end

% what to sort
% 1. no sort
% 2. SNR
% 3. CC
% 4. magnitude
% 5. epicentral distance
% 6. event depth
% 7. multiple sorting criteria

bandwidth = obs_struct.fcorners(:,2) - obs_struct.fcorners(:,1);
lowcorner = obs_struct.fcorners(:,1);
midband = (obs_struct.fcorners(:,2) + obs_struct.fcorners(:,1)) / 2;
highcorner = obs_struct.fcorners(:,2);
[~,~,BAZ_bin] = histcounts(obs_struct.metadata.BAZ, 'BinWidth', 45);
[~,~,MAG_bin] = histcounts(obs_struct.metadata.MAG, 'BinWidth', 1);
snr = obs_struct.snr;
cc = obs_struct.CCmaxs(:,2);

polarity = nan(size(obs_struct.snr));
for ii = 1:length(obs_struct.snr)
    try
        cmtp = cmtpolarity([obs_struct.cmt.Mrr(ii), ...
            obs_struct.cmt.Mtt(ii), obs_struct.cmt.Mpp(ii), ...
            obs_struct.cmt.Mrt(ii), obs_struct.cmt.Mrp(ii), ...
            obs_struct.cmt.Mtp(ii)], obs_struct.cmt.Dep(ii), ...
            obs_struct.metadata.AZ(ii), obs_struct.metadata.USER9(ii), ...
            'ak135', strcmp('P', indeks(obs_struct.metadata.KT0{ii})));
        if abs(imag(cmtp)) < 1e-8
            polarity(ii) = cmtp;
        else
            polarity(ii) = nan;
        end
    catch
        polarity(ii) = nan;
        continue
    end
end
[~,~,POL_bin] = histcounts(polarity, [-1.0 -0.6 -0.2 0.2 0.6 1.0]);

sorting_criteria = [bandwidth, ...
                    lowcorner, ...
                    midband, ...
                    highcorner, ...
                    BAZ_bin, ...
                    MAG_bin, ...
                    POL_bin, ...
                    snr, ...
                    cc];

[~,I] = sortrows(sorting_criteria, [1 4]);
plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'occupied bandwidth', 'occupied_bandwidth');

[~,I] = sortrows(sorting_criteria, [1 3 8]);
plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'occupied bandwidth - midband - snr', 'occupied_bandwidth-midband-snr');
    
[~,I] = sortrows(sorting_criteria, [1 3 9]);
plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'occupied bandwidth - midband - cc', 'occupied_bandwidth-midband-cc');
    
[~,I] = sort(obs_struct.snr);
plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
    'signal-to-noise ratio', 'snr');

[~,I] = sortrows(sorting_criteria, 5);
plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'back azimuth', 'back_azimuth');

[~,I] = sortrows(sorting_criteria, 7);
plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'polarity', 'polarity');
    
if true
    plotter(OB_snr, OB_cc, fc_bin_mid, (1:size(obs_struct.fcorners,1))', ...
        'event id and station id', 'nosort');

    [~,I] = sort(obs_struct.snr);
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'signal-to-noise ratio', 'snr');

    [~,I] = sort(obs_struct.CCmaxs(:,2));
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'maximum correlation coefficient', 'cc');

    [~,I] = sort(obs_struct.metadata.MAG);
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'event magnitude', 'mag');

    [~,I] = sort(obs_struct.metadata.GCARC);
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'epicentral distance', 'gcarc');

    [~,I] = sort(obs_struct.metadata.EVDP);
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'event depth', 'evdp');

    [~,I] = sort(obs_struct.metadata.MAG / (obs_struct.metadata.GCARC .^ 2));
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'event magnitude / distance^2', 'mag_by_gcarc_squared');

    [~,I] = sort(obs_struct.metadata.USER9);
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'ray parameter', 'rayparam');

    [~,I] = sort(-obs_struct.metadata.STEL);
    plotter(OB_snr, OB_cc, fc_bin_mid, I, ...
        'seafloor depth (shallow to deep)', 'stel');
end

end

% make plots
function plotter(OB_snr, OB_cc, fc_bin_mid, sortindex, sortedby, savename)
    % create a plot
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 1 10 8])
    
    % SNR
    ax1 = subplot('Position', [0.08 0.54 0.88 0.41]);
    imagesc([1 size(OB_snr,1)], [fc_bin_mid(1) fc_bin_mid(end)], ...
        log10(OB_snr(sortindex,:)'));
    setimagenan(ax1, ax1.Children, [0 0 0]);
    ylabel('frequency (Hz)')
    xlabel('trace rank')
    c1 = colorbar('SouthOutside');
    c1.Limits = log10([3 10000]);
    c1.Ticks = log10([3 10 30 100 300 1000 3000 10000]);
    c1.TickLabels = [3 10 30 100 300 1000 3000 10000];
    c1.TickDirection = 'out';
    axis xy

    % CC
    ax2 = subplot('Position', [0.08 0.04 0.88 0.41]);
    imagesc([1 size(OB_cc,1)], [fc_bin_mid(1) fc_bin_mid(end)], ...
        OB_cc(sortindex,:)');
    setimagenan(ax2, ax2.Children, [0 0 0]);
    ylabel('frequency (Hz)')
    xlabel('trace rank')
    c2 =  colorbar('SouthOutside');
    c2.Limits = [0 1];
    c2.TickDirection = 'out';
    axis xy

    set(ax1.Parent.Children(1), 'FontSize', 11)
    set(ax1.Parent.Children(2), 'FontSize', 11)
    set(ax1.Parent.Children(3), 'FontSize', 11)
    set(ax1.Parent.Children(4), 'FontSize', 11)
    
    set(ax1.Parent.Children(2), 'TickDir', 'out', 'Box', 'on')
    set(ax1.Parent.Children(4), 'TickDir', 'out', 'Box', 'on')

    c1.Label.String = 'SNR';
    c2.Label.String = 'maximum correlation coefficient';

    % throwaway title axes
    axt = subplot('Position', [0 0.96 1 0.01]);
    title(sprintf('Sorted by %s', sortedby), 'FontSize', 12)
    [axt.Title.Position(1), axt.Title.Position(2)] = ...
            norm2trueposition(axt, 0.5, 1.06);
    axt.XAxis.Visible = 'off';
    axt.YAxis.Visible = 'off';
    axt.Color = 'none';

    % save the figure
%     set(gcf, 'Renderer', 'opengl')
%     filename = sprintf('%s_%s_opengl.eps', mfilename, savename);
%     figdisp(filename, [], '-r1200', 2, [], 'epstopdf');
    
    set(gcf, 'Renderer', 'painters')
    filename = sprintf('%s_%s.eps', mfilename, savename);
    figdisp(filename, [], '-r1200', 2, [], 'epstopdf');
end