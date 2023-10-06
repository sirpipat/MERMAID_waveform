function plotmeasurementstats(obs_struct, min_cc, min_snr, min_gcarc)
% PLOTMEASUREMENTSTATS(obs_struct, min_cc, min_snr, min_gcarc)
%
% Plots various travel time measurement statistics and the metadata for
% further analysis.
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
% min_cc            Correlation coefficient cut-off [default: -1]
% min_snr           Signal-to-noise ratio cut-off   [default: 0]
% min_gcarc         Epicentral distance cut-off     [default: 0]
%
% Last modified by sirawich-at-princeton.edu: 10/06/2023

defval('min_cc', -1)
defval('min_snr', 0)
defval('min_gcarc', 0)

%% calculate the derived variables
% relative travel time from correlation travel time
dlnt = obs_struct.t_shifts(:,2) ./ ...
    (obs_struct.metadata.T0 - obs_struct.metadata.USER8);

% travel time
ttravel = obs_struct.metadata.USER5;

% corrected time shift
t_shift_corrected = obs_struct.t_shifts(:,2) + obs_struct.presiduals;

% relative travel time from Joel's pick and corrected correlation travel
% time
dlnt_joel = obs_struct.metadata.USER4 ./ ttravel;
dlnt_corrected = t_shift_corrected ./ ttravel;

%% filter out for some variables
% remove outliers for SNR
[~, i_snr] = rmoutliers(obs_struct.snr);
i_snr = ~i_snr;

% remove outliers for relative travel time
i_dlnt2 = and(dlnt >= -0.03, dlnt <= 0.03);

% remove negative travel time
i_ttravel = (ttravel > 0);

% limit the presiduals to [-5 5] seconds
i_presiduals = and(obs_struct.presiduals >= -4, obs_struct.presiduals <= 4);

% limit the relative travel time to [-40 40] percents
i_dlnt_joel2 = and(dlnt_joel >= -0.4, dlnt_joel <= 0.4);
i_dlnt_corrected2 = and(dlnt_corrected >= -0.4, dlnt_corrected <= 0.4);

% limit the relative travel time to [-10 10] percents
i_dlnt_joel3 = and(dlnt_joel >= -0.1, dlnt_joel <= 0.1);
i_dlnt_corrected3 = and(dlnt_corrected >= -0.1, dlnt_corrected <= 0.1);

% limit the relative travel time to [-3 3] percents
i_dlnt_joel4 = and(dlnt_joel >= -0.03, dlnt_joel <= 0.03);
i_dlnt_corrected4 = and(dlnt_corrected >= -0.03, dlnt_corrected <= 0.03);

% ANOTHER WAY TO LIMIT: using percentile maybe from 5th to 95th

% limit to only good match
i_cc = (obs_struct.CCmaxs(:,2) >= min_cc);
i_snr2 = (obs_struct.snr >= min_snr);
i_gcarc = (obs_struct.metadata.GCARC >= min_gcarc);
i_dlnt = and(dlnt >= prctile(dlnt, 5), dlnt <= prctile(dlnt, 95));
i_mask = and(and(i_cc, i_snr2), i_gcarc);

%% list of things to plot
variables = [...
    variableconstructor('t_shift', obs_struct.t_shifts(:,2), 'time shift (s)', [], 1, []);
    variableconstructor('cc', obs_struct.CCmaxs(:,2), 'correlation coefficient', [0 1], 0.05, []);
    variableconstructor('log10snr', log10(obs_struct.snr), 'log_{10} signal-to-noise ratio', [], 0.1, []);
    variableconstructor('snr', obs_struct.snr, 'signal-to-noise ratio', [], 10, i_snr);
    variableconstructor('dlnt', dlnt * 100, 'relative time shift (%)', [], 10, []);
    variableconstructor('dlnt2', dlnt * 100, 'relative time shift, outliers removed (%)', [-3 3], 0.2, i_dlnt2);
    variableconstructor('ttravel', ttravel, 'travel time (s)', [], 100, i_ttravel);
    variableconstructor('t_res_joel', obs_struct.metadata.USER4, 'travel time residual from Simon et al. 2022 (s)', [-20 20], 1, i_ttravel);
    variableconstructor('t_res_correct', t_shift_corrected, 'adjusted correlation travel time residual (s)', [-20 20], 1, i_ttravel);
    variableconstructor('presidual', obs_struct.presiduals, 'InstaSeis - ray theory prediction of P-wave arrival on AK135 model (s)', [], 1, []);
    variableconstructor('presidual_limit', obs_struct.presiduals, 'InstaSeis - ray theory prediction of P-wave arrival on AK135 model (s)', [-4 4], 0.2, i_presiduals);
    variableconstructor('fc_lower', obs_struct.fcorners(:,1), 'lower corner frequency (Hz)', [0.375 1.525], 0.375:0.05:1.525, []);
    variableconstructor('fc_upper', obs_struct.fcorners(:,2), 'upper corner frequency (Hz)', [0.875 2.025], 0.875:0.05:2.025, []);
    variableconstructor('fc_mid', (obs_struct.fcorners(:,1) + obs_struct.fcorners(:,2)) / 2, 'median of frequency band (Hz)', [0.6375 1.7625], 0.6325:0.025:1.7625, []);
    variableconstructor('bandwidth', obs_struct.fcorners(:,2) - obs_struct.fcorners(:,1), 'bandwidth (Hz)', [0.475 1.625], 0.475:0.05:1.625, []);
    variableconstructor('t_rel_joel', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-40 40], 2, i_ttravel);
    variableconstructor('t_rel_correct', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-40 40], 2, i_ttravel);
    variableconstructor('t_rel_joel2', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-20 20], 1, and(i_ttravel, i_dlnt_joel2));
    variableconstructor('t_rel_correct2', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-20 20], 1, and(i_ttravel, i_dlnt_corrected2));
    variableconstructor('t_rel_joel3', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-10 10], 0.5, and(i_ttravel, i_dlnt_joel3));
    variableconstructor('t_rel_correct3', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-10 10], 0.5, and(i_ttravel, i_dlnt_corrected3));
    variableconstructor('t_rel_joel4', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-3 3], 0.2, and(i_ttravel, i_dlnt_joel4));
    variableconstructor('t_rel_correct4', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-3 3], 0.2, and(i_ttravel, i_dlnt_corrected4));
    variableconstructor('gcarc', obs_struct.metadata.GCARC, 'great-circle epicentral distance (degree)', [0 180], 5, []);
    variableconstructor('log10gcarc', log10(obs_struct.metadata.GCARC), 'log_{10}great-circle epicentral distance (degree)', [], 0.1, []);
    variableconstructor('baz', obs_struct.metadata.BAZ, 'back azimuth (degree)', [0 360], 0:30:360, []);
    variableconstructor('evdp', obs_struct.metadata.EVDP', 'event depth (km)', [0 700], 25, []);
];

variable_pairs = [...
    1 2 nan;
    3 2 nan;
    4 2 nan;
    1 5 nan;
    3 5 nan;
    4 5 nan;
    1 6 nan;
    3 6 nan;
    4 6 nan;
    7 5 nan;
    7 6 nan;
    9 8 2;
    9 8 3;
    9 8 4;
    9 8 24;
    10 12 nan;
    10 13 nan;
    10 14 nan;
    10 15 nan;
    11 12 nan;
    11 13 nan;
    11 14 nan;
    11 15 nan;
    16 17 2;
    18 19 2;
    20 21 2;
    22 23 2;
    16 17 3;
    18 19 3;
    20 21 3;
    22 23 3;
    16 17 4;
    18 19 4;
    20 21 4;
    22 23 4;
    16 17 24;
    18 19 24;
    20 21 24;
    22 23 24;
    16 17 25;
    18 19 25;
    20 21 25;
    22 23 25;
    24 16 2;
    24 17 2;
    24 7 1;
    24 7 2;
    24 7 26;
    1 24 2;
    9 24 2;
    1 7 2;
    9 7 2;
];

%% make histograms of time shifts, maximum correlation
for ii = 1:length(variables)
    figure(1)
    set(gcf, 'Units', 'inches', 'Position', [0 1 8 5])
    clf
    if length(variables(ii).BinWidth) == 1
        histogram(variables(ii).value(variables(ii).indices), ...
            'BinWidth', variables(ii).BinWidth);
    else
        histogram(variables(ii).value(variables(ii).indices), ...
            'BinEdges', variables(ii).BinWidth);
    end
    ax = gca;
    grid on
    if ~isempty(variables(ii).axlimit)
        xlim(variables(ii).axlimit)
    end
    xlabel(variables(ii).label)
    ylabel('counts')
    set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
    varmed = median(variables(ii).value(variables(ii).indices));
    [~, vl] = vline(gca, varmed, 'Color', 'k', ...
        'LineWidth', 2, 'LineStyle', '-.');
    uistack(vl, 'bottom')
    
    title(gca, sprintf('n = %d, x = %s (median = %.2f)', ...
        sum(variables(ii).indices), variables(ii).label, ...
        varmed), 'FontSize', 12)
    [ax.Title.Position(1), ax.Title.Position(2)] = ...
        norm2trueposition(ax, 0.5, 1.06);
    
    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%s_histogram.eps', mfilename, ...
        variables(ii).name);
    figdisp(savename, [], [], 2, [], 'epstopdf');
end

%% make scatter plots between "all" continuous, quantitative variables
for ii = 1:size(variable_pairs, 1)
    var1 = variables(variable_pairs(ii, 1));
    var2 = variables(variable_pairs(ii, 2));
    if ~isnan(variable_pairs(ii, 3))
        var3 = variables(variable_pairs(ii, 3));
    else
        var3 = [];
    end
    
    % check if filtering is needed
    i_var1 = and(var1.indices, i_mask);
    i_var2 = and(var2.indices, i_mask);
    i_var = and(i_var1, i_var2);
    if ~isempty(var3)
        i_var3 = and(var3.indices, i_mask);
        i_var = and(i_var, i_var3);
    end
    % change i_var from logical array to indices array
    % e.g. from [1 0 1 1 0 0 1] to [1 3 4 7]
    i_all = (1:length(var1.value))';
    i_var = i_all(i_var);
        
    % set the savename to the variables being plotted
    if isempty(var3)
        savename = sprintf('%s-v-%s', var1.name, var2.name);
    else
        savename = sprintf('%s-v-%s-w-%s', var1.name, var2.name, var3.name);
    end
    
    if length(var1.BinWidth) == 1
        histx_arg = {'BinWidth', var1.BinWidth};
    else
        histx_arg = {'BinEdges', var1.BinWidth};
    end
    if length(var2.BinWidth) == 1
        histy_arg = {'BinWidth', var2.BinWidth};
    else
        histy_arg = {'BinEdges', var2.BinWidth};
    end
    
    % reorder the dots accordingly
    if ~isempty(var3)
        switch var3.name
            case {'cc', 'log10snr', 'snr', 'gcarc', 'log10gcarc'}
                [~, i_sort] = sort(var3.value(i_var), 'ascend');
                i_var = i_var(i_sort);
            otherwise
        end
    end
    
    if isempty(var3)
        [~, ~, ~, ax_scat] = scathistplot(var1.value(i_var), ...
            var2.value(i_var), [], [], var1.label, var2.label, [], ...
            var1.axlimit, var2.axlimit, [], histx_arg, histy_arg, ...
            {'SizeData', 9});
    else
        [~, ~, ~, ax_scat] = scathistplot(var1.value(i_var), ...
            var2.value(i_var), var3.value(i_var), [], var1.label, ...
            var2.label, var3.label, var1.axlimit, var2.axlimit, ...
            var3.axlimit, histx_arg, histy_arg, {'SizeData', 9});
        if strcmp(var3.name, 'baz')
            colormap(ax_scat, 'hsv');
        else
            colormap(ax_scat, 'parula');
        end
    end
    
    % store i_var in scatter plot UserData for GETDATAPOINT.m
    for jj = 1:length(ax_scat.Children)
        if isa(ax_scat.Children(jj), ...
                'matlab.graphics.chart.primitive.Scatter')
            set(ax_scat.Children(jj), 'UserData', i_var);
            break
        end
    end
    
    % ADD lines paralleled to the refline (1 seconds apart)
    % COMPUTE the percentange of points fall within X seconds away from the
    % refline
    % MAYBE compute correlation between x and y
    rfweight = 0.5;
    rfcolor = [0.75 0.75 0.75];
    if strcmp(var2.name, 't_res_joel') && strcmp(var1.name, 't_res_correct')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        rf = refline(ax_scat, 1, 2);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        rf = refline(ax_scat, 1, 4);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        rf = refline(ax_scat, 1, 6);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        rf = refline(ax_scat, 1, -2);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        rf = refline(ax_scat, 1, -4);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        rf = refline(ax_scat, 1, -6);
        set(rf, 'LineWidth', rfweight, 'Color', rfcolor)
        uistack(rf, 'bottom')
        set(ax_scat, 'XLim', var1.axlimit, 'YLim', var2.axlimit)
        
        % add # of observations within the regions
        diff = var2.value(i_var) - var1.value(i_var);
        p2 = sum(and(diff > 0, diff <= 2));
        p4 = sum(and(diff > 2, diff <= 4));
        p6 = sum(and(diff > 4, diff <= 6));
        pp = sum(diff > 6);
        n2 = sum(and(diff > -2, diff <= 0));
        n4 = sum(and(diff > -4, diff <= -2));
        n6 = sum(and(diff > -6, diff <= -4));
        nn = sum(diff <= -6);
        
        text(ax_scat, 10, -1, sprintf('pp = %.2f %%', pp * 100 / length(diff)))
        text(ax_scat, 10, -3, sprintf('p6 = %.2f %%', p6 * 100 / length(diff)))
        text(ax_scat, 10, -5, sprintf('p4 = %.2f %%', p4 * 100 / length(diff)))
        text(ax_scat, 10, -7, sprintf('p2 = %.2f %%', p2 * 100 / length(diff)))
        text(ax_scat, 10, -9, sprintf('n2 = %.2f %%', n2 * 100 / length(diff)))
        text(ax_scat, 10, -11, sprintf('n4 = %.2f %%', n4 * 100 / length(diff)))
        text(ax_scat, 10, -13, sprintf('n6 = %.2f %%', n6 * 100 / length(diff)))
        text(ax_scat, 10, -15, sprintf('nn = %.2f %%', nn * 100 / length(diff)))
        
        % Stats
        stats = regstats(var2.value(i_var), var1.value(i_var), ...
            'linear', {'beta', 'rsquare', 'fstat', 'tstat'});
        
        text(ax_scat, 10, -17, sprintf('R^2 = %.2f', stats.rsquare))
        text(ax_scat, 10, -19, sprintf('pval = %.2g', stats.fstat.pval))
    elseif strcmp(var1.name, 't_rel_joel') && strcmp(var2.name, 't_rel_correct')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
        uistack(rf, 'bottom')
    elseif strcmp(var1.name, 't_rel_joel2') && strcmp(var2.name, 't_rel_correct2')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
        uistack(rf, 'bottom')
    elseif strcmp(var1.name, 't_rel_joel3') && strcmp(var2.name, 't_rel_correct3')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
        uistack(rf, 'bottom')
    elseif strcmp(var1.name, 't_rel_joel4') && strcmp(var2.name, 't_rel_correct4')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
        uistack(rf, 'bottom')
    end
    
    fname = sprintf('%s_%s.eps', 'scathistplot', savename);
    figdisp(fname, [], [], 2, [], 'epstopdf');
    
    % save the figure for later call
    fname = fullfile(getenv('IFILES'), 'FIGURES', ...
        sprintf('%s_%s.fig', mfilename, savename));
    savefig(fname);
end
end

% Constructs a struct storing variables for plotting/saving
% histograms/scatterplots
%
% INPUT:
% name          variable name
% value         value of the variable
% label         axes label when plotting this variable
% axlimit       axes limit when plotting this variable
% BinWidth      histogram bin width (if given as a scalar) or
%               histogram bin edges (if given as a vector)
% indices       logical array whether to use each element in value or not
%               (it has to be the same size as value)
%
% OUTPUT:
% variable      a struct
function variable = variableconstructor(name, value, label, axlimit, BinWidth, indices)
defval('name', 'unnamed')
defval('label', 'unnamed')
defval('indices', true(size(value)))

variable = struct('name', name, 'value', value, 'label', label, ...
    'axlimit', axlimit, 'BinWidth', BinWidth, 'indices', indices);
end