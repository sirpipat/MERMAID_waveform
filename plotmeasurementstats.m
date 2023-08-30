function plotmeasurementstats(obs_struct)
% PLOTMEASUREMENTSTATS(obs_struct)
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
%   - presiduals        InstaSeis arrival - TauP prediction for rirst P
%                       arrival
%
% Last modified by sirawich-at-princeton.edu: 08/29/2023

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
    variableconstructor('t_rel_joel', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-80 80], 10, i_ttravel);
    variableconstructor('t_rel_correct', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-80 80], 10, i_ttravel);
    variableconstructor('t_rel_joel2', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-40 40], 1, and(i_ttravel, i_dlnt_joel2));
    variableconstructor('t_rel_correct2', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-40 40], 1, and(i_ttravel, i_dlnt_corrected2));
    variableconstructor('t_rel_joel3', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-10 10], 0.5, and(i_ttravel, i_dlnt_joel3));
    variableconstructor('t_rel_correct3', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-10 10], 0.5, and(i_ttravel, i_dlnt_corrected3));
    variableconstructor('t_rel_joel4', dlnt_joel * 100, 'relative travel time residual from Simon et al. 2022 (%)', [-3 3], 0.2, and(i_ttravel, i_dlnt_joel4));
    variableconstructor('t_rel_correct4', dlnt_corrected * 100, 'adjusted relative correlation travel time residual (%)', [-3 3], 0.2, and(i_ttravel, i_dlnt_corrected4));
    variableconstructor('gcarc', obs_struct.metadata.GCARC, 'great-circle epicentral distance (degree)', [0 180], 5, []);
    variableconstructor('log10gcarc', log10(obs_struct.metadata.GCARC), 'log_{10}great-circle epicentral distance (degree)', [], 0.1, []);
    variableconstructor('baz', obs_struct.metadata.BAZ, 'back azimuth (degree)', [0 360], 0:30:360, []);
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
    8 9 2;
    8 9 3;
    8 9 4;
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
];

%% make histograms of time shifts, maximum correlation
for ii = 100:length(variables)
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
    vline(gca, varmed, 'Color', 'k', ...
        'LineWidth', 2, 'LineStyle', '-.');
    
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
    i_var1 = var1.indices;
    i_var2 = var2.indices;
    i_var = and(i_var1, i_var2);
    if ~isempty(var3)
        i_var3 = var3.indices;
        i_var = and(i_var, i_var3);
    end
        
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
    if isempty(var3)
        scathistplot(var1.value(i_var), var2.value(i_var), [], [], ...
            var1.label, var2.label, [], var1.axlimit, var2.axlimit, [], ...
            histx_arg, histy_arg, {'SizeData', 9});
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
    
    if strcmp(var1.name, 't_res_joel') && strcmp(var2.name, 't_res_correct')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
    elseif strcmp(var1.name, 't_rel_joel') && strcmp(var2.name, 't_rel_correct')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
    elseif strcmp(var1.name, 't_rel_joel2') && strcmp(var2.name, 't_rel_correct2')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
    elseif strcmp(var1.name, 't_rel_joel3') && strcmp(var2.name, 't_rel_correct3')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
    elseif strcmp(var1.name, 't_rel_joel4') && strcmp(var2.name, 't_rel_correct4')
        rf = refline(ax_scat, 1, 0);
        set(rf, 'LineWidth', 1, 'Color', 'k')
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
% incdices      logical array whether to use each element in value or not
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