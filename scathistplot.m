function [fig, ax_histx, ax_histy, ax_scat, cb_scat] = ...
    scathistplot(x, y, c, savename, labelx, labely, labelc, ...
    limx, limy, limc, histx_arg, histy_arg, scat_arg, cb_arg)
% [fig, ax_histx, ax_histy, ax_scat] = ....
%     scathistplot(x, y, c, savename, labelx, labely, labelc, ...
%                  limx, limy, limc, histx_arg, histy_arg, scat_arg, ...
%                  cb_arg)
%
% Plots a scatter plot of (X,Y) values as well as histograms of X and Y
% with the X=Y reference line and the medians of X and Y. You may specify
% the keyword arguments for histograms and a scatter plot using the last
% three arguments HISTX_ARG, HISTY_ARG, and SCAT_ARG as cell arrays
% containing pairs of keywords and values.
%
%   +----------------------------+
%   | +--------------+           |
%   | |              |           |
%   | |   ax_histx   |           |
%   | +--------------+           |
%   | +--------------+ +-------+ |
%   | |              | |       | |
%   | |              | | ax_   | |
%   | |   ax_scat    | |       | |
%   | |              | | histy | |
%   | |              | |       | |
%   | +--------------+ +-------+ |
%   +----------------------------+
%
% INPUT:
% x             x-values [required]
% y             y-values [required]
% c             scatter plot face color [optional: a color name, 
%               a RGB triplet, or an array with the same size as x and y]
%               Note: If you want to plot 3 data points each with differnt
%               color values, consider appending nan to the end of x, y, 
%               and, c, so that the function does not confuse c values as 
%               an RGB triplet.
% savename      name of the file [optinal; if it is not specified, the
%               figure is not saved.]
% labelx        x-label of the scatter plot
% labely        y-label of the scatter plot
% labelc        colorbar label of the scatter plot (only works when c is an
%               array of the same size as x and y)
% limx          x-limit of the plots
% limy          y-limit of the plots
% limc          c-limit of the colorbar
% histx_arg     keyword arguments for the histogram of X
% histy_arg     keyword arguments for the histogram of Y
% scat_arg      keyword arguments for the scatter plot
% cb_arg        keyword arguments for the colorbar
%
% OUTPUT:
% fig           figure handle
% ax_histx      axes handle of the histogram of X
% ax_histy      axes handle of the histogram of Y
% ax_scat       axes handle of the scatter plot
% cb_scat       colorbar handle of the scatter plot (return empty if the
%               scatter plot color is empty or an RGB triplet)
%
% SEE ALSO:
% COLORBAR, HISTOGRAM, SCATTER
%
% EXAMPLE:
% x = 5 * randn([100,1]) - 0.05;
% y = 5 * randn([100,1]) + 0.05;
% c = randn([100,1]);
% scathistplot(x, y, c, 'demo', 'x', 'y', [-20 20], [-20 20], 
%     {'BinWidth', 0.25}, {'BinWidth', 0.25}, ...
%     {'SizeData', 16, 'MarkerFaceColor', 'r'});
%
% Last modified by sirawich-at-princeton.edu: 07/24/2023

defval('c', rgbcolor('1'))
defval('savename', false)
defval('labelx', 'x')
defval('labely', 'y')
defval('labelc', 'value')
defval('limx', nan)
defval('limy', nan)
defval('limc', nan)
defval('histx_arg', {'FaceAlpha', 1})
defval('histy_arg', {'FaceAlpha', 1})
defval('scat_arg', {'Marker', 'o'})
defval('cb_arg', {})

if ischar(c)
    cname = c;
    c = csscolor(c);
    fprintf('CONVERTING COLOR ''%s'' TO [%.4f %.4f %.4f]\n', cname, ...
        c(1), c(2), c(3));
end

x_med = median(x);
y_med = median(y);

% create figure
fig = figure(2);
clf
set(fig, 'Units', 'inches', 'Position', [0 1 10 10])

% figure layout
%tiledlayout(3, 3, 'TileSpacing', 'compact');

% histogram of x
ax_histx = subplot('Position', [0.1300 0.7000 0.5050 0.2350]);%0.2484nexttile([1 2]);%
histogram(ax_histx, x, histx_arg{:})
grid on
box on
hold on
if ~isnan(limx)
    xlim(limx)
end
ylabel('counts')
set(ax_histx, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', ...
    'XAxisLocation', 'top')
vline(ax_histx, x_med, 'Color', 'k', 'LineWidth', 2, ...
    'LineStyle', '-.');

% histogram of y
ax_histy = subplot('Position', [0.6700 0.1600 0.2350 0.5050]);%0.5317nexttile(6, [2 1]);%
histogram(ax_histy, y, 'Orientation', 'horizontal', histy_arg{:})
grid on
box on
if ~isnan(limy)
    ylim(limy)
end
xlabel('counts')
set(ax_histy, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on', ...
    'YAxisLocation', 'right')
hline(ax_histy, y_med, 'Color', 'k', 'LineWidth', 2, ...
    'LineStyle', '-.');

% scatter plot
ax_scat = subplot('Position', [0.1300 0.1600 0.5050 0.5050]);
scatter(ax_scat, x, y, [], c, 'filled', scat_arg{:})
grid on
box on
hold on
xlim(ax_histx.XLim)
ylim(ax_histy.YLim)
vline(ax_scat, x_med, 'Color', 'k', 'LineWidth', 2, ...
    'LineStyle', '-.');
hline(ax_scat, y_med, 'Color', 'k', 'LineWidth', 2, ...
    'LineStyle', '-.');
if limx == limy
    rl = refline(ax_scat, 1, 0);
    set(rl, 'Color', [0 0 0], 'LineWidth', 1)
end
xlabel(labelx)
ylabel(labely)
set(ax_scat, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')

% add colorbar
if all(size(c) == [1 3]) || all(size(c) == [3 1])
    cb_scat = [];
else
    ax_scat_position = ax_scat.Position;
    colormap('parula')
    cb_scat = colorbar(ax_scat, 'Location', 'southoutside', cb_arg{:});
    cb_scat.TickDirection = 'out';
    cb_scat.Label.String = labelc;
    if ~isnan(limc)
        cb_scat.Limits = limc;
    end
    ax_scat.Position = ax_scat_position;
end

% add plot labels
ax_histxb = boxedlabel(ax_histx, 'northwest', 0.027, 'norm', 'a', ...
    'FontSize', 14);
axes(ax_histxb)

ax_histyb = boxedlabel(ax_histy, 'northeast', 0.027, 'norm', 'b', ...
    'FontSize', 14);
axes(ax_histyb)

ax_scatb = boxedlabel(ax_scat, 'northwest', 0.027, 'norm', 'c', ...
    'FontSize', 14);
axes(ax_scatb)

% title
if all(size(c) == [1 3]) || all(size(c) == [3 1]) 
    title(ax_histx, sprintf('n = %d, x = %s (median = %.2f), y = %s (median = %.2f)', ...
        length(x), labelx, x_med, labely, y_med), 'FontSize', 12)
else
    title(ax_histx, sprintf('n = %d, x = %s (median = %.2f), y = %s (median = %.2f), c = %s', ...
        length(x), labelx, x_med, labely, y_med, labelc), 'FontSize', 12)
end
[ax_histx.Title.Position(1), ax_histx.Title.Position(2)] = ...
        norm2trueposition(ax_histx, (0.5 - ax_histx.Position(1)) / ...
        ax_histx.Position(3), 1.18);
    
set(fig, 'Renderer', 'painters')
if savename
    fname = sprintf('%s_%s.eps', mfilename, savename);
    figdisp(fname, [], [], 2, [], 'epstopdf');
end
end

