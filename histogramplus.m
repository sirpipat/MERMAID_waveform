function histogramplus(ax, x, savename, varargin)
% HISTOGRAMPLUS(ax, x, savename, varargin)
%
% Plots a histogram of X on an axes AX. Then, added mean(X), median(X) and
% one standard deviation from the median. The values are included in the
% axes title.
%
% INPUT:
% ax            target axes
% x             values to plot
% savename      figure file name saved to $EPS
% varargin      keyword arguments for HISTOGRAM
%
% OUTPUT:
% no output besides a figure saved at $EPS. Make sure that the enviornment
% variable $EPS is set to a directory.
%
% Last modified by sirawich-at-princeton.edu, 01/21/2022

axes(ax)
cla

% plot the histogram
histogram(x, varargin{:})
hold on
% add the mean line
[~,v1] = vline(ax, mean(x), 'Color', 'k', 'LineWidth', 1.5, ...
    'LineStyle', '-.');
% add the median line
[~,v2] = vline(ax, median(x), 'Color', 'r', 'LineWidth', 1.5, ...
    'LineStyle', '--');
% add 1-std from the median
[~,v3] = vline(ax, median(x) + std(x) * [-1 1], 'Color', [0.1 0.4 0.9], ...
    'LineWidth', 1.5, 'LineStyle', '--');
grid on
set(gca, 'FontSize', 12, 'TickDir', 'both');
xlabel('residual (s)')
ylabel('counts')
title(sprintf('n = %d, mean = %.2f, median = %.2f, std = %.2f', ...
    length(x), mean(x), median(x), std(x)));
legend([v1 v2 v3(1)], 'mean', 'median', '1 std from median', ...
    'Location', 'best')

% save the figure
set(gcf, 'Renderer', 'painters')
figdisp(savename,[],[],2,[],'epstopdf');
end