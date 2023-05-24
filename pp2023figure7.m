function varargout = pp2023figure7
% fig = PP2023FIGURE7
%
% Makes figure 7 of Pipatprathanporn+2023
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 05/04/2023
%% load data

% SPECFEM-vertical displacement (template)
ddir = sprintf('%sDATA/Figure7/bathymetry/', getenv('MERMAID2'));
[t_s, x_s] = getarrivaltemplate(ddir, 'flat_10996154_P0009', 'bottom');

% SPECFEM-pressure
file_p = [ddir 'OUTPUT_FILES/AA.S0001.PRE.semp'];
[t_p, x_p] = read_seismogram(file_p);

%% calculations
fs = 1 / (t_s(2) - t_s(1));

% response function
[t_cc, cc, t_rf, rf, d] = cctransplot(ddir, ddir, 'flat_10996154_P0009', ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, [], fs, false);

% convolve
x_p_conv = conv(x_s, rf);
x_p_conv = x_p_conv(1:length(x_p));

%% plot
figure(7)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 4])
clf

ax1 = subplot('Position', [0.09 0.74 0.88 0.24]);
plot(ax1, t_p, x_p_conv / max(abs(x_p)), 'r', 'LineWidth', 1.3)
hold on
plot(ax1, t_p, x_p / max(abs(x_p)), 'k', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.5, 1.5])
nolabels(ax1, 1)

% label graph
text(ax1, 19.5, 0.9, 'pressure at hydrophone ($$p$$)', ...
    'FontSize', 14, 'Interpreter', 'latex')
text(ax1, 13, -1.0, 'estimated pressure at hydrophone ($$\hat{p} = s*r$$)', ...
    'FontSize', 14, 'Color', 'r', 'Interpreter', 'latex')

% label scale
scalelabel = sprintf('%.2g}', max(abs(x_p)));
scalelabel = replace(scalelabel, 'e', '\times10^{');
scalelabel = replace(scalelabel, 'x', '\times');
yticklabels({['-' scalelabel], 0 ,scalelabel})

set(ax1, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax1, 'northwest', 0.28, [], 'a', 'FontSize', 14);

ax2 = subplot('Position', [0.09 0.44 0.88 0.24]);
plot(ax2, t_s, x_s / max(abs(x_s)), 'k', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.2, 1.2])
nolabels(ax2, 1)

% label graph
text(ax2, 17, 0.4, 'displacement at ocean bottom ($$s$$)', ...
    'FontSize', 14, 'Interpreter', 'latex')

% label scale
scalelabel = sprintf('%.2g}', max(abs(x_s)));
scalelabel = replace(scalelabel, 'e', '\times10^{');
scalelabel = replace(scalelabel, 'x', '\times');
yticklabels({['-' scalelabel], 0 ,scalelabel})

set(ax2, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax2, 'northwest', 0.28, [], 'b', 'FontSize', 14);

ax3 = subplot('Position', [0.09 0.14 0.88 0.24]);
plot(ax3, t_rf, rf / max(abs(rf)), 'k', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.5, 1.5])
xlabel('time (s)')

% label graph
text(ax3, 19.3, 0.7, 'response function ($$r = s/p$$)', ...
    'FontSize', 14, 'Interpreter', 'latex')

% label scale
scalelabel = sprintf('%.2g}', max(abs(rf)));
scalelabel = replace(scalelabel, 'e', '\times10^{');
scalelabel = replace(scalelabel, '+0', '');
scalelabel = replace(scalelabel, 'x', '\times');
yticklabels({['-' scalelabel], 0 ,scalelabel})

set(ax3, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax3, 'northwest', 0.28, [], 'c', 'FontSize', 14);

set(gcf, 'Renderer', 'painters')

%% save figure
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');

if nargout > 0
    varargout{1} = fig;
end
end