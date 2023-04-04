function varargout = pp2023figure7
% fig = PP2023FIGURE7
%
% Makes figure 7 of Pipatprathanporn+2023
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 04/04/2023
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

ax1 = subplot('Position', [0.08 0.54 0.88 0.44]);
plot(ax1, t_s, x_s / max(abs(x_s)), 'k', 'LineWidth', 1)
hold on
plot(ax1, t_p, x_p_conv / max(abs(x_p)) - 2, 'r', 'LineWidth', 1.3)
plot(ax1, t_p, x_p / max(abs(x_p)) - 2, 'k', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-3.6, 1.2])
xlabel('time (s)')
yticks([])

% label graph
text(ax1, 20, 0.4, 's: displacement at ocean bottom', 'FontSize', 12)
text(ax1, 20, -1.1, 'p: pressure at hydrophone', 'FontSize', 12)
text(ax1, 15, -2.9, 'p* = s*r: estimated pressure at hydrophone', 'FontSize', 12, 'Color', 'r')

% label scale
scalelabel = sprintf('( x %.2g} )', max(abs(x_s)));
scalelabel = replace(scalelabel, 'e', '\times 10^{');
text(ax1, 3, 0.4, scalelabel, 'FontSize', 12)
scalelabel = sprintf('( x %.2g} )', max(abs(x_p)));
scalelabel = replace(scalelabel, 'e', '\times 10^{');
text(ax1, 3, -3.2, scalelabel, 'FontSize', 12)

set(ax1, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax1, 'northwest', 0.28, [], 'a', 'FontSize', 14);

ax2 = subplot('Position', [0.08 0.14 0.88 0.22]);
plot(ax2, t_rf, rf / max(abs(rf)), 'k', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.2, 1.2])
xlabel('time (s)')
yticks([])

% label graph
text(ax2, 20, 0.7, 'r = s/p: response function', 'FontSize', 12)

% label scale
scalelabel = sprintf('( x %.2g} )', max(abs(rf)));
scalelabel = replace(scalelabel, 'e', '\times 10^{');
scalelabel = replace(scalelabel, '+0', '');
text(ax2, 3, 0.8, scalelabel, 'FontSize', 12)

set(ax2, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax2, 'northwest', 0.28, [], 'b', 'FontSize', 14);

set(gcf, 'Renderer', 'painters')

%% save figure
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');

if nargout > 0
    varargout{1} = fig;
end
end