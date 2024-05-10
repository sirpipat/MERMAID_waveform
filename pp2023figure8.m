function varargout = pp2023figure8
% fig = PP2023FIGURE8
%
% Makes figure 8 of Pipatprathanporn+2024
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 05/10/2024
%% Great-circle path bathymetry near the receiver
ddir = sprintf('%sDATA/Figure8/bathymetry/', getenv('MERMAID2'));
[t_s, x_s, t_p, x_p, t_rf, rf, x_p_conv] = responsedemodata(ddir);

%% plot
figure(8)
set(gcf, 'Units', 'inches', 'Position', [8 1 8 4])
clf

ax1 = subplot('Position', [0.09 0.74 0.88 0.24]);
plot(ax1, t_p, x_p_conv / max(abs(x_p)), 'r', 'LineWidth', 1.4)
hold on
plot(ax1, t_p, x_p / max(abs(x_p)), 'b', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.5, 1.5])
nolabels(ax1, 1)
ylabel('press. (Pa)')

% label graph
text(ax1, 29.5, 0.9, 'pressure at hydrophone ($$p$$)', ...
    'FontSize', 14, 'Color', 'b', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right')
text(ax1, 29.5, -1.0, 'estimated pressure at hydrophone ($$\hat{p} = s*r$$)', ...
    'FontSize', 14, 'Color', 'r', 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right')

% label scale
% scalelabel = sprintf('%.2g}', max(abs(x_p)));
% scalelabel = replace(scalelabel, 'e', '\times10^{');
% scalelabel = replace(scalelabel, 'x', '\times');
% yticklabels({['-' scalelabel], 0 ,scalelabel})
yticklabels({'', '', ''})

set(ax1, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax1, 'northwest', 0.28, [], 'a', 'FontSize', 14);

ax2 = subplot('Position', [0.09 0.44 0.88 0.24]);
hold on
plot(ax2, t_s, x_s / max(abs(x_s)), 'Color', [0 0.4 0], 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.2, 1.2])
nolabels(ax2, 1)
ylabel('disp. (m)')

% label graph
text(ax2, 29.5, 0.4, 'displacement at ocean bottom ($$s$$)', ...
    'FontSize', 14, 'Color', [0 0.4 0], 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right')

% label scale
% scalelabel = sprintf('%.2g}', max(abs(x_s)));
% scalelabel = replace(scalelabel, 'e', '\times10^{');
% scalelabel = replace(scalelabel, 'x', '\times');
% yticklabels({['-' scalelabel], 0 ,scalelabel})
yticklabels({'', '', ''})

set(ax2, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax2, 'northwest', 0.28, [], 'b', 'FontSize', 14);

ax3 = subplot('Position', [0.09 0.14 0.88 0.24]);
hold on
plot(ax3, t_rf, rf / max(abs(rf)), 'k', 'LineWidth', 1)
grid on
xlim([0, 30])
ylim([-1.5, 1.5])
xlabel('time (s)')
ylabel('resp. (Pa/m)')

% label graph
text(ax3, 29.5, 1.05, ['response function ' ...
    '($$r =\mathcal{F}^{-1}(\tilde{p}/\tilde{s})$$)'], ...
    'FontSize', 14, 'Interpreter', 'latex', 'HorizontalAlignment', 'right')
yticklabels({'', '', ''})

set(ax3, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')

boxedlabel(ax3, 'northwest', 0.28, [], 'c', 'FontSize', 14);

%% Loop over false bathymetry
if true
    dirlist = {...
        fullfile(getenv('MERMAID2'), 'DATA', 'Figure8', 'bath_10936816_P0009_shift045deg/'), ...
        fullfile(getenv('MERMAID2'), 'DATA', 'Figure8', 'bath_10936816_P0009_shift315deg/') ...
        };
    COLOR_FOR_FALSES = [0.7 0.7 0.7; 0.7 0.7 0.7] .^ 1;

    for ii = 1:length(dirlist)
        [t_s, x_s, t_p, x_p, t_rf, rf, x_p_conv] = responsedemodata(dirlist{ii});
        
        plot(ax1, t_p, x_p / max(abs(x_p)), 'Color', COLOR_FOR_FALSES(ii,:), 'LineWidth', 0.5)
        ax1.Children = ax1.Children([2:end 1]);
        
        plot(ax2, t_s, x_s / max(abs(x_s)), 'Color', COLOR_FOR_FALSES(ii,:), 'LineWidth', 0.5)
        ax2.Children = ax2.Children([2:end 1]);

        plot(ax3, t_rf, rf / max(abs(rf)), 'Color', COLOR_FOR_FALSES(ii,:), 'LineWidth', 0.5)
        ax3.Children = ax3.Children([2:end 1]);
    end
end

set(gcf, 'Renderer', 'painters')

%% save figure
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');

if nargout > 0
    varargout{1} = fig;
end
end

%% [t_s, x_s, t_p, x_p, t_rf, rf, x_p_conv] = RESPONSEDEMODATA(ddir)
% 
% Generate time-series data for response function calculation demo for
% FLUIDSOLID SPECFEM2D simulations
%
% INPUT:
% ddir          FLUIDSOLID SPECFEM2D simulation directory
%
% OUTPUT:
% t_s           seismogram time at the ocean bottom receiver
% x_s           vertical displacement template waveform
%               at the ocean bottom receiver
% t_p           seismogram time at the hydrophone receiver
% x_p           acoustic pressure at the hydrophone receiver
% t_rf          time for the response function
% rf            response function
% x_p_conv      estimated acoustic pressure at the hydrophone receiver
function [t_s, x_s, t_p, x_p, t_rf, rf, x_p_conv] = ...
    responsedemodata(ddir)

%SPECFEM-vertical displacement (template)
[t_s, x_s] = getarrivaltemplate(ddir, 'bottom');

% SPECFEM-pressure
try
    file_p = [ddir 'OUTPUT_FILES/AA.S0001.PRE.semp'];
    [t_p, x_p] = read_seismogram(file_p);
catch ME
    file_p = [ddir 'OUTPUT_FILES/Up_file_single_p.bin'];
    par_file = [ddir 'DATA/Par_file'];
    seiss = freadseismograms(file_p, par_file);
    x_p = seiss(:,1);
    t_p = t_s;
end

% calculations
fs = 1 / (t_s(2) - t_s(1));

% response function
[~, ~, t_rf, rf, ~] = cctransplot(ddir, ddir, 'flat_10936816_P0009', ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, [], fs, false);

% convolve
x_p_conv = conv(x_s, rf);
x_p_conv = x_p_conv(1:length(x_p));
end