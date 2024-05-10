function varargout = pp2023figure9
% fig = PP2023FIGURE9
%
% Makes figure 9 of Pipatprathanporn+2024
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 05/10/2024

%% load data

% load INSTASEIS displacement seismogram at the ocean bottom
[seis_s, hdr_s] = readsac(sprintf(...
    '%sDATA/Figure9/20180817T153501_09_0_SYNTHETIC.sac', ...
    getenv('MERMAID2')));
[dt_ref_s, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);

% read seismograms from SPECFEM2D simulation
ddir = sprintf('%sDATA/Figure8/bathymetry/', getenv('MERMAID2'));
file_s = [ddir 'OUTPUT_FILES/AB.S0001.BXZ.semd'];
[t_s, x_s] = read_seismogram(file_s);

% read the observed pressure record
[seis_o, hdr_o] = readsac(sprintf(...
    '%sDATA/Figure9/20180817T154351.09_5B77394A.MER.DET.WLT5.sac', ...
    getenv('MERMAID2')));
[dt_ref_o, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);

%% calculations

% sampling rate of the seiemograms from SPECFEM2D simulation
fs = 1 / (t_s(2) - t_s(1));

% deconvolve for the response function due to the ocean layer
[~, ~, t_r, seis_r] = cctransplot(ddir, ddir, 'flat_10936816_P0009', ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, [], fs, false);

% time relative to the first P-wave arrival
t_relative = seconds(dts_o - dt_ref_o) - hdr_o.T0;


% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);

% determine the corner frequency
%[fcorners, ~] = freqselect(t_relative, pres_o, fs_o, false);

% filter
fcorners = [0.40 2.00];
pres_o = bandpass(pres_o, fs_o, fcorners(1), fcorners(2), 4, 2, 'butter', 'linear');

% timeshift
timeshift = 4.8501;

% resample to MERMAID datetimes
seis_s_interp = shannon(dts_s, seis_s, dts_o);
t_r_interp = 0:(1/fs_o):t_r(end);
seis_r = shannon(t_r, seis_r, t_r_interp);

% convolve for the synthetic pressure
pres_s = conv(seis_s_interp, seis_r);
pres_s = pres_s(1:length(seis_o), 1);
pres_s = bandpass(pres_s, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

% filter
seis_s_interp = bandpass(seis_s_interp, fs_o, fcorners(1), fcorners(2), ...
    4, 2, 'butter', 'linear');

%% plot
figure(9)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 4])
clf

ax1 = subplot('Position', [0.09 0.74 0.88 0.24]);
plot(ax1, t_relative, pres_o / max(abs(pres_o)), 'LineWidth', 1, 'Color', 'b')
hold on
plot(ax1, t_relative + 0 * timeshift, pres_s / max(abs(pres_s)), 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5])
plot(ax1, t_relative + timeshift, pres_s / max(abs(pres_s)), 'LineWidth', 1, 'Color', 'r')

grid on
xlim([-10 20])
ylim([-1.5, 1.5])
nolabels(ax1, 1)
ylabel('press. (Pa)')

% label graph
text(ax1, 19, -0.9, 'synthetic pressure ($$\hat{p} = s*r$$)', ...
    'FontSize', 14, 'HorizontalAlignment', 'right', 'Color', 'r', ...
    'Interpreter', 'latex')
text(ax1, 19, 1.0, 'observed pressure ($$p$$)', ...
    'FontSize', 14, 'HorizontalAlignment', 'right', 'Color', 'b', ...
    'Interpreter', 'latex')

% label scale
% scalelabel = sprintf('%.2g', max(abs(pres_o)));
% if contains(scalelabel, 'e')
%     if str2double(scalelabel(end-1:end)) > 0
%         scalelabel = replace(scalelabel, '+', '');
%     end
%     if abs(str2double(scalelabel(end-1:end))) < 10
%         scalelabel(end-1) = '';
%     end
%     scalelabel = [scalelabel '}'];
%     scalelabel = replace(scalelabel, 'e', '\times10^{');
%     scalelabel = replace(scalelabel, 'x', '\times');
% end
% yticklabels({['-' scalelabel], 0 ,scalelabel})
yticklabels({'', '', ''})

set(ax1, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')
boxedlabel(ax1, 'northwest', 0.28, [], 'a', 'FontSize', 14);

ax2 = subplot('Position', [0.09 0.44 0.88 0.24]);
plot(ax2, t_relative, seis_s_interp / max(abs(seis_s_interp)), 'LineWidth', 1, 'Color', [0 0.4 0])
grid on
xlim([-10 20])
ylim([-1.5, 1.5])
nolabels(ax2, 1)
ylabel('disp. (m)')

% label graph
text(ax2, 19, 1.0, 'synthetic displacement ($$s$$)', 'FontSize', 14, ...
    'HorizontalAlignment', 'right', 'Color', [0 0.4 0], ...
    'Interpreter', 'latex')

% label scale
% scalelabel = sprintf('%.2g', max(abs(seis_s_interp)));
% if contains(scalelabel, 'e')
%     if str2double(scalelabel(end-1:end)) > 0
%         scalelabel = replace(scalelabel, '+', '');
%     end
%     if abs(str2double(scalelabel(end-1:end))) < 10
%         scalelabel(end-1) = '';
%     end
%     scalelabel = [scalelabel '}'];
%     scalelabel = replace(scalelabel, 'e', '\times10^{');
%     scalelabel = replace(scalelabel, 'x', '\times');
% end
% yticklabels({['-' scalelabel], 0 ,scalelabel})
yticklabels({'', '', ''})

set(ax2, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')
boxedlabel(ax2, 'northwest', 0.28, [], 'b', 'FontSize', 14);

ax3 = subplot('Position', [0.09 0.14 0.88 0.24]);
hold on
plot(ax3, [-10 -5 -1 t_r_interp 33 35], ...
    [0; 0; 0; seis_r / max(abs(seis_r)); 0; 0], 'k', 'LineWidth', 1)
grid on
xlim([0 30])
ylim([-1.5, 1.5])
xlabel('time (s)')
ylabel('resp. (Pa/m)')

% label graph
text(ax3, 29, 1.05, 'response function ($$r$$)', 'FontSize', 14, ...
    'HorizontalAlignment', 'right', 'Interpreter', 'latex')

% label scale
% scalelabel = sprintf('%.2g', max(abs(seis_r)));
% if contains(scalelabel, 'e')
%     if str2double(scalelabel(end-1:end)) > 0
%         scalelabel = replace(scalelabel, '+', '');
%     end
%     if abs(str2double(scalelabel(end-1:end))) < 10
%         scalelabel(end-1) = '';
%     end
%     scalelabel = [scalelabel '}'];
%     scalelabel = replace(scalelabel, 'e', '\times10^{');
%     scalelabel = replace(scalelabel, 'x', '\times');
% end
% yticklabels({['-' scalelabel], 0 ,scalelabel})
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
    Y_SHIFT = [0 0] * 1.0;

    for ii = 1:length(dirlist)
        [~, ~, t_r, seis_r] = cctransplot(dirlist{ii}, dirlist{ii}, ...
            'bath_10936816_P0009', {'bottom', 'displacement'}, ...
            {'hydrophone', 'pressure'}, [], fs, false);
        plot(ax3, [-10; -5; -1; t_r; 33; 35], ...
            [0; 0; 0; seis_r / max(abs(seis_r)); 0; 0] + Y_SHIFT(ii), ...
            'Color', COLOR_FOR_FALSES(ii,:), 'LineWidth', 0.5)
        ax3.Children = ax3.Children([2:end 1]);
        
%         % resample to MERMAID datetimes
%         t_r_interp = (0:(1/fs_o):t_r(end))';
%         seis_r = shannon(t_r, seis_r, t_r_interp);
% 
%         % convolve for the synthetic pressure
%         pres_s = conv(seis_s_interp, seis_r);
%         pres_s = pres_s(1:length(seis_o), 1);
%         pres_s = bandpass(pres_s, fs_o, fcorners(1), fcorners(2), 4, 2, ...
%             'butter', 'linear');
%         
%         plot(ax1, t_relative + timeshift, ...
%             pres_s / max(abs(pres_s)) + Y_SHIFT(ii), ...
%             'LineWidth', 0.5, 'Color', COLOR_FOR_FALSES(ii,:))
%         ax1.Children = ax1.Children([2:end 1]);
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
