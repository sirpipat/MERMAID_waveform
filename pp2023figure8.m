function varargout = pp2023figure8
% fig = PP2023FIGURE8
%
% Makes figure 8 of Pipatprathanporn+2023
%
% OUTPUT:
% fig       figure handle to the plots
%
% Last modified by sirawich-at-princeton.edu: 04/04/2023

%% load data

% load INSTASEIS displacement seismogram at the ocean bottom
[seis_s, hdr_s] = readsac(sprintf(...
    '%sDATA/Figure8/20190115T180634_09_0_SYNTHETIC.sac', ...
    getenv('MERMAID2')));
[dt_ref_s, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);

% read seismograms from SPECFEM2D simulation
ddir = sprintf('%sDATA/Figure7/bathymetry/', getenv('MERMAID2'));
file_s = [ddir 'OUTPUT_FILES/AB.S0001.BXZ.semd'];
[t_s, x_s] = read_seismogram(file_s);

% read the observed pressure record
[seis_o, hdr_o] = readsac(sprintf(...
    '%sDATA/Figure8/20190115T181021.09_5C3E73D6.MER.DET.WLT5.sac', ...
    getenv('MERMAID2')));
[dt_ref_o, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);

%% calculations

% sampling rate of the seiemograms from SPECFEM2D simulation
fs = 1 / (t_s(2) - t_s(1));

% deconvolve for the response function due to the ocean layer
[~, ~, t_r, seis_r] = cctransplot(ddir, ddir, 'flat_10996154_P0009', ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, [], fs, false);

% time relative to the first P-wave arrival
t_relative = seconds(dts_o - dt_ref_o) - hdr_o.T0;


% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);

% determine the corner frequency
%[fcorners, ~] = freqselect(t_relative, pres_o, fs_o, false);

% filter
fcorners = [0.85 1.45];
pres_o = bandpass(pres_o, fs_o, fcorners(1), fcorners(2), 4, 2, 'butter', 'linear');

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
figure(8)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 3.5])
clf

ax1 = subplot('Position', [0.04 0.20 0.92 0.76]);
plot(ax1, t_relative, seis_s_interp / max(abs(seis_s_interp)) + 4, 'LineWidth', 1, 'Color', 'k')
hold on
plot(ax1, t_relative, pres_s / max(abs(pres_s)) + 2, 'LineWidth', 1, 'Color', 'k')
plot(ax1, t_relative, pres_o / max(abs(pres_o)), 'LineWidth', 1, 'Color', 'k')
grid on

% label traces
text(ax1, -9, 4.6, 'synthetic displacement (s)', 'FontSize', 14)
text(ax1, -9, 2.6, 'synthetic pressure (p* = s*r)', 'FontSize', 14)
text(ax1, -9, 0.8, 'observed pressure', 'FontSize', 14)

% label scale
scalelabel = sprintf('( x %.2g} m )', max(abs(seis_s_interp)));
scalelabel = replace(scalelabel, 'e', '\times 10^{');
scalelabel = replace(scalelabel, '-0', '-');
text(ax1, 5, 4.6, scalelabel, 'FontSize', 12)
scalelabel = sprintf('( x %.2g Pa )', max(abs(pres_s)));
text(ax1, 5, 2.6, scalelabel, 'FontSize', 12)
scalelabel = sprintf('( x %.2g Pa )', max(abs(pres_o)));
text(ax1, 5, 0.8, scalelabel, 'FontSize', 12)

xlim([-10 25])
ylim([-1.5 5.5])
yticks([0 2 4])
yticklabels({'', '', ''})
xlabel('time (s)')

set(ax1, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)

set(gcf, 'Renderer', 'painters')

%% save figure
figname = sprintf('%s.eps', mfilename);
figdisp(figname, [], [], 2, [], 'epstopdf');

if nargout > 0
    varargout{1} = fig;
end
end