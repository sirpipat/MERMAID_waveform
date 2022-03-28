function compareresponsefunctions(obsfile, synfile, ddir1, ddir2)
% compareresponsefunctions(obsfile, synfile, ddir1, ddir2)
%
% Plot two synthetic pressure records obtained by convolving the synthetic
% vertical displacement at the ocean bottom with the two response functions
% that convert the displacement at the ocean bottom to acoustic pressure at
% the hydrophone. The goal is to illustrate the difference between the two
% response functions and its consequence.
%
% INPUT
% obsfile       SAC file containing pressure recorded by a hydrophone
% synfile       SAC file containing the synthetic vertical displacement
%               from INSTASEIS
% ddir1         directory to first  SPECFEM2D fluid-solid simulation
% ddir2         directory to second SPECFEM2D fluid-solid simulation
%
% OUTPUT
%
% Last modified by sirawich-at-princeton.edu, 03/28/2022

fcorners = [0.4 1];
window_envelope = [-10 20];
window_waveform = [-5 5];

%% read the data
% read the observed data
[seis_o, hdr_o] = readsac(obsfile);
[dt_ref_o, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);

% time relative to the first P-wave arrival
t_relative = seconds(dts_o - dt_ref_o) - hdr_o.T0;

% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);
pres_o = bandpass(pres_o, fs_o, fcorners(1), fcorners(2), 4, 2, 'butter', 'linear');

% read the synthetic data
[seis_s, hdr_s] = readsac(synfile);
[~, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);


% obtain the response function
[~, ~, t_r1, seis_r1, d] = cctransplot(ddir1, ddir1, [], ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, ...
    [], fs_o, false);

[~, ~, t_r2, seis_r2] = cctransplot(ddir2, ddir2, [], ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, ...
    [], fs_o, false);

% resample to MERMAID datetimes
seis_s_interp = shannon(dts_s, seis_s, dts_o);
t_r1_interp = 0:(1/fs_o):t_r1(end);
seis_r1 = shannon(t_r1, seis_r1, t_r1_interp);
t_r2_interp = 0:(1/fs_o):t_r2(end);
seis_r2 = shannon(t_r2, seis_r2, t_r2_interp);

% read bathymetry
fname1 = cindeks(ls2cell(sprintf('%sDATA/interfaces*.dat', ddir1), 1), 1);
[itfs1, ~] = loadinterfacefile(fname1);

fname2 = cindeks(ls2cell(sprintf('%sDATA/interfaces*.dat', ddir2), 1), 1);
[itfs2, ~] = loadinterfacefile(fname2);

%% convolve
pres_s1 = conv(seis_s_interp, seis_r1);
pres_s1 = pres_s1(1:length(seis_o), 1);
pres_s1 = bandpass(pres_s1, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

pres_s2 = conv(seis_s_interp, seis_r2);
pres_s2 = pres_s2(1:length(seis_o), 1);
pres_s2 = bandpass(pres_s2, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

%% determine the best timeshift
[t_shift1, CCmax1, lags1, cc1, s1] = comparepressure(seis_s, hdr_s, ...
    seis_o, hdr_o, seis_r1, t_r1, window_envelope, window_waveform, ...
    fcorners, false, 1, false);
[t_shift2, CCmax2, lags2, cc2, s2] = comparepressure(seis_s, hdr_s, ...
    seis_o, hdr_o, seis_r2, t_r2, window_envelope, window_waveform, ...
    fcorners, false, 1, ...
    false);

%% correlation coefficient between two synthetics
[t_shift12, CCmax12, lags12, CC12] = ccshift(pres_s1, pres_s2, dt_ref_o, ...
    dt_ref_o, fs_o, seconds(100), 'soft');

%% plot
% color
red   = [0.9 0.0 0.0];
blue  = [0.0 0.1 0.9];
black = [0.4 0.4 0.4];

figure(3)
set(gcf, 'Units', 'inches', 'Position', [9 8 6 8])
clf

% report
ax0 = subplot(27,1,1);
% criteria for CMT solution searching
dt_origin = dt_ref_o + seconds(hdr_o.USER8);
monthname = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
    'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};
tbeg = datenum(dt_origin - minutes(1));
tend = datenum(dt_origin + minutes(1));
mblo = hdr_o.MAG - 0.5;
mbhi = hdr_o.MAG + 0.5;
depmin = hdr_o.EVDP - 50;
depmax = hdr_o.EVDP + 50;

% get CMT solution to report
fname = sprintf('%s%02d.ndk', monthname{dt_origin.Month}, ...
    mod(dt_origin.Year, 100));
[~,~,CMT] = readCMT(fname, strcat(getenv('IFILES'),'CMT'), tbeg, ...
    tend, mblo, mbhi, depmin, depmax, 'hypocenter');

titlestr = ['$$ \textnormal{CMT}-\textnormal{' ...
            CMT.EventName ...
            sprintf('}, M_w=%.2f', hdr_o.MAG) ...
            ',~\Delta=' ...
            sprintf('%.2f', hdr_o.GCARC) ...
            '^{\circ},~\textnormal{depth} = ' ...
            sprintf('%.2f', hdr_o.EVDP) ...
            '~\textnormal{km} $$'];
title(titlestr, 'Interpreter', 'latex')

% LaTeX string reporting the waveform timeshift finding
txtstr = ['$$ \theta =' ...
          sprintf('%.2f', asin(hdr_s.USER9 * 3400 / 6365000) * 180/pi) ...
          '^{\circ}~\textnormal{, avg depth} = ' ...
          sprintf('%.2f', 9600-mean(itfs2{2}.pts(:,2))) ...
          '~\textnormal{m, std depth} = ' ...
          sprintf('%.2f', std(itfs2{2}.pts(:,2))) ...
          '~\textnormal{m, slope} = ' ...
          sprintf('%.2f^{\\circ} $$', ...
          asin(indeks(polyfit(itfs2{2}.pts(:,1), ...
          itfs2{2}.pts(:,2), 1), 1)) * 180/pi)];
[x_pos, y_pos] = norm2trueposition(ax0, -3/40, 3/4);
text(x_pos, y_pos, txtstr, 'FontSize', 12, 'Interpreter', 'latex');

set(ax0, 'FontSize', 12, 'Color', 'none');
ax0.XAxis.Visible = 'off';
        ax0.YAxis.Visible = 'off';

% bathymetry
ax1 = subplot(27,1,[2,8]);
[ax1, hs] = drawbackground(fname2, ax1);
ax1.YTick = 9600 - (8000:-2000:0);
ax1.YTickLabel = string(ax1.YTick-9600);
hold on
vline(ax1, 5000*(1:3), 'LineWidth', 0.5, 'LineStyle', ':', 'Color', [1 1 0.5]);
hline(ax1, ax1.YTick, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', [1 1 0.5]);
plot(itfs2{2}.pts(:,1), itfs2{2}.pts(:,2), 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 1 1])
plot(itfs1{2}.pts(:,1), itfs1{2}.pts(:,2), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r')
xlabel('x (m)')
ylabel('elevation (m)')
%title('bathymetry')

% response function
ax2 = subplot(27,1,[11,14]);
plot(t_r1, seis_r1 / max(abs(seis_r1)) + 1, 'LineWidth', 1, 'Color', red)
hold on
grid on
plot(t_r2, seis_r2 / max(abs(seis_r1)) - 1, 'LineWidth', 1, 'Color', blue)
ylim([-2 2])
legend('flat', 'bathymetry', 'location', 'east')
xlabel('time (s)')
ylabel('response')
title('response function')
set(ax2, 'Box', 'on')

ax3 = subplot(27,1,[17,18]);
plot(t_relative, bandpass(seis_s_interp, fs_o, fcorners(1), ...
    fcorners(2), 4, 2, 'butter', 'linear'), 'LineWidth', 1, 'Color', black)
hold on
grid on
xlim([-10 25])
vline(ax3, 0, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
xlabel('time since first picked arrival (s)')
ylabel('u_z (m)')
title(sprintf('synthetic vertical displacement: bp%.1f-%.1f', fcorners(1), ...
    fcorners(2)))
set(ax3, 'Box', 'on')

ax4 = subplot(27,1,[22,27]);
plot(t_relative, pres_o, 'LineWidth', 1, 'Color', black)
hold on
grid on
plot(t_relative + t_shift1, pres_s1 * s1, 'LineWidth', 1, 'Color', red)
plot(t_relative + t_shift2, pres_s2 * s2, 'LineWidth', 1, 'Color', blue)
xlim([-10 25])
%ylim([-1.1 1.1] * max(max(abs(pres_o)), max(abs(s_e * pres_s1))))
vline(ax4, 0, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
legend('observed', ...
    sprintf('flat : $\\tau^W = %.2f$ s, $X^W = %.2f$', t_shift1, CCmax1), ...
    sprintf('bathymetry : $\\tau^W = %.2f$ s, $X^W = %.2f$', t_shift2, CCmax2), ...
    'Location', 'southoutside', 'Interpreter', 'latex')
xlabel('time since first picked arrival (s)')
ylabel('P (Pa)')
title(sprintf(['acoustic pressure record: bp%.1f-%.1f W^E[%d %d]' ...
    ' W^W[%d %d]'], fcorners(1), fcorners(2), window_envelope(1), ...
    window_envelope(2), window_waveform(1), window_waveform(2)))
set(ax4, 'Box', 'on')

set(gcf, 'Renderer', 'painters')
savename = sprintf('%s_%d_%s.eps', mfilename, hdr_o.USER7, ...
    replace(hdr_o.KSTNM, ' ', ''));
figdisp(savename, [], [], 2, [], 'epstopdf');
end