function [t_shift1, t_shift2, CCmax1, CCmax2, bath1, bath2, fcorners, ...
    snr, s1, s2] = ...
    compareresponsefunctions(obsfile, synfile, ddir1, ddir2, opt, plt, la)
% [t_shift1, t_shift2, CCmax1, CCmax2, bath1, bath2, fcorners, ...
%     snr, s1, s2] = ...
%     compareresponsefunctions(obsfile, synfile, ddir1, ddir2, opt, plt, la)
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
% opt           corner frequency options
%                   1             fixed at 1-2 Hz
%                   2             selected by FREQSELECT  [Defatult]
%                   [fc1 fc2]     user defined corner frequencies
% plt           whether to plot or not [Defatult: true]
% la            label in the plot referred as {first_record, second_record}
%               [Default: {'flat', 'bathymetry'}
%
% OUTPUT
% t_shift1      best time shifts where CCs are at maximum peak for first  
%               SPECFEM2D fluid-solid simulation
% t_shift2      best time shifts where CCs are at maximum peak for second 
%               SPECFEM2D fluid-solid simulation
% CCmax1        maximum correlation coefficients for first  SPECFEM2D 
%               fluid-solid simulation
% CCmax2        maximum correlation coefficients for second SPECFEM2D 
%               fluid-solid simulation
% bath1         bathymetry profile [x, z] for first  SPECFEM2D 
%               fluid-solid simulation
% bath2         bathymetry profile [x, z] for second SPECFEM2D 
%               fluid-solid simulation
% fcorners      corner frequencies used for comparing synthetic and
%               observed acoustic pressures
% snr           best-signal to noise ratio given fcorners (only available
%               for fcorners chosen by FREQSELECT; opt==2. Otherwise, snr 
%               is set to NaN)
% s1            scaling to minimize the misfit for first  SPECFEM2D 
% s2            scaling to minimize the misfit for second SPECFEM2D 
%
% Last modified by sirawich-at-princeton.edu, 12/04/2023

defval('fopt', 2)
defval('plt', true)
defval('la', {'flat', 'bathymetry'})

window_envelope = [-20 20];
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

% determine corner frequency
if length(opt) == 1 && opt == 1
    fcorners = [1 2];
    snr = nan;
elseif length(opt) == 1 && opt == 2
    [fcorners, snr] = freqselect(t_relative, pres_o, fs_o, false);
    fcorners(1) = max(fcorners(1), 0.05);
    fcorners(2) = min(fcorners(2), 2);
elseif length(opt) == 2
    fcorners = opt;
    snr = nan;
else
    fprintf('ERROR: invalid corner frequency option\n');
    return
end

% filter
pres_o = bandpass(pres_o, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

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
bath1 = itfs1{2}.pts;

fname2 = cindeks(ls2cell(sprintf('%sDATA/interfaces*.dat', ddir2), 1), 1);
[itfs2, ~] = loadinterfacefile(fname2);
bath2 = itfs2{2}.pts;

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
    fcorners, false, 1, false);

%% correlation coefficient between two synthetics
[t_shift12, CCmax12, lags12, CC12] = ccshift(pres_s1, pres_s2, dt_ref_o, ...
    dt_ref_o, fs_o, seconds(100), 'soft');

%% plot
if plt
    % color
    red   = [0.9 0.0 0.0];
    blue  = [0.0 0.1 0.9];
    black = [0.4 0.4 0.4];

    figure(3)
    set(gcf, 'Units', 'inches', 'Position', [0 8 6 8])
    clf

    %% report
    ax0 = subplot('Position', [0.1300 0.9400 0.7750 0.0220]);
    % criteria for CMT solution searching
    dt_origin = dt_ref_o + seconds(hdr_o.USER8);
    monthname = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};
    tbeg = datenum(dt_origin - minutes(1));
    tend = datenum(dt_origin + minutes(1));
    mblo = hdr_o.MAG - 0.51;
    mbhi = hdr_o.MAG + 0.51;
    depmin = hdr_o.EVDP - 50;
    depmax = hdr_o.EVDP + 50;

    % get CMT solution to report
    fname = sprintf('%s%02d.ndk', monthname{dt_origin.Month}, ...
        mod(dt_origin.Year, 100));
    [~,~,CMT] = readCMT(fname, strcat(getenv('IFILES'),'CMT'), tbeg, ...
        tend, mblo, mbhi, depmin, depmax, 'hypocenter');

    txtstr0 = ['$$ \textnormal{' ...
               CMT.EventName ...
               sprintf(', lat} = %.2f', hdr_o.EVLA) ...
               '^{\circ},~\textnormal{lon} = ' ...
               sprintf('%.2f', hdr_o.EVLO) ...
               '^{\circ},~\textnormal{depth} = ' ...
               sprintf('%.2f', hdr_o.EVDP) ...
               '~\textnormal{km},~' ...
               sprintf('M_w=%.2f', hdr_o.MAG) ...
               ' $$'];
    [x_pos0, y_pos0] = norm2trueposition(ax0, -4/40, 1.5);
    text(x_pos0, y_pos0, txtstr0, 'FontSize', 11, 'Interpreter', 'latex');
    % LaTeX string reporting the receiver
    txtstr1 = ['$$ \textnormal{MH-' ...
               replace(hdr_o.KSTNM, ' ', '') ...
               ', lat} = ' ...
               sprintf('%.2f', hdr_o.STLA) ...
               '^{\circ},~\textnormal{lon} = ' ...
               sprintf('%.2f', hdr_o.STLO) ...
               '^{\circ},~\textnormal{elev} = ' ...
               sprintf('%.2f', -hdr_o.STDP) ...
               '~\textnormal{m},~\Delta=' ...
               sprintf('%.2f', hdr_o.GCARC) ...
               '^{\circ} $$'];
    [x_pos1, y_pos1] = norm2trueposition(ax0, -4/40, 0.7);
    text(x_pos1, y_pos1, txtstr1, 'FontSize', 11, 'Interpreter', 'latex');
    % LaTeX string reporting the waveform timeshift finding
    txtstr2 = ['$$ \theta_i =' ...
               sprintf('%.2f', asin(hdr_s.USER9 * 3400 / 6365000) * 180/pi) ...
               '^{\circ}~\textnormal{, avg depth} = ' ...
               sprintf('%.2f', 9600-mean(itfs2{2}.pts(:,2))) ...
               '~\textnormal{m, std depth} = ' ...
               sprintf('%.2f', std(itfs2{2}.pts(:,2))) ...
               '~\textnormal{m, slope} = ' ...
               sprintf('%.2f^{\\circ} $$', ...
               atan(indeks(polyfit(itfs2{2}.pts(:,1), ...
               itfs2{2}.pts(:,2), 1), 1)) * 180/pi)];
    [x_pos2, y_pos2] = norm2trueposition(ax0, -4/40, -0.1);
    text(x_pos2, y_pos2, txtstr2, 'FontSize', 11, 'Interpreter', 'latex');

    set(ax0, 'FontSize', 11, 'Color', 'none');
    set(ax0.Title, 'FontSize', 11)
    ax0.XAxis.Visible = 'off';
    ax0.YAxis.Visible = 'off';

    %% bathymetry
    ax1 = subplot('Position', [0.1300 0.7100 0.7750 0.2050]);
    [ax1, hs] = drawbackground(fname2, ax1);
    ax1.YTick = 9600 - (8000:-2000:0);
    ax1.YTickLabel = string(ax1.YTick-9600);
    hold on
    vline(ax1, 5000*(1:3), 'LineWidth', 0.5, 'LineStyle', ':', ...
        'Color', [1 1 0.5]);
    hline(ax1, ax1.YTick, 'LineWidth', 0.5, 'LineStyle', ':', ...
        'Color', [1 1 0.5]);
    plot(itfs2{2}.pts(:,1), itfs2{2}.pts(:,2), 'LineStyle', '-', ...
        'LineWidth', 2, 'Color', [0.6 1 1])
    plot(itfs1{2}.pts(:,1), itfs1{2}.pts(:,2), 'LineStyle', '--', ...
        'LineWidth', 2, 'Color', 'r')
    xlabel('x (m)')
    ylabel('elevation (m)')
    set(ax1, 'FontSize', 8, 'Box', 'on', 'TickDir', 'both')
    
    axb1 = boxedlabel(ax1, 'northwest', 0.18, [], 'a');

    %% response function
    ax2 = subplot('Position', [0.1300 0.5175 0.7750 0.1135]);
    plot(t_r1, seis_r1 / max(abs(seis_r1)) + 1, 'LineWidth', 1, ...
        'Color', red)
    hold on
    grid on
    plot(t_r2, seis_r2 / max(abs(seis_r1)) - 1, 'LineWidth', 1, ...
        'Color', blue)
    ylim([-2 2])
    legend(la{1}, la{2}, 'location', 'east')
    xlabel('time (s)')
    ylabel('response')
    title('response function', 'Interpreter', 'latex', 'FontSize', 11)
    set(ax2, 'FontSize', 8, 'Box', 'on', 'TickDir', 'both')
    axb2 = boxedlabel(ax2, 'northwest', 0.18, [], 'b');

    %% synthetic vertical displacement at the ocean bottom
    ax3 = subplot('Position', [0.1300 0.3970 0.7750 0.0525]);
    plot(t_relative, bandpass(seis_s_interp, fs_o, fcorners(1), ...
        fcorners(2), 4, 2, 'butter', 'linear'), 'LineWidth', 1, ...
        'Color', black)
    hold on
    grid on
    xlim([-10 25])
    ylabel('u_z (m)')
    set(ax3, 'FontSize', 8, 'Box', 'on', 'TickDir', 'both', ...
        'XTickLabel', {})
    title(sprintf('synthetic vertical displacement: bp%.1f-%.1f', ...
        fcorners(1), fcorners(2)), 'Interpreter', 'latex', 'FontSize', 11)
    axb3 = boxedlabel(ax3, 'northwest', 0.18, [], 'c');

    %% pressure recorded by the hydrophone
    ax4 = subplot(27,1,[20,25]);
    p1 = plot(t_relative, pres_o, 'LineWidth', 1, 'Color', black);
    hold on
    grid on
    p2 = plot(t_relative + t_shift1, pres_s1 * s1, 'LineWidth', 0.5, ...
        'Color', red);
    p3 = plot(t_relative + t_shift2, pres_s2 * s2, 'LineWidth', 1, ...
        'Color', blue);
    xlim([-10 25])
    wh = and(t_relative >= window_waveform(1), ...
        t_relative <= window_waveform(2));
    wh1 = and(t_relative + t_shift1 >= window_waveform(1), ...
        t_relative + t_shift1 <= window_waveform(2));
    wh2 = and(t_relative + t_shift2 >= window_waveform(1), ...
        t_relative + t_shift2 <= window_waveform(2));
    ylim([-1.2 1.2] * max([max(abs(pres_o(wh))), ...
                           max(abs(s1 * pres_s1(wh1))), ...
                           max(abs(s2 * pres_s2(wh2))) ...
                          ]));
    
    label2 = sprintf('$$ \\textnormal{%s} : \\tau^W = %.2f~\\textnormal{s, X}^W = %.2f $$', la{1}, t_shift1, CCmax1);
    label3 = sprintf('$$ \\textnormal{%s} : \\tau^W = %.2f~\\textnormal{s, X}^W = %.2f $$', la{2}, t_shift2, CCmax2);
    legend([p1, p2, p3], {'observed', label2, label3}, ...
        'Location', 'southoutside', 'Interpreter', 'latex')
    xlabel('time since first picked arrival (s)')
    ylabel('P (Pa)')
    set(ax4, 'FontSize', 8, 'Box', 'on', 'TickDir', 'both', ...
        'Color', 'none', 'Position', [0.1300 0.3080 0.7750 0.0525])
    
    % high light the waveform window
    ax4s = doubleaxes(ax4);
    [xbox, ybox] = boxcorner(window_waveform, 0.98 * ax4.YLim);
    pgon = polyshape(xbox, ybox);
    bx = plot(ax4s, pgon, 'FaceColor', [1 0.9 0.4], 'FaceAlpha', 0.4, ...
        'EdgeAlpha', 0);
    ax4s.XAxis.Visible = 'off';
    ax4s.YAxis.Visible = 'off';
    set(ax4s, 'Box', 'on', 'TickDir', 'both', 'XLim', ax4.XLim, 'YLim', ...
        ax4.YLim, 'Position', ax4.Position);
    
    % fix the title
    axes(ax4);
    title(sprintf(['acoustic pressure record: bp%.1f-%.1f W$^E$[%d %d]' ...
        ' W$^W$[%d %d]'], fcorners(1), fcorners(2), window_envelope(1), ...
        window_envelope(2), window_waveform(1), window_waveform(2)), ...
        'Interpreter', 'latex', 'FontSize', 11)
    
    axb4 = boxedlabel(ax4, 'northwest', 0.18, [], 'd');
    
    %% cross correlation between the observed and first synthetic
    ax5 = subplot('Position', [0.1300 0.1180 0.7750 0.0525]);
    plot(lags1, cc1, 'Color', red, 'LineWidth', 0.5)
    grid on
    ylim([-1 1])
    ylabel('X^W')
    set(ax5, 'FontSize', 8, 'Box', 'on', 'TickDir', 'both', ...
        'XTickLabel', {})
    title(sprintf('correlation coefficients (red - %s, blue - %s)', ...
        la{1}, la{2}), 'Interpreter', 'latex', 'FontSize', 11)

    
    % cross correlation between the observed and second synthetic
    ax6 = subplot('Position', [0.1300 0.0500 0.7750 0.0525]);
    plot(lags2, cc2, 'Color', blue, 'LineWidth', 0.5)
    grid on
    ylim([-1 1])
    xlabel('time shift (s)')
    ylabel('X^W')
    set(ax6, 'FontSize', 8, 'Box', 'on', 'TickDir', 'both')
    
    % align the two cross correlation plots
    extend = [min(ax5.XLim(1), ax6.XLim(1)) max(ax5.XLim(2), ax6.XLim(2))];
    ax5.XLim = extend;
    ax6.XLim = extend;
    
    vline(ax5, t_shift1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.);
    hline(ax5, CCmax1, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.);
    
    vline(ax6, t_shift2, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.);
    hline(ax6, CCmax2, 'LineStyle', ':', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.);
    
    axb5 = boxedlabel(ax5, 'northwest', 0.18, [], 'e');
    axb6 = boxedlabel(ax6, 'northwest', 0.18, [], 'f');

    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%d_%s.eps', mfilename, hdr_o.USER7, ...
        replace(hdr_o.KSTNM, ' ', ''));
    figdisp(savename, [], [], 2, [], 'epstopdf');
end
end