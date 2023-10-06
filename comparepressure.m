function [t_shift, CCmax, lag, cc, s] = ...
    comparepressure(seis_s, hdr_s, seis_o, hdr_o, seis_r, t_r, ...
    envelope_window, waveform_window, fcorners, plt, numpicks, is_rmnoise)
% [t_shift, CCmax, lag, cc, s]  = ...
%     COMPAREPRESSURE(seis_s, hdr_s, seis_o, hdr_o, seis_r, t_r, ...
%                     [E_B E_E], [W_B W_E], [lo hi], plt, numpicks, ...
%                     is_rmnoise)
%
% interpolates synthetic and response. Then, convolves the two and compares
% with the observed with instrument response removed.
%
% INPUT:
% seis_s        synthetic seismogram
% hdr_s         SAC header of the synthetic seismogram
% seis_o        observed seismogram (pressure record at hydrophone)
% hdr_o         SAC header of the observed seismogram
% seis_r        receiver function
% t_r           time of the receiver function
% [E_B E_E]     begin and end of the window for envelope correlation
%               [Default: [-20 20]]
% [W_B W_E]     begin and end of the window for waveform correlation
%               [Default: [-5 5]]
% [lo hi]       corner frequencies for bandpass filtering 
%               [Default: [0.4 1] Hz]
% plt           whether to plot or not [Default: true]
% numpicks      the number of bestshift candidates. Only work with plots.
%               [Default: 1]
% is_rmnoise    whether to use REMOVENOISE to remove the noise or not
%               This is experimental. [default: false]
% method        options for method to determine the starting point for
%               waveform correlation [default: 1]
%               1 - rising signal in the synthetic pressure waveform
%               2 - ray-theory arrival time with transfer function included
%               3 - envelope correlation
%
% OUTPUT:
% t_shift       N best time shifts where CCs are at maximum peak
% CCmax         N maximum correlation coefficients
% lag           Vector of all time shifts
% CC            Vector of CC for every time shift in lag
% s             Scaling to minimize the misfit
%
% Last modified by sirawich-at-princeton.edu, 10/06/2023

defval('envelope_window', [-20 20])
defval('waveform_window', [-5 5])
defval('fcorners', [0.4 1])
defval('plt', true)
defval('numpicks', 1)
defval('is_rmnoise', false)
defval('method', 2)

%% Part 1: compute the synthetic pressure
% sampling rate
[dt_ref_s, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);
[dt_ref_o, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);
fs_r = 1 / (t_r(2) - t_r(1));

% downsample the response to about 100 Hz
d_factor = ceil(fs_r / 100);
if d_factor > 1
    seis_r = lowpass(seis_r, fs_r, fs_r / d_factor, 2, 2, 'butter', 'linear');
    seis_r = downsample(seis_r, d_factor);
    t_r = downsample(t_r, d_factor);
    fs_r = fs_r / d_factor;
end

% resample to MERMAID datetimes
seis_s = shannon(dts_s, seis_s, dts_o);
tr_r_interp = 0:(1/fs_o):t_r(end);
seis_r = shannon(t_r, seis_r, tr_r_interp);

% convolve
pres_s = conv(seis_s, seis_r) ;
pres_s = pres_s(1:length(seis_o), 1);
pres_s = bandpass(pres_s, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

% remove instrument response and filter
pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
pres_o = real(pres_o);
pres_o = bandpass(pres_o, fs_o, fcorners(1), fcorners(2), 4, 2, ...
    'butter', 'linear');

% EXPERIMENT: use removenoise to remove noise prior to the arrival
if is_rmnoise
    dt_arrival = dt_ref_o + seconds(hdr_o.T0);
    pres_o_ori = pres_o;
    pres_o = removenoise(dts_o, pres_o, fs_o, dt_arrival, 100);
end

%% Part 2: compute the cross correlation
% Epsilon for comparison (i.e. <, <=, >, >= etc.)
ep = seconds((1 / fs_o) / 100);

% cut the short windows to do envelope cross correlation
dt_start1 = dt_ref_o + seconds(hdr_o.T0 + envelope_window(1));
dt_end1 = dt_ref_o + seconds(hdr_o.T0 + envelope_window(2));
pres_o1 = pres_o(and(geq(dts_o, dt_start1, ep), leq(dts_o, dt_end1, ep)));
pres_s1 = pres_s(and(geq(dts_o, dt_start1, ep), leq(dts_o, dt_end1, ep)));

% CURRENT METHOD:
% try figuring the starting point for waveform cross-correlation using
% 5% of the maximum
if method == 1
    max_amplitude = max(abs(pres_s(dts_o - dt_ref_o > seconds(50))));
    dt_arrival = indeks(dts_o(and(abs(pres_s) >= 0.05 * max_amplitude, ...
        dts_o - dt_ref_o > seconds(50))), 1);
    best_lags_time_e = seconds(dt_ref_o - dt_arrival) + hdr_o.T0;
elseif method == 2
    % EXPERIMENTAL METHOD: 
    % use expected arrival time as the starting point
    % Try to align the starting piont with the timestamp of the observed
    % pressure waveform
    
    % arrival times with respect to beginning of the observed waveform
    t_P_arrival_guess = seconds(dt_ref_s - dts_o(1)) + hdr_s.T0 + ...
        tr_r_interp(seis_r == max(seis_r));
    t_P_arrival_record = seconds(dt_ref_o - dts_o(1)) + hdr_o.T0;
    % round to the exact timestamp of the observed pressure waveform
    idt_guess = knnsearch(seconds(dts_o - dts_o(1)), ...
        t_P_arrival_guess);
    idt_record = knnsearch(seconds(dts_o - dts_o(1)), ...
        t_P_arrival_record);
    
    best_lags_time_e = seconds(dts_o(idt_record) - dts_o(idt_guess));
else
    % assigned the starting point to be a very large number to force
    % envelope correlation
    best_lags_time_e = inf;
end
% determine maximum accepted best lag time for envelope cross-correlation
% which is the maximum of 5 seconds or 2 percent of the travel time
maxmargin_e = seconds(max(5, 0.02 * hdr_o.USER5));

if abs(best_lags_time_e) > seconds(maxmargin_e)
% first correlate by the envelope
    [best_lags_time_e, ~, lags_time_e, cc_e, s_e] = ccscale(pres_o1, pres_s1, ...
        dt_begin_o, dt_begin_o, fs_o, maxmargin_e, 'soft', true);
else
    [~, ~, lags_time_e, cc_e, s_e] = ccscale(pres_o1, pres_s1, ...
        dt_begin_o, dt_begin_o, fs_o, maxmargin_e, 'soft', true);
end
% cut the short windows to do waveform cross correlation
dt_start2 = dt_ref_o + seconds(hdr_o.T0 + waveform_window(1));
dt_end2 = dt_ref_o + seconds(hdr_o.T0 + waveform_window(2));
pres_o2 = pres_o(and(geq(dts_o, dt_start2, ep), leq(dts_o, dt_end2, ep)));
pres_s2 = pres_s(and(geq(dts_o, ...
    dt_start2 - seconds(best_lags_time_e) + seconds(waveform_window(1)), ep), ...
    leq(dts_o, dt_end2 - seconds(best_lags_time_e) + seconds(waveform_window(2)), ep)));

% then correlate by the waveform
maxmargin_w = seconds(5);%seconds(waveform_window(2) - waveform_window(1)) / 2;
[best_lags_time, ~, lags_time, cc, s] = ccscale(pres_o2, pres_s2, ...
    dt_start2, dt_start2 + seconds(waveform_window(1)), fs_o, ...
    maxmargin_w, 'soft', false);

best_lags_time = best_lags_time + best_lags_time_e;
lags_time = lags_time + best_lags_time_e;

% identify NUMPICKS of peaks
[pks, locs] = findpeaks(cc, lags_time);
[pks, isort] = sort(pks, 'descend');
% list of time shifts that give maximum CC peaks (sorted by CC)
locs = locs(isort);

% variables for the output
t_shift = locs(1:numpicks);
CCmax = indeks(cc(logical(sum(lags_time == locs',1))), isort(1:numpicks))';
lag = lags_time;

%% Part 2B: track correlation coefficient as the waveform window grows
if false
    ww = [0, 5; 0, 10; 0, 15; 0, 20; 0, 25] + waveform_window;
    win = [waveform_window(2); ww(:,2)] - waveform_window(1);
    CCmax_after = zeros(size(ww, 1), numpicks);
    for ii = 1:size(ww, 1)
        dt_start3 = dt_ref_o + seconds(hdr_o.T0 + ww(ii, 1));
        dt_end3 = dt_ref_o + seconds(hdr_o.T0 + ww(ii, 2));
        pres_o3 = pres_o(and(geq(dts_o, dt_start3, ep), leq(dts_o, dt_end3, ep)));
        pres_s3 = pres_s(and(geq(dts_o, ...
            dt_start3 - seconds(best_lags_time_e), ep), ...
            leq(dts_o, dt_end3 - seconds(best_lags_time_e), ep)));

        [blt3, ~, lt3, cc3, s3] = ccscale(pres_o3, pres_s3, ...
            dt_start3, dt_start3, fs_o, ...
            seconds(ww(ii, 2) - ww(ii, 1)) / 2, 'soft', false);

        lt3 = lt3 + best_lags_time_e;
        for jj = 1:numpicks
            CCmax_after(ii, jj) = cc3(lt3 == t_shift(jj));
        end
    end

    CCmax_all = [CCmax; CCmax_after];

    figure(2)
    clf
    plot(win, CCmax_all, '-o', 'LineWidth', 1)
    grid on
    xlabel('waveform window length')
    ylabel('correlation coefficient')
    set(gca, 'TickDir', 'both', 'FontSize', 12)
    legend(sprintf('%.2f s', t_shift(1)), sprintf('%.2f s', t_shift(2)), ...
        sprintf('%.2f s', t_shift(3)), sprintf('%.2f s', t_shift(4)), ...
        sprintf('%.2f s', t_shift(5)));
    title(sprintf('Event ID: %d, Magnitude: %.2f, Distance: %.2f\\circ, Station: %s', ...
        hdr_o.USER7, hdr_o.MAG, hdr_o.GCARC, replace(hdr_o.KSTNM, ' ', '')));

    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%d_%s_extended.eps', mfilename, hdr_o.USER7, ...
        replace(hdr_o.KSTNM, ' ', ''));
    figdisp(savename, [], [], 2, [], 'epstopdf');
end

%% Part 3: plot the result
if plt
    % criteria for CMT solution searching
    dt_origin = dt_ref_o + seconds(hdr_o.USER8);
    monthname = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};
    tbeg = datenum(dt_origin - minutes(1));
    tend = datenum(dt_origin + minutes(1));
    mblo = hdr_o.MAG - 1.0;
    mbhi = hdr_o.MAG + 1.0;
    depmin = hdr_o.EVDP - 50;
    depmax = hdr_o.EVDP + 50;
    
    % get CMT solution to report
    fname = sprintf('%s%02d.ndk', monthname{dt_origin.Month}, ...
        mod(dt_origin.Year, 100));
    [~,~,CMT] = readCMT(fname, strcat(getenv('IFILES'),'CMT'), tbeg, ...
        tend, mblo, mbhi, depmin, depmax, 'hypocenter');
    
    for ii = 1:min(numpicks, length(pks))
        figure(1)
        clf
        set(gcf, 'Units', 'inches', 'Position', [0 6 6 8])
        
        % axes handling the title
        ax0 = subplot('Position', [0.08 0.92 0.86 0.03]);
        % LaTeX string describing the information of the earthquake
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
        txtstr = ['$$ \tau^W =' ...
                  sprintf('%.2f', locs(ii)) ...
                  '~\textnormal{s},~X^W = ' ...
                  sprintf('%.2f', pks(ii)) ...
                  ',~s^E = ' ...
                  sprintf('%.2f', s_e) ...
                  ',~s^W = ' ...
                  sprintf('%.2f $$', s)];
        [x_pos, y_pos] = norm2trueposition(ax0, 1/10, 3/4);
        text(x_pos, y_pos, txtstr, 'FontSize', 12, 'Interpreter', 'latex');
        
        set(ax0, 'FontSize', 12, 'Color', 'none');
        ax0.XAxis.Visible = 'off';
        ax0.YAxis.Visible = 'off';

        % plot two pressure records: observed vs synthetic
        ax1 = subplot('Position', [0.08 0.74 0.86 0.16]);
        cla
        plot(seconds(dts_o' - dt_ref_o) - hdr_o.T0, pres_o, 'k', 'LineWidth', 0.5)
        hold on
        plot(seconds(dts_o' - dt_ref_o) - hdr_o.T0 + locs(ii), ...
            s_e * pres_s, 'b', 'LineWidth', 1)
        grid on
        xlim([-15 35])
        ylim([-1.1 1.1] * max(max(abs(pres_o1)), max(abs(s_e * pres_s1))))
        vline(ax1, 0, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
        legend('observed', 'synthetic', ...
            'Location', 'northwest')
        xlabel('time since first picked arrival (s)')
        ylabel('acoustic pressure (Pa)')
        title('pressure record: waveform')
        set(ax1, 'Box', 'on')

        % plot two envelope records: observed vs synthetic
        ax2 = subplot('Position', [0.08 0.50 0.86 0.16]);
        cla
        plot(seconds(dts_o' - dt_ref_o) - hdr_o.T0, envelope(pres_o), ...
            'Color', 'k', 'LineWidth', 0.5)
        hold on
        plot(seconds(dts_o' - dt_ref_o) - hdr_o.T0 + locs(ii), ...
            s_e * envelope(pres_s), 'Color', [0 0 1], 'LineWidth', 1)
        grid on
        xlim([-15 35])
        ylim([-0.1 1.1] * max(max(envelope(pres_o1)), ...
            max(s_e * envelope(pres_s1))))
        vline(ax2, 0, 'LineWidth', 2, ...
            'LineStyle', '--', 'Color', [0.1 0.8 0.1]);
        legend('observed', 'synthetic', ...
            'Location', 'northwest')
        xlabel('time since first picked arrival (s)')
        ylabel('acoustic pressure (Pa)')
        title('pressure record: envelope')
        set(ax2, 'Box', 'on')
        
        % correlation coefficient: waveform
        ax3 = subplot('Position', [0.08 0.28 0.86 0.14]);
        cla
        hold on
        xlim([-20 20])
        ylim([-1 1] * 1.25 * max(abs(cc)))
        vline(ax3, locs(ii), 'LineStyle', '--', ...
            'Color', rgbcolor('2'), 'LineWidth', 2);
        hline(ax3, pks(ii), 'LineStyle', '--', ...
            'Color', rgbcolor('2'), 'LineWidth', 2);
        plot(lags_time, cc, 'k')
        grid on
        title('correlation coefficient: waveform')
        xlabel('timeshift (s)')
        ylabel('correlation coefficient')
        set(ax3, 'Box', 'on')
        
        % correlation coefficient: envelope
        ax4 = subplot('Position', [0.08 0.06 0.86 0.14]);
        cla
        hold on
        xlim([-20 20])
        ylim([-0.05 1.05])
        vline(ax4, locs(ii), 'LineStyle', '--', ...
            'Color', rgbcolor('2'), 'LineWidth', 2);
        plot(lags_time_e, cc_e, 'k')
        grid on
        title('correlation coefficient: envelope')
        xlabel('timeshift (s)')
        ylabel('correlation coefficient')
        set(ax4, 'Box', 'on')

        % save figure
        set(gcf, 'Renderer', 'painters')
        if ii == 0
            savename = sprintf('%s_%d_%s.eps', mfilename, hdr_o.USER7, ...
                replace(hdr_o.KSTNM, ' ', ''));
        else
            savename = sprintf('%s_%d_%s_rank%d.eps', mfilename, hdr_o.USER7, ...
                replace(hdr_o.KSTNM, ' ', ''), ii);
        end
        figdisp(savename, [], [], 2, [], 'epstopdf');
    end
end
end

function r = leq(a, b, ep)
r = or(a < b, abs(a - b) < ep);
end

function r = geq(a, b, ep)
r = or(a > b, abs(a - b) < ep);
end