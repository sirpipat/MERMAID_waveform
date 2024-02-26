function plotfreqbands(sacfiles, is_pass)
% PLOTFREQBANDS(sacfiles)
%
% plot the seismograms from the sacfiles in various frequency bands to
% figure out which one give the best signal to noise ratio (SNR). SNR is
% computed from the variance of the signal divided the variance of the
% noise.
%
% INPUT:
% sacfiles      cell array of fullfile paths to those SAC files
% is_pass       true  : bandpass                [default]
%               false : original - bandpass
%
% Last modified by sirawich-at-princeton.edu, 02/26/2024

defval('is_pass', true)

% Parameters
SNR_THRESHOLD = 10;

for ii = 1:length(sacfiles)
    try
        % read the seismograms
        [seisdata, hdrdata] = readsac(sacfiles{ii});
        [dt_ref, ~, ~, fs, ~, dts] = gethdrinfo(hdrdata);

        % time relative to picked first arrival time
        t_relative = seconds(dts - dt_ref) - hdrdata.T0;

        % remove instrument response
        pa = real(counts2pa(seisdata, fs, [0.01 0.02 10 20], '', 'sacpz'));

        % create a figure
        figure(7);
        clf
        set(gcf, 'Unit', 'inches', 'Position', [17 12 8 12])

        for jj = 1:7
            % corner frequencies
            fc1 = round(5  / (2 ^ (jj-1)), 2);
            fc2 = round(10 / (2 ^ (jj-1)), 2);

            % filter
            pa_1 = bandpass(pa, fs, fc1, fc2, 2, 2, 'butter', 'linear');
            if ~is_pass
                pa_1 = pa - pa_1;
            end

            % signal-to-noise ratio]
            if is_pass
                halfwin = 2 / fc1;
            else
                halfwin = 2 / 0.05;
            end
            [snr, tt_max] = snrvar(t_relative, pa_1, [-1 1] * halfwin/2, ...
                -60, 60, 1 * halfwin);

            % plot the filtered seismogram
            ax = subplot('Position', [0.13 1-jj*0.118 0.7750 0.085]);
            if snr >= SNR_THRESHOLD
                plot(t_relative, pa_1 .* shanning(length(pa), 0.05, 0), ...
                    'Color', [0.3 0.3 0.3], 'LineWidth', 0.8)
                hold on
                
                % plot the noise window
                t_start = max(-60, tt_max - 1 * halfwin);
                t_end = min(60, tt_max + 1 * halfwin);
                pan = pa_1(and(t_relative >= t_start, t_relative <= tt_max));
                tn = t_relative(and(t_relative >= t_start, t_relative <= tt_max));
                plot(tn, pan, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
                
                % plot the signal window
                pas = pa_1(and(t_relative >= tt_max, t_relative < t_end));
                ts = t_relative(and(t_relative >= tt_max, t_relative < t_end));
                plot(ts, pas, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
            else
                plot(t_relative, pa_1 .* shanning(length(pa), 0.05, 0), ...
                    'Color', [0.7 0.7 0.7], 'LineWidth', 0.8)
                hold on
                
                % plot the noise window
                t_start = max(-60, tt_max - 1 * halfwin);
                t_end = min(60, tt_max + 1 * halfwin);
                pan = pa_1(and(t_relative >= t_start, t_relative <= tt_max));
                tn = t_relative(and(t_relative >= t_start, t_relative <= tt_max));
                plot(tn, pan, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5)
                
                % plot the signal window
                pas = pa_1(and(t_relative >= tt_max, t_relative < t_end));
                ts = t_relative(and(t_relative >= tt_max, t_relative < t_end));
                plot(ts, pas, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5)
            end
            if is_pass
                xlim([-200 200] / 2 ^ (8-jj))
            else
                xlim([-160 160])
            end
            % limit y-axis using only samples within x-limit
            wh = and(t_relative >= ax.XLim(1), t_relative <= ax.XLim(2));
            ylim([-1.2 1.2] * max(abs(pa_1(wh)) .* indeks(shanning(length(pa), 0.05, 0), wh)))
            ax.XTick = (-1:(1/3):1) * ax.XLim(2);
            ax.XTickLabel = {round(ax.XTick(1), 1), [], [], 0, [], [], ...
                round(ax.XTick(end), 1)};
            ylabel('pressure (Pa)')
            grid on
            %nolabels(ax, 1)
            hold on
            vline(ax, tt_max, 'LineStyle', '-', 'LineWidth', 1, ...
                'Color', [0 0.6 1]);
            if jj == 1
                title(sprintf('Event ID: %d, Magnitude: %.2f, Distance: %.2f^{\\circ}, Station: %s', ...
                    hdrdata.USER7, hdrdata.MAG, hdrdata.GCARC, ...
                    hdrdata.KSTNM(ismember(hdrdata.KSTNM, 33:126))));
            end

            % label frequency bands
            axb1 = addbox(ax, [0 0 0.25 0.2]);
            text(0.15, 0.45, sprintf('%.2f--%.2f Hz', fc1, fc2), 'FontSize', 12);

            % label SNR
            axb2 = addbox(ax, [0.75 0 0.25 0.2]);
            if snr > 5
                text(0.2, 0.45, sprintf('SNR = %.d', round(snr)), 'FontSize', 12);
            else
                text(0.2, 0.45, sprintf('SNR = %.2f', snr), 'FontSize', 12);
            end

            axes(axb1)
            set(ax, 'TickDir', 'out', 'FontSize', 12)
        end
        % corner frequencies
        fc1 = 0.05;
        fc2 = 0.10;

        % filter
        pa_1 = bandpass(pa, fs, fc1, fc2, 2, 2, 'butter', 'linear');
        if ~is_pass
            pa_1 = pa - pa_1;
        end

        % signal-to-noise ratio (Adaptive)]
        halfwin = 2 / fc1;
        [snr, tt_max] = snrvar(t_relative, pa_1, [-1 1] * halfwin/2, ...
            -60, 60, 1 * halfwin);

        % plot the filtered seismogram
        ax = subplot('Position', [0.13 0.056 0.7750 0.085]);
        if snr >= SNR_THRESHOLD
            plot(t_relative, pa_1 .* shanning(length(pa), 0.05, 0), ...
                'Color', [0.3 0.3 0.3], 'LineWidth', 0.8)
            hold on
            % plot the noise window
            t_start = max(-60, tt_max - 1 * halfwin);
            t_end = min(60, tt_max + 1 * halfwin);
            pan = pa_1(and(t_relative >= t_start, t_relative <= tt_max));
            tn = t_relative(and(t_relative >= t_start, t_relative <= tt_max));
            plot(tn, pan, 'Color', rgbcolor('1'), 'LineWidth', 1.5)
            % plot the signal window
            pas = pa_1(and(t_relative >= tt_max, t_relative < t_end));
            ts = t_relative(and(t_relative >= tt_max, t_relative < t_end));
            plot(ts, pas, 'Color', rgbcolor('2'), 'LineWidth', 1.5)
        else
            plot(t_relative, pa_1 .* shanning(length(pa), 0.05, 0), ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 0.8)
            hold on
            % plot the noise window
            t_start = max(-60, tt_max - 1 * halfwin);
            t_end = min(60, tt_max + 1 * halfwin);
            pan = pa_1(and(t_relative >= t_start, t_relative <= tt_max));
            tn = t_relative(and(t_relative >= t_start, t_relative <= tt_max));
            plot(tn, pan, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5)
            % plot the signal window
            pas = pa_1(and(t_relative >= tt_max, t_relative < t_end));
            ts = t_relative(and(t_relative >= tt_max, t_relative < t_end));
            plot(ts, pas, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5)
        end
        xlim([-160 160])
        % limit y-axis using only samples within x-limit
        wh = and(t_relative >= ax.XLim(1), t_relative <= ax.XLim(2));
        ylim([-1.2 1.2] * max(abs(pa_1(wh) .* indeks(shanning(length(pa), 0.05, 0), wh))))
        ax.XTick = (-1:(1/3):1) * ax.XLim(2);
        ax.XTickLabel = {round(ax.XTick(1), 1), [], [], 0, [], [], ...
            round(ax.XTick(end), 1)};
        xlabel('time since picked first arrival (s)')
        ylabel('pressure (Pa)')
        grid on
        hold on
        vline(ax, tt_max, 'LineStyle', '-', 'LineWidth', 1, ...
            'Color', [0 0.6 1]);

         % label frequency bands
        axb1 = addbox(ax, [0 0 0.25 0.2]);
        text(0.15, 0.45, sprintf('%.2f--%.2f Hz', fc1, fc2), 'FontSize', 12);

        % label SNR
        axb2 = addbox(ax, [0.75 0 0.25 0.2]);
        if snr > 5
            text(0.2, 0.45, sprintf('SNR = %.d', round(snr)), 'FontSize', 12);
        else
            text(0.2, 0.45, sprintf('SNR = %.2f', snr), 'FontSize', 12);
        end

        axes(axb1)
        set(ax, 'TickDir', 'out', 'FontSize', 12)

        % save a figure
        set(gcf, 'Renderer', 'painters')
        savename = sprintf('%s_%d_%s', mfilename, hdrdata.USER7, ...
            hdrdata.KSTNM(ismember(hdrdata.KSTNM, 33:126)));
        figdisp(savename, [], [], 2, [], 'epstopdf');
    catch ME
        fprintf('%s\n', ME.getReport);
        fprintf('Error occured. Move on to the next iteration.\n');
    end
end
end

% compute signal-to-noise ratio
%
% INPUT:
% t             time
% x             signal
% win_select    time window for time search
% t_begin       begining time for noise window
% t_end         ending time for signal window
% t_length      length of noise and signal windows
%
% OUTPUT:
% s             signal-to-noise ratio
% t_max         time that give the maximum signal-to-noise ratio
function [s, t_max] = snrvar(t, x, win_select, t_begin, t_end, t_length)
tt = t(and(t >= win_select(1), t <= win_select(2)));
ss = zeros(size(tt));
for ii = 1:length(tt)
    ss(ii) = var(x(and(t >= tt(ii), t < min(t_end, tt(ii) + t_length)))) / ...
        var(x(and(t >= max(t_begin, tt(ii) - t_length), t < tt(ii))));
end
s = max(ss);
t_max = tt(ss == s);
end