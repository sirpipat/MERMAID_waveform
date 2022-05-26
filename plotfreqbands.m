function plotfreqbands(sacfiles)
% PLOTFREQBANDS(sacfiles)
%
% plot the seismograms from the sacfiles in various frequency bands to
% figure out which one give the best signal to noise ratio (SNR). SNR is
% computed from the variance of the signal divided the variance of the
% noise.
%
% INPUT:
% sacfiles      cell array of fullfile paths to those SAC files
%
% Last modified by sirawich-at-princeton.edu, 05/18/2022

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

            % signal-to-noise ratio
            % snr = var(pa_1(t_relative >= 0)) / var(pa_1(t_relative < 0));
            tt = t_relative(and(t_relative >= -2/sqrt(fc1), t_relative <= 2/sqrt(fc1)));
            snrs = zeros(size(tt));
            for kk = 1:length(tt)
                snrs(kk) = var(pa_1(and(t_relative >= tt(kk), t_relative < 60))) / ...
                    var(pa_1(and(t_relative < tt(kk), t_relative >= -80)));
            end
            snr = max(snrs);
            tt_max = tt(snrs == snr);

            % plot the filtered seismogran
            ax = subplot('Position', [0.13 1-jj*0.118 0.7750 0.085]);
            plot(t_relative, pa_1 .* shanning(length(pa), 0.05, 0), 'k', ...
                'LineWidth', 0.8)
            xlim([-100 100] / 2 ^ (4-jj/2))
            ylim([-1 1] * max(abs(pa_1 .* shanning(length(pa), 0.05, 0))))
            ylabel('pressure (Pa)')
            grid on
            %nolabels(ax, 1)
            hold on
            vline(ax, tt_max, 'LineStyle', '-', 'LineWidth', 1, ...
                'Color', [0 0.6 1]);
            if jj == 1
                title(sprintf('Event ID: %d, Magnitude: %.2f, Distance: %.2f^{\\circ}, Station: %s', ...
                    hdrdata.USER7, hdrdata.MAG, hdrdata.GCARC, ...
                    replace(hdrdata.KSTNM, ' ', '')));
            end

            % label frequency bands
            axb1 = addbox(ax, [0 0 0.25 0.2]);
            text(0.15, 0.45, sprintf('%.2f--%.2f Hz', fc1, fc2), 'FontSize', 12);

            % label SNR
            axb2 = addbox(ax, [0.75 0 0.25 0.2]);
            text(0.2, 0.45, sprintf('SNR = %.2f', snr), 'FontSize', 12);

            axes(axb1)
            set(ax, 'TickDir', 'both', 'FontSize', 12)
        end
        % corner frequencies
        fc1 = 0.05;
        fc2 = 0.10;

        % filter
        pa_1 = bandpass(pa, fs, fc1, fc2, 2, 2, 'butter', 'linear');

        % signal-to-noise ratio (Adaptive)
        tt = t_relative(and(t_relative >= -20, t_relative <= 20));
        snrs = zeros(size(tt));
        for kk = 1:length(tt)
            snrs(kk) = var(pa_1(and(t_relative >= tt(kk), t_relative < 60))) / ...
                var(pa_1(and(t_relative < tt(kk), t_relative >= -60)));
        end
        snr = max(snrs);
        tt_max = tt(snrs == snr);

        % plot the filtered seismogran
        ax = subplot('Position', [0.13 0.056 0.7750 0.085]);
        plot(t_relative, pa_1 .* shanning(length(pa), 0.05, 0), 'k', ...
            'LineWidth', 0.8)
        xlim([-100 100])
        ylim([-1 1] * max(abs(pa_1 .* shanning(length(pa), 0.05, 0))))
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
        text(0.2, 0.45, sprintf('SNR = %.2f', snr), 'FontSize', 12);

        axes(axb1)
        set(ax, 'TickDir', 'both', 'FontSize', 12)

        % save a figure
        set(gcf, 'Renderer', 'painters')
        savename = sprintf('%s_%d_%s', mfilename, hdrdata.USER7, ...
            replace(hdrdata.KSTNM, ' ', ''));
        figdisp(savename, [], [], 2, [], 'epstopdf');
    catch ME
        fprintf('%s\n', ME.getReport);
        fprintf('Error occured. Move on to the next iteration.\n');
    end
end