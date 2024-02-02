function noise = mermaidnoise(sacfiles, cutoff, fcorners, plt)
% noise = MERMAIDNOISE(sacfiles, cutoff, fcorners, plt)
%
% Computes the noise level in the MERMAID observations in two ways: (1) 
% standard deviation of the time-series data and (2) power spectral density
% of (a) instrument response removed acoustic pressure seismogram from the
% start to 2 seconds prior to the picked arrival time, (b) the same as (a)
% but highpass filtered at 0.4 Hz, and (c) the entire seismogram from (a)
% hipass at 0.4 Hz and bandstop at the specified corner frequency pair.
%
% INPUT:
% sacfiles      cell array contain the SAC filenames
% cutoff        cut off time relative to picked arrival time in seconds
%               It must be an [n 2] matrix where n is either
%                   0 for default cut off itme at -2 seconds (2 seconds
%                     before the picekd arrival time)
%                   1 for the same cut off time
%                   length(sacfiles) for a cut off time for each 
%                         observaitons
% fcorners      corner frequency pairs. It must be an [n 2] matrix where n
%               is either
%                   0 for default corner frequency pairs of [0.4 2] Hz
%                   1 for the same corner frequency pair for all
%                         observations
%                   length(sacfiles) for a corner frequency pairs for each
%                         observaitons
% plt           whether to plot or not [default: true]
%
% OUTPUT:
% noise         struct with a following fields
%       std_all - standard deviation of the instrument response 
%                 removed, acoustic pressure seismogram
%       std_before - standard deviation before the cut off time
%       std_after - standard deviation after the cut off time
%       psd_section_raw - power spectral density of instrument response 
%                         removed, acoustic pressure seismogram (a)
%       psd_section_hipass - power spectral density of hipass seismogram (b)
%       psd_section_bandstop - power spectral density of bandstop seismogram (c)
%       F - list of frequencies for the power spectral density outputs
%
% SEE ALSO:
% FREQSELECT, PCHAVE
%
% Last modified by sirawich-at-princeton.edu, 02/02/2024

% bad value for a SAC header field
badval = -12345;

% cut off time for the noise before earthquake arrival only
defval('cutoff', -2)
defval('fcorners', [0.4 2])

if isscalar(cutoff)
    cutoff = repmat(cutoff, length(sacfiles), 1);
end
if size(fcorners, 1) == 1
    fcorners = repmat(fcorners, length(sacfiles), 1);
end

% initialize the output variables
noise.std_all = nan(length(sacfiles), 1);
noise.std_before = nan(length(sacfiles), 1);
noise.std_after = nan(length(sacfiles), 1);
noise.psd_section_raw = cell(length(sacfiles), 1);
noise.psd_section_hipass = cell(length(sacfiles), 1);
noise.psd_bandstop = cell(length(sacfiles), 1);

for ii = 1:length(sacfiles)
    % read the seismograms and its header
    [seis, hdr] = readsac(sacfiles{ii});
    [dt_ref, ~, ~, fs, ~, dts] = gethdrinfo(hdr);
    
    % determine arrival time if it does not exist
    if isnan(hdr.T0) || hdr.T0 == badval
        % first, try to get travel time if earthquake info is available
        if ~isnan(hdr.USER7) && hdr.USER7 ~= badval
            try
                ev = irisFetch.Events('eventID', string(hdr.USER7));
                tt = tauptime('mod', 'ak135', 'dep', hdr.EVDP, ...
                    'ph', 'p,P,Pdiff,PKP,PKIKP', 'deg', hdr.GCARC);
                hdr.T0 = seconds(datetime(ev(1).PreferredTime, 'Format', ...
                    'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', 'TimeZone', 'UTC') + ...
                    seconds(tt(1).time) - dt_ref);
            catch ME
                fprintf('%s\n', ME.getReport);
                keyboard
            end
        % else, pick the arrival time 
        else
            % TODO: Implement
            continue
        end
        
        if isnan(hdr.T0) || hdr.T0 == badval
            keyboard
        end
    end
    
    % remove the instruemnt response
    pa = real(counts2pa(seis, fs, [0.05 0.1 10 20], [], [], false));
    
    % time relative to the picekd first arrival time
    t = seconds(dts - dt_ref) - hdr.T0;
    
    % compute standard deviations of the noise
    noise.std_all(ii) = std(pa);
    noise.std_before(ii) = std(pa(t < cutoff(ii)));
    noise.std_after(ii) = std(pa(t >= cutoff(ii)));
    
    % compute the power spectral densities
    try
        [SD, F] = pchave({pa(t < cutoff(ii))}, 500, 70, 500, fs);
        Flog = logspace(log10(F(2)), log10(F(end)), length(F) -1);
        SDlog = interp1(F, SD, Flog);
        noise.F = Flog;
        noise.psd_section_raw{ii} = 10 * log10(SDlog);
        
        % highpass the signal like "all" in freqselect.m
        pah = detrend(pa .* shanning(length(pa), 0.05, 0), 1);
        pah = hipass(pah, fs, 0.4, 4, 2, 'butter', 'linear');
        
        [SD, F] = pchave({pah(t < cutoff(ii))}, 500, 70, 500, fs);
        Flog = logspace(log10(F(2)), log10(F(end)), length(F) -1);
        SDlog = interp1(F, SD, Flog);
        noise.psd_section_hipass{ii} = 10 * log10(SDlog);
    catch ME
        fprintf('%s\n', ME.getReport);
        noise.psd_section_raw{ii} = nan;
        noise.psd_section_hipass{ii} = nan;
        continue
    end
    try
        % bandpass the signal
        pas = bandstop(pah, fs, fcorners(ii, 1), fcorners(ii, 2), 4, 2, ...
            'butter', 'linear');
        
        [SD, F] = pchave(pas, 500, 70, 500, fs);
        Flog = logspace(log10(F(2)), log10(F(end)), length(F) -1);
        SDlog = interp1(F, SD, Flog);
        noise.psd_bandstop{ii} = 10 * log10(SDlog);
    catch ME
        fprintf('%s\n', ME.getReport);
        pah = nan .* pa;
        pas = nan .* pah;
        noise.psd_bandstop{ii} = nan;
        continue
    end
    
    %% making figure
    if plt
        figure(2)
        clf
        set(gcf, 'Units', 'inches', 'Position', [0 2 6 6])
        
        % title of the figure
        ax0 = subplot('Position', [0.15 0.93 0.70 0.01]);
        title_str = sprintf(['Event ID: %d, Station: %s, f_{corner}: ' ...
            '%.2f--%.2f Hz'], hdr.USER7, replace(hdr.KSTNM, ' ', ''), ...
            fcorners(ii, 1), fcorners(ii, 2));
        title(ax0, title_str, 'FontSize', 12)
        set(ax0, 'Color', 'none')
        set(ax0.XAxis, 'Visible', 'off')
        set(ax0.YAxis, 'Visible', 'off')
        
        % time-series data: (a) raw, (b) hipass 0.4 Hz, (c) bandstop
        ax1 = subplot('Position', [0.19 0.69 0.72 0.24]);
        maxpa = max(abs(pa));
        maxpah = max(abs(pah));
        plot(ax1, t, pa / maxpa, 'LineWidth', 0.5, 'Color', [1 1 1] * 0.75)
        hold on
        plot(ax1, t(t < cutoff(ii)), pa(t < cutoff(ii)) / maxpa, 'LineWidth', 1, ...
            'Color', 'r')
        plot(ax1, t, pah / maxpah - 2, 'LineWidth', 1, 'Color', ...
            [1 1 1] * 0.75)
        plot(ax1, t(t < cutoff(ii)), pah(t < cutoff(ii)) / maxpah - 2, ...
            'LineWidth', 1, 'Color', [0 0.5 0])
        plot(ax1, t, pas / max(abs(pas)) - 4, 'LineWidth', 1, 'Color', 'b')
        grid on
        xlabel(ax1, 'time since picked arrival (s)')
        set(ax1, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12, ...
            'XLim', [min(t) max(t)], 'YLim', [-5.2 1.2], ...
            'YTick', [-4 -2 0], 'YTickLabel', ...
            {sprintf('hp 0.4 bs %.2f-%.2f', fcorners(ii, 1), ...
            fcorners(ii, 2)), 'hp 0.4', 'raw'})
        set(ax1.YAxis, 'TickLabelRotation', 45)
        
        % power spectral densities
        ax2 = subplot('Position', [0.19 0.11 0.72 0.36]);
        semilogx(ax2, noise.F, noise.psd_section_raw{ii}, ...
            'LineWidth', 1, 'Color', 'r')
        hold on
        semilogx(ax2, noise.F, noise.psd_section_hipass{ii}, ...
            'LineWidth', 1, 'Color', [0 0.5 0])
        semilogx(ax2, noise.F, noise.psd_bandstop{ii}, ...
            'LineWidth', 1, 'Color', 'b')
        grid on
        xlabel(ax2, 'frequency (Hz)')
        ylabel(ax2, '10 log_{10} spectral density (Pa^2/Hz)')
        set(ax2, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)
        legend('raw', 'hp 0.4', sprintf('hp 0.4 bs %.2f-%.2f', ...
            fcorners(ii, 1), fcorners(ii, 2)), 'Location', 'southwest')
        
        % add period as the second x-axis labels
        ax2d = doubleaxes(ax2);
        inverseaxis(ax2d.XAxis, 'period (s)');
        
        % save the figure
        set(gcf, 'Renderer', 'painters')
        savename = sprintf('%s_%d_%s', mfilename, hdr.USER7, ...
            replace(hdr.KSTNM, ' ', ''));
        figdisp(savename, [], [], 2, [], 'epstopdf')
    end
end
end