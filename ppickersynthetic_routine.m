function ppickersynthetic_routine(sacfiles, option, method)
% ppickersynthetic_routine(sacfiles, option)
%
% A script to test how good a P-phase picker is on some synthetic seismic traces.
%
% INPUT
% sacfiles      cell array containing sacfiles
% option        1 -- without filtering
%               2 -- bandpass 1--2 Hz
% method        method for the pickers
%               1 -- PphasePicker
%               2 -- simplepicker : signal value threshold
%               3 -- simplepicker : optimal signal-to-noise ratio
%
% PphasePicker see:
% Kalkan, E. (2016). "An automatic P-phase arrival time picker", Bull. of
%   Seismol. Soc. of Am., 106, No. 3, doi: 10.1785/0120150111
%
% simplepicker
% uses STA/LTA algorithm to find the P-wave arrival window. Then, from the
% window, determine either the first instance where the value exceeds the
% threshold (5% of the max of the window) or the instance at which the
% signal (any data points after the instance) to noise (any data points
% before the instance) is the greatest.
%
% Last modified by sirawich-at-princeton.edu, 01/13/2022

% option 1: picking from raw seismogram
if option == 1
    for ii = 1:length(sacfiles)
        % read the sac file
        [SeisData,HdrData,~,~,tims] = readsac(sacfiles{ii});
        dt_ref = gethdrinfo(HdrData);
        
        % pick the P-wave arrival
        if method == 1
            [loc,snr_db] = PphasePicker(SeisData, HdrData.DELTA, 'sm', ...
                'n', 5, 0.005, 8, 'full');
        else
            loc = simplepicker(SeisData, HdrData.DELTA, tims, method);
        end
        
        % make a figure
        figure
        set(gcf, 'Units', 'inches', 'Position', [18 12 7 6.6])
        
        % plot the seismogram within 5 seconds from the picked arrival
        subplot(3,1,1)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(tims-loc, SeisData, 'LineWidth', 1, 'Color', 'k')
        xlim([-5 5])
        fprintf('loc = %f\n', loc);
        ylim([-1.1 1.1] * max(abs(SeisData(and(tims > loc(1) - 5, ...
            tims < loc(1) + 5)))))
        vline(gca,0,'Color',[1 0 0],'LineWidth',1);
        grid on
        xlabel('time since the pick (s)')
        title(sprintf('Origin: %s, ID: %d, Mw = %5.2f', ...
            string(dt_ref), HdrData.USER7, HdrData.MAG));

        % plot the seismogram within 60 seconds from the picked arrival
        subplot(3,1,2)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(tims-loc, SeisData, 'LineWidth', 1, 'Color', 'k')
        xlim([-60 60])
        ylim([-1.1 1.1] * max(abs(SeisData(and(tims > loc(1)-60, ...
            tims < loc(1) + 60)))))
        vline(gca,0,'Color',[1 0 0],'LineWidth',1);
        grid on
        xlabel('time since the pick (s)')

        % plot the whole seismogram
        subplot(3,1,3)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(tims,SeisData,'LineWidth',1,'Color','k')
        xlim([tims(1) tims(end)])
        ylim([-1.1 1.1]*max(abs(SeisData)))
        vline(gca,loc(1),'Color',[1 0 0],'LineWidth',1);
        grid on
        xlabel('time (s)')
        
        % save the figure
        if method == 1
            method_phrase = 'PphasePicker';
        elseif method == 2
            method_phrase = 'simplepicker_threshold';
        else
            method_phrase = 'simplepicker_snr';
        end
        savename = sprintf('%s_raw_%s_%d_%s.eps', mfilename, ...
            method_phrase, HdrData.USER7, replace(HdrData.KSTNM, ' ', ''));
        figdisp(savename, [], [], 2, [], 'epstopdf');
        delete(gcf)
    end
% option 2: bandpass 1--2 Hz before picking
else
    for ii = 1:length(sacfiles)
        % read the sac file
        [SeisData,HdrData,~,~,tims] = readsac(sacfiles{ii});
        dt_ref = gethdrinfo(HdrData);
        
        % resampled to 20.00703 Hz
        fs = 20.00703;
        t = tims(1):1/fs:tims(end);
        x = shannon(tims, SeisData,t);
        
        % filter 1--2 Hz
        x = bandpass(x, fs, 1, 2, 2, 2, 'butter', 'linear');
        
        % hanning window
        w = shanning(length(t), 0.01, 0);
        x = x .* w;
        
        % pick the P-wave arrival
        if method == 1
            [loc,snr_db] = PphasePicker(x, 1/fs, 'sm', 'n', 5, 0.005, ...
                20, 'full');
        else
            loc = simplepicker(x, 1/20.00703, t, method, ...
                0.01 * (t(end) - t(1)));
        end
        
        % make a figure
        figure
        set(gcf, 'Units', 'inches', 'Position', [18 12 7 8.8])
        
        % plot the seismogram within 5 seconds from the picked arrival
        subplot(4,1,1)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(t-loc, x, 'LineWidth', 1, 'Color', 'k')
        xlim([-5 5])  
        ylim([-1.1 1.1]*max(abs(x(and(t>loc(1)-5,t<loc(1)+5)))))      
        vline(gca,0,'Color',[1 0 0],'LineWidth',1);                 
        grid on
        xlabel('time since the pick (s)')
        %nolabels(gca, 2)
        title(sprintf('Origin: %s, ID: %d, Mw = %5.2f', ...
            string(dt_ref), HdrData.USER7, HdrData.MAG));
        
        % plot the seismogram within 60 seconds from the picked arrival
        subplot(4,1,2)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(t-loc, x, 'LineWidth', 1, 'Color', 'k')
        xlim([-60 60])
        ylim([-1.1 1.1]*max(abs(x(and(t>loc(1)-60,t<loc(1)+60)))))
        vline(gca,0,'Color',[1 0 0],'LineWidth',1);
        grid on
        xlabel('time since the pick (s)')
        
        % plot the whole seismogram
        subplot(4,1,3)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(t,x,'LineWidth',1,'Color','k')
        xlim([t(1) t(end)])
        ylim([-1.1 1.1]*max(abs(x)))
        vline(gca,loc(1),'Color',[1 0 0],'LineWidth',1);
        grid on
        xlabel('time (s)')
        
        % plot the whole raw seismogram
        subplot(4,1,4)
        set(gca, 'FontSize', 12, 'TickDir', 'both')
        plot(tims,SeisData,'LineWidth',1,'Color','k')
        xlim([tims(1) tims(end)])
        ylim([-1.1 1.1]*max(abs(SeisData)))
        vline(gca,loc(1),'Color',[1 0 0],'LineWidth',1);
        xlabel('time since the pick (s)')
        grid on

        % save the figure
        if method == 1
            method_phrase = 'PphasePicker';
        elseif method == 2
            method_phrase = 'simplepicker_threshold';
        else
            method_phrase = 'simplepicker_snr';
        end
        savename = sprintf('%s_filtered_%s_%d_%s.eps', mfilename, ...
            method_phrase, HdrData.USER7, replace(HdrData.KSTNM, ' ', ''));
        figdisp(savename, [], [], 2, [], 'epstopdf');
        delete(gcf)
    end
end

end

% Picks the first P-wave arrival. It uses STA/LTA algorithm to find the
% arrival window then picks based on either the signal goes above the
% threshold or the optial signal-to-noise ratio.
%
% INPUT:
% x         signal
% dt        sampling interval
% t         time for each sample
% method    picking method:
%           2 -- simplepicker : signal value threshold
%           3 -- simplepicker : optimal signal-to-noise ratio
% mintime   when in the signal at which the picker start to look at
%           [default: 0]
%
% OUPUT:
% loc       the picked arrival time
function loc = simplepicker(x, dt, t, method, mintime)
defval('mintime', 0)

% find the window using STA/LTA algorithm
trigt = stalta(x, dt, [t(1) t(end)], 10, 30, 1.2, 1.2, 10, 5, 20, 5);

% slice for a window
i_window = 1;
wh = and(t > trigt(1,1), t < trigt(1,2) - 10);
x_cut = x(wh);
t_cut = t(wh);

% method 2 -- simplepicker : signal value threshold
if method == 2
    % check whether the window contains a peak
    while (mean(x_cut) == 0 || std(x_cut) == 0) && i_window < ...
            size(trigt, 1)
        i_window = i_window + 1;
        wh = and(t > trigt(i_window,1), t < trigt(i_window,2) - 10);
        x_cut = x(wh);
        t_cut = t(wh);
    end
    
    % identify the rise of the peak
    loc = indeks(t_cut(abs(x_cut) > 2e-2 * max(abs(x_cut))), 1);
% method 3 -- simplepicker : optimal signal-to-noise ratio
else
    tt = t_cut(2:end-1)';
    snr = zeros(size(tt));

    for ii = 1:length(tt)
        sig = mean(x_cut(t_cut > tt(ii)) .^ 2);
        noi = mean(x_cut(t_cut < tt(ii)) .^ 2);
        snr(ii) = 10 * log10(sig / (noi + 1e-4));
    end

    loc = tt(snr == max(snr));
end

% reject the pick if the pick time is not unique and look at the next
% window
i_window = 1;
while (length(loc) > 1 || t_cut(1) < mintime) && i_window < size(trigt, 1)
    i_window = i_window + 1;
    wh = and(t > trigt(i_window,1), t < trigt(i_window,2) - 10);
    x_cut = x(wh);
    t_cut = t(wh);
    
    % method 2 -- simplepicker : signal value threshold
    if method == 2
        while (mean(x_cut) == 0 || std(x_cut) == 0) && i_window < ...
                size(trigt, 1)
            i_window = i_window + 1;
            wh = and(t > trigt(i_window1,1), t < trigt(i_window,2) - 10);
            x_cut = x(wh);
            t_cut = t(wh);
        end
        loc = indeks(t_cut(abs(x_cut) > 2e-2 * max(abs(x_cut))), 1);
    % method 3 -- simplepicker : optimal signal-to-noise ratio
    else
        tt = t_cut(2:end-1)';
        snr = zeros(size(tt));

        for ii = 1:length(tt)
            sig = mean(x_cut(t_cut > tt(ii)) .^ 2);
            noi = mean(x_cut(t_cut < tt(ii)) .^ 2);
            snr(ii) = 10 * log10(sig / (noi + 1e-4));
        end

        loc = tt(snr == max(snr));
    end
end
end
