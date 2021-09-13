function [ax1, ax2] = plotrecords(allfiles, op1, op2, op3, op4, op5)
% [ax1, ax2] = PLOTRECORDS(allfiles, op1)
%
% Plots the seismic traces of all stations (sorted by names?) and the map
% containing all stations and the source. All first P-wave arrivals of 
% seismic traces are aligned.
%
% INPUT:
% allfiles      full path to all SAC files
% op1           options for time axis
%               1  --   relative to picked first P arrival
%               2  --   relative to model ak135
%               3  --   absolute time
%
% OUTPUT:
% ax1           axes handle to the seismic trace plot
% ax2           axes handle to the map
%
% Last modified by sirawich-at-princeton.edu, 09/13/2021

%% create figures
fig1 = figure(1);
set(gcf, 'Units', 'inches', 'Position', [2 2 8 10])
clf;
% plot seismic from all stations
ax1 = subplot(1,1,1, 'Box', 'on', 'TickDir', 'both');
cla
hold on

% plot stations on the map with the source
fig2 = figure(2);
set(gcf, 'Units', 'inches', 'Position', [2 12 8 5])
clf
ax2 = subplot(1,1,1, 'Box', 'on', 'TickDir', 'both');
cla
hold on
grid on

%% plot the coastlines and plate boundaries on the map
% plot tcoastlines
[axlim,handl,XYZ] = plotcont([0 90], [360 -90], 1, 0);
% plot plate boundaries
[handlp, XYp] = plotplates([0 90], [360 -90], 1);
handlp.Color = 'r';

%% plot all seismograms and source, paths, and stations on the map
% limits of the time window of the seismograms
window_left = -10;
window_right = 5;
for ii = 1:length(allfiles)
    [SeisData, HdrData, ~, ~, tims]=readsac(allfiles{ii});
    [dt_ref, dt_B, dt_E, fs, npts, dts, ~] = gethdrinfo(HdrData);
    eqid = HdrData.USER7;
    if ii == 1
        % set the overall y-limit
        dist_limit = [HdrData.GCARC HdrData.GCARC];
    else
        % set the overall y-limit
        if HdrData.GCARC < dist_limit(1)
            dist_limit(1) = HdrData.GCARC;
        end
        if HdrData.GCARC > dist_limit(2)
            dist_limit(2) = HdrData.GCARC;
        end
    end
    if isnan(HdrData.T0) || HdrData.T0 == -12345
        continue
    end
    % convert digital counts to pressure
    x = real(counts2pa(SeisData, fs));
    % filter out the ambient noise below 1 Hz
    x = bandpass(x, fs, 1, 2, 2, 1, 'butter', 'linear');
    % figure out time axis given the option
    switch op1
        case 1
            wh = and(tims >= HdrData.T0 + window_left, ...
                tims <= HdrData.T0 + window_right);
            x = x(wh);
            t = tims(wh) - HdrData.T0;
            idx = length(t) + 1 - sum(t >= - HdrData.USER4);
        case 2
            wh = and(tims >= HdrData.T0 - HdrData.USER4 + window_left, ...
                tims <= HdrData.T0 - HdrData.USER4 + window_right);
            x = x(wh);
            t = tims(wh) - (HdrData.T0 - HdrData.USER4);
            idx = length(t) + 1 - sum(t >= HdrData.USER4);
        otherwise
            t = seconds(tims) + dt_ref;
            idx = -1;
            if ii == 1
                time_limit = [t(1) t(end)];
            else
                if t(1) < time_limit(1)
                    time_limit(1) = t(1);
                end
                if t(end) > time_limit(2)
                    time_limit(2) = t(end);
                end
            end
    end
    x = x / max(abs(x));
    if and(idx > 1, idx < length(t))
        signalplot(x(1:idx) + HdrData.GCARC, fs, t(1), ax1, [], [], ...
            rgbcolor(string(mod(ii-1,7)+1)), ...
            'LineWidth', 1);
        signalplot(x(idx:end) + HdrData.GCARC, fs, t(idx), ax1, [], [], ...
            rgbcolor(string(mod(ii-1,7)+1)), ...
            'LineWidth', 2);
    else
        signalplot(x + HdrData.GCARC, fs, t(1), ax1, [], [], ...
            rgbcolor(string(mod(ii-1,7)+1)), ...
            'LineWidth', 1);
    end
    text(ax1, t(40), HdrData.GCARC + 0.7, HdrData.KSTNM, ...
        'Color', rgbcolor(string(mod(ii-1,7)+1)), ...
        'FontSize', 12);
    % plot the path from source to station
    plottrack(ax2, [HdrData.EVLO HdrData.EVLA], [HdrData.STLO HdrData.STLA], 0, ...
          100, 'LineWidth', 0.5, 'Color', [0 0.5 0.9]);
    % plot the station
    scatter(ax2, mod(HdrData.STLO,360), HdrData.STLA, 20, ...
        'Marker', 'v', ...
        'MarkerEdgeColor', rgbcolor('k'), ...
        'MarkerEdgeAlpha', 0.7, ...
        'MarkerFaceColor', rgbcolor(string(mod(ii-1,7)+1)), ...
        'MarkerFaceAlpha', 0.7);
end
% plot the source
addfocalmech(ax2, 'PublicID', sprintf('%d', HdrData.USER7));

%% seismogram plot decoration
if op1 == 1
    ax1.XLim = [window_left window_right];
    ax1.XLabel.String = 'time relative to picked first P-wave arrival (s)';
elseif op1 == 2
    ax1.XLim = [window_left window_right];
    ax1.XLabel.String = 'time relative to predicted first P-wave arrival (s)';    
else
    ax1.XLim = time_limit;
    ax1.XLabel.String = 'absolute time (hh:mm)';
end
ax1.YLim = dist_limit + [-1 2];
ax1.YLabel.String = 'distance (degrees)';
ax1.Title.String = sprintf('Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', eqid, HdrData.MAG, HdrData.EVDP);
ax1.FontSize = 12;

%% map decoration
% ticks label
ax2.XTick = 0:30:360;
ax2.XTickLabel = {'0', '30', '60', '90', '120', '150', '180', '-150', ...
                  '-120', '-90', '-60', '-30', '0'};
ax2.YTick = -90:15:90;
ax2.YTickLabel = {'-90', '-75', '-60', '-45', '-30', '-15', '0', '15', ...
                  '30', '45', '60', '75', '90'};
ax2.XLim = [0 360];
ax2.YLim = [-90 90];
ax2.Title.String = sprintf('Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', eqid, HdrData.MAG, HdrData.EVDP);
ax2.FontSize = 12;

%% save figures
figure(fig1);
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_%d_arrivals.eps', mfilename, eqid), [], [], 2, [], 'epstopdf')
figure(fig2);
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_%d_map.eps', mfilename, eqid), [], [], 2, [], 'epstopdf')
end