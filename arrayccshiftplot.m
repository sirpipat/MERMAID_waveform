function arrayccshiftplot(CCmaxs, t_shifts, metadata)
% arrayccshiftplot(CCmaxs, t_shifts, metadata)
%
% Generates plots of source-receivers with MERMAID station numbers, maximum
% correlation coefficients, and best time-shifts for each number.
%
% INPUT:
% t_shifts          Best time shift where CC is maximum
% CCmaxs            Maximum correlation coefficient
% metadata          SAC header variables sorted by variable names
%
% SEE ALSO:
% COMPAREPRESSURE_ROUTINE
%
% Last modified by sirawich-at-princeton.edu, 03/01/2022

% get station number
metadata.STNM = zeros(size(metadata.T0));
for ii = 1:length(metadata.STNM)
    metadata.STNM(ii) = str2double(indeks(metadata.KSTNM{ii},4:5));
end

[uniqevent, ~, ~] = unique(metadata.USER7);
for ii = 1:length(uniqevent)
    whevent = (metadata.USER7 == uniqevent(ii));
    if sum(whevent) >= 2
        % list of source and receriver locations for an event
        stlo = mod(metadata.STLO(whevent), 360);
        stla = metadata.STLA(whevent);
        evlo = indeks(unique(mod(metadata.EVLO(whevent), 360)), 1);
        evla = indeks(unique(metadata.EVLA(whevent)), 1);
        
        % MERMAID number
        n = length(stlo);
        
        % map extent
        latmin = min(min(stla), evla);
        latmax = max(max(stla), evla);
        lonmin = min(min(stlo), evlo);
        lonmax = max(max(stlo), evlo);
        
        % extend the max extent by 10%
        latmid = (latmin + latmax) / 2;
        halfheight = (latmax - latmin) / 2;
        lonmid = (lonmin + lonmax) / 2;
        halfwidth = (lonmax - lonmin) / 2;
        
        latmin = latmid - 1.1 * halfheight;
        latmax = latmid + 1.1 * halfheight;
        lonmin = lonmid - 1.1 * halfwidth;
        lonmax = lonmid + 1.1 * halfwidth;
        
        figure(3)
        clf
        set(gcf, 'Units', 'inches', 'Position', [18 8 6 8]);
        ax1 = subplot('Position', subplotposition(3, 1, 1, ...
            [0.08 0.12 0.02 0.2], [0.05 0.05 0.01 -0.04]));
        scatter(mod(metadata.STLO(whevent), 360), ...
            metadata.STLA(whevent), 50, (1:n)', ...
            'filled', 'Marker', 'v', 'MarkerEdgeColor', 'k');
        grid on
        hold on
        box on
        % TODO: plotcont plotplates (doubleaxes, addfocalmech)
        [~, cont] = plotcont();
        plate = plotplates();
        plate.Color = 'r';
        
        % zoom in the map
        original_x2y_ratio = (ax1.XLim(2)-ax1.XLim(1))/(ax1.YLim(2)-ax1.YLim(1));
        new_x2y_ratio = (lonmax-lonmin)/(latmax-latmin);
        if new_x2y_ratio > original_x2y_ratio
            latmid = (latmin + latmax) / 2;
            latmin = latmid - (lonmax - lonmin) / original_x2y_ratio / 2;
            latmax = latmid + (lonmax - lonmin) / original_x2y_ratio / 2;
        else
            lonmid = (lonmin + lonmax) / 2;
            lonmin = lonmid - (latmax - latmin) * original_x2y_ratio / 2;
            lonmax = lonmid + (latmax - latmin) * original_x2y_ratio / 2;
        end
        
        
        maxSTNM = max(metadata.STNM(whevent));
        minSTNM = min(metadata.STNM(whevent));
        colormap(gca, jet(n));
        colorbar('Ticks', 1:n, 'TickLabels', metadata.STNM(whevent))
        xlim([lonmin lonmax])
        ylim([latmin latmax])
        caxis([0.5 n+0.5])
        xlabel('longitude (degrees)')
        ylabel('latitude (degrees)')
        set(gca, 'TickDir', 'both', 'FontSize', 12)
        
        % add event's focal mechanism
        ax1s = doubleaxes(ax1);
        axes(ax1s);
        ax1s.XAxisLocation = 'bottom';
        ax1s.PlotBoxAspectRatio = ax1.PlotBoxAspectRatio;
        addfocalmech(ax1s, [evlo evla], 'PublicID', string(uniqevent(ii)));
        ax1s.XLim = ax1.XLim;
        ax1s.YLim = ax1.YLim;
        ax1s.Visible = 'off';
        ax1s.XAxis.Visible = 'off';
        ax1s.YAxis.Visible = 'off';
        ax1s.TickDir = 'both';
        title(ax1, sprintf('Event ID: %d, MERMAID number', uniqevent(ii)))
        
        ax2 = subplot('Position', subplotposition(3, 1, 2, ...
            [0.08 0.12 0.02 0.2], [0.05 0.05 0.01 -0.04]));
        scatter(mod(metadata.STLO(whevent), 360), ...
            metadata.STLA(whevent), 50, CCmaxs(whevent,1), ...
            'filled', 'Marker', 'v', 'MarkerEdgeColor', 'k');
        grid on
        hold on
        box on
        [~, cont] = plotcont();
        plate = plotplates();
        plate.Color = 'r';
        
        colormap(gca, bwrmap(255, 'rwb'));
        colorbar
        caxis([0 1])
        xlim([lonmin lonmax])
        ylim([latmin latmax])
        xlabel('longitude (degrees)')
        ylabel('latitude (degrees)')
        set(gca, 'TickDir', 'both', 'FontSize', 12)
        
        % add event's focal mechanism
        ax2s = doubleaxes(ax2);
        axes(ax2s);
        ax2s.XAxisLocation = 'bottom';
        ax2s.PlotBoxAspectRatio = ax2.PlotBoxAspectRatio;
        addfocalmech(ax2s, [evlo evla], 'PublicID', string(uniqevent(ii)));
        ax2s.XLim = ax2.XLim;
        ax2s.YLim = ax2.YLim;
        ax2s.Visible = 'off';
        ax2s.XAxis.Visible = 'off';
        ax2s.YAxis.Visible = 'off';
        ax2s.TickDir = 'both';
        title(ax2, sprintf('Event ID: %d, maximum correlation coefficients', uniqevent(ii)))
        
        ax3 = subplot('Position', subplotposition(3, 1, 3, ...
            [0.08 0.12 0.02 0.2], [0.05 0.05 0.01 -0.04]));
        scatter(mod(metadata.STLO(whevent), 360), ...
            metadata.STLA(whevent), 50, t_shifts(whevent, 1), ...
            'filled', 'Marker', 'v', 'MarkerEdgeColor', 'k');
        grid on
        hold on
        box on
        [~, cont] = plotcont();
        plate = plotplates();
        plate.Color = 'r';
        
        colormap(gca, bwrmap(255, 'rwb'));
        colorbar
        xlim([lonmin lonmax])
        ylim([latmin latmax])
        xlabel('longitude (degrees)')
        ylabel('latitude (degrees)')
        set(gca, 'TickDir', 'both', 'FontSize', 12)
        
        % add event's focal mechanism
        ax3s = doubleaxes(ax3);
        axes(ax3s);
        ax3s.XAxisLocation = 'bottom';
        ax3s.PlotBoxAspectRatio = ax3.PlotBoxAspectRatio;
        addfocalmech(ax3s, [evlo evla], 'PublicID', string(uniqevent(ii)));
        ax3s.XLim = ax3.XLim;
        ax3s.YLim = ax3.YLim;
        ax3s.Visible = 'off';
        ax3s.XAxis.Visible = 'off';
        ax3s.YAxis.Visible = 'off';
        ax3s.TickDir = 'both';
        title(ax3, sprintf('Event ID: %d, best time-shift (s)', uniqevent(ii)))
        
        [azDeg,distDeg] = azimproj([evlo evla], [stlo stla]);
        
        figure(2)
        clf
        set(gcf, 'Units', 'inches', 'Position', [12 8 6 8]);
        ax4 = polaraxes('ThetaZeroLocation', 'top', 'ThetaDir', ...
            'clockwise', 'Position', [0.0875 0.5444 0.8265 0.3966]);
        hold on
        polarscatter(ax4, azDeg*pi/180, distDeg, 50, (1:n)', 'filled', ...
            'Marker', 'v', 'MarkerEdgeColor', 'k')
        colormap(ax4, jet(n))
        c4 = colorbar('Ticks', 1:n, 'TickLabels', ...
            metadata.STNM(whevent), 'Position', ...
            [0.89 0.5442 0.0370 0.3968]);
        caxis([0.5 n+0.5])
        rlim([0 30*ceil(max(distDeg)/30)])
        rticks(0:30:180)
        
        ax5 = subplot('Position', subplotposition(2, 1, 2, ...
            [0.08 0.12 0.02 0.2], [0.05 0.05 0.01 -0.04]));
        scatter(ax5, azDeg, distDeg, 50, (1:n)', 'filled', ...
            'Marker', 'v', 'MarkerEdgeColor', 'k')
        grid on
        box on
        colormap(ax5, jet(n))
        colorbar('Ticks', 1:n, 'TickLabels', metadata.STNM(whevent))
        caxis([0.5 n+0.5])
        
        keyboard
    end
end
end