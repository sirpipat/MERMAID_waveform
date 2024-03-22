function plotmermaid(obsmasterdir, synmasterdir, fcorners, CCmaxs, ...
    t_shifts, metadata, op1, op2, op3, op4, op5)
% PLOTMERMAID(obsmasterdir, synmasterdir, fcorners, CCmaxs, ...
%     t_shifts, metadata, op1, op2, op3, op4, op5)
%
% Plots a source-receiver map with the focal mechanism and the radiation
% pattern as well as the pressure records from MERMAIDs
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% fcorners          corner frequencies used for comparing synthetic and
%                   observed acoustic pressures
% CCmaxs            Maximum correlation coefficient
% t_shifts          Best time shift where CC is maximum
% metadata          SAC header variables sorted by variable names
% op1               options for y-axis sorting [Default: 2]
%                   1  --   by distance (degrees)
%                   2  --   by azimuth
% op2               options for y-axis scaling and positioning [Default: 2]
%                   1  --   actual value from op1
%                   2  --   equal spreading for readability
% op3               options for zoom-in verions [Default: 2]
%                   1  --   no zoom in
%                   2  --   zoom in
% op4               options for filter [Default: 2]
%                   1  --   no filter
%                   2  --   filter with chosen corner frequencies
% op5               options for graying out traces if CCmaxs <= 0.6
%                   1  --   no graying out
%                   2  --   graying out
%
% SEE ALSO:
% PLOTINSTASEIS, PLOTSYNTHETICS, PLOTRECORDS, ARRAYCCSHIFTPLOT
%
% Last modified by sirawich-at-princeton.edu, 03/20/2024

defval('op1', 2)
defval('op2', 2)
defval('op3', 2)
defval('op4', 2)
defval('op5', 2)

tsmul = 0;

if op3 == 1
    window_plot = [-40 60];
else
    window_plot = [-20 40];
end

%% compute for extra metadata
% get station number
metadata.STNM = zeros(size(metadata.T0));
for ii = 1:length(metadata.STNM)
    % remove whitespace first
    metadata.KSTNM{ii} = indeks(metadata.KSTNM{ii},1:5);
    metadata.STNM(ii) = str2double(indeks(metadata.KSTNM{ii},4:5));
end

% compute relative travel-time difference
metadata.DLNT = t_shifts ./ (metadata.T0 - metadata.USER8);

[uniqevent, ~, ~] = unique(metadata.USER7);

% iterate over all events
for ii = 1:length(uniqevent)
    whevent = (metadata.USER7 == uniqevent(ii));
    
    if uniqevent(ii) == 11618973
        continue
    end
    
    % only make a plot when there are more than one station
    if sum(whevent) >= 2
        %% list used metadata for making map
        % list of source and receriver locations for an event
        stlo = mod(metadata.STLO(whevent), 360);
        stla = metadata.STLA(whevent);
        evlo = indeks(unique(mod(metadata.EVLO(whevent), 360)), 1);
        evla = indeks(unique(metadata.EVLA(whevent)), 1);
        evdp = indeks(unique(metadata.EVDP(whevent)), 1);
        
        % MERMAID number
        n = length(stlo);
        
        %% plot the map of source and receivers
        % map extent
        latmin = min(min(stla), evla);
        latmax = max(max(stla), evla);
        lonmin = min(min(stlo), evlo);
        lonmax = max(max(stlo), evlo);
        
        % extend the max extent by 20%
        latmid = (latmin + latmax) / 2;
        halfheight = (latmax - latmin) / 2;
        lonmid = (lonmin + lonmax) / 2;
        halfwidth = (lonmax - lonmin) / 2;
        
        latmin = latmid - 1.2 * halfheight;
        latmax = latmid + 1.2 * halfheight;
        lonmin = lonmid - 1.2 * halfwidth;
        lonmax = lonmid + 1.2 * halfwidth;
        
        fig = figure(3);
        clf
        set(gcf, 'Units', 'inches', 'Position', [0 2 9 9])
        
        ax1 = subplot('Position', [0.08 0.72 0.89 0.26]);
        scatter(mod(metadata.STLO(whevent), 360), ...
            metadata.STLA(whevent), 80, (1:n)', ...
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
        cmap = jet(n);
        colormap(gca, cmap);
        c = colorbar('Ticks', 1:n, 'TickLabels', metadata.STNM(whevent));
        c.Label.String = 'MERMAID number';
        c.Label.FontSize = 12;
        xlim([lonmin lonmax])
        ylim([latmin latmax])
        caxis([0.5 n+0.5])
        xlabel('longitude (degrees)')
        ylabel('latitude (degrees)')
        
        % add countour lines
        addcontourlines(ax1, evlo, evla);
        
        set(gca, 'TickDir', 'out', 'FontSize', 12, 'GridAlpha', 0.5, ...
            'GridLineStyle', ':')
        
        %% add P-wave polarity
        % plot lon/lat grid
        resolution_multiplier = 5;
        nptx = 14 * resolution_multiplier + 1;
        npty = 4 * resolution_multiplier + 1;
        pllo = linspace(ax1.XLim(1), ax1.XLim(2), nptx)';
        plla = linspace(ax1.YLim(1), ax1.YLim(2), npty);
        [pllo, plla] = meshgrid(pllo, plla);
        
        % focal mechanism
        [~, ~, CMT] = getfocalmech('PublicID', string(uniqevent(ii)));
        if ~isempty(CMT)
            M = CMT.M;
        
            % polarity
            fp = cmtpolarity(M, evla, evlo, evdp, plla, pllo, 'ak135');

            % colormap for polarity plot
            c0 = [1 1 1];
            cn = [1 0.83 0.88];
            icmap = (0:2)' / 2;
            cmap_polarity = icmap * cn + (1 - icmap) * c0;

            ax1p = doubleaxes(ax1);
            axes(ax1p)
            imagesc(ax1p,[floor(ax1.XLim(1)) ceil(ax1.XLim(2))], ...
                [floor(ax1.YLim(1)) ceil(ax1.YLim(2))],fp);
            axis xy
            grid on
            set(ax1p, 'CLim', [-1 1], 'TickDir', 'out', 'Box', 'on', ...
                'FontSize', 11)
            colormap(ax1p, cmap_polarity)
            ax1p.XAxis.Visible = 'off';
            ax1p.YAxis.Visible = 'off';
            set(ax1p, 'XLim', ax1.XLim, 'YLim', ax1.YLim, ...
                'Position', ax1.Position, 'DataAspectRatio', [1 1 1])
            axes(ax1)
            set(ax1, 'Color', 'None')
        end
        
        %% add event's focal mechanism
        ax1s = doubleaxes(ax1);
        axes(ax1s);
        ax1s.XAxisLocation = 'bottom';
        ax1s.PlotBoxAspectRatio = ax1.PlotBoxAspectRatio;
        addfocalmech(ax1s, [evlo evla], 'PublicID', string(uniqevent(ii)), 20);
        ax1s.XLim = ax1.XLim;
        ax1s.YLim = ax1.YLim;
        ax1s.Visible = 'off';
        ax1s.XAxis.Visible = 'off';
        ax1s.YAxis.Visible = 'off';
        ax1s.TickDir = 'both';
        
        %% added plot label
        ax1b = boxedlabel(ax1s, 'northwest', 0.03, 'norm', 'a', ...
            'FontSize', 12);
        axes(ax1b)
        ax1b.Position = [0.08 0.941 0.03 0.03];
        
        %% list metadata for plotting traces
        % filter out the data from other events
        stationid = metadata.STNM(whevent);
        stationname = metadata.KSTNM(whevent);
        CCmax = CCmaxs(whevent);
        t_shift = t_shifts(whevent);
        fcs = fcorners(whevent, :);
        gcarc = metadata.GCARC(whevent);
        azim = metadata.AZ(whevent);
        dlnt = metadata.DLNT(whevent);
        
        % pressure amplitude
        amp = zeros(size(stationid));
        
        if op1 == 1
            % sort everything by epicentral distance
            [gcarc, i_gcarc] = sort(gcarc);
            stationid = stationid(i_gcarc);
            stationname = stationname(i_gcarc);
            CCmax = CCmax(i_gcarc);
            t_shift = t_shift(i_gcarc);
            fcs = fcs(i_gcarc, :);
            azim = azim(i_gcarc);
            dlnt = dlnt(i_gcarc);
        else
            % sort everything by azimuth
            [azim, i_azim] = sort(azim);
            gcarc = gcarc(i_azim);
            stationid = stationid(i_azim);
            stationname = stationname(i_azim);
            CCmax = CCmax(i_azim);
            t_shift = t_shift(i_azim);
            fcs = fcs(i_azim, :);
            dlnt = dlnt(i_azim);
        end
        
        %% Plot the seimograms comparisons
        % the second plot consists of layers of axes [from front to back]
        % - plot label [b]
        %   [assigned property value: 4]
        % - label boxes sorted by CCMaxs from largest to smallest
        %   [assigned priority value: 2 + CCMaxs]
        % - seismograms of predicted and observed with ticks for round trip
        %   [assigned priority value: CCMaxs]
        % - waveform window highlights
        %   [assinged priority value: -1]
        % - epicentral distances y-label
        %   [assigned priority value: -2]
        % - azimuth y-label
        %   [assigned priority value: -3]
        %   
        % Then the layers are sorted by priority from largest to smallest
        
        % base layer axes
        ax2 = subplot('Position', [0.08 0.08 0.84 0.56]);
        priority_values = -1;
        axes_collection = ax2;
        
        % determine the y-limit
        if op2 == 1
            if op1 == 1
                ymid = (gcarc(end) + gcarc(1)) / 2;
                ywidth = (gcarc(end) - gcarc(1));
                ylimit = ymid + 1.12 * ywidth/2 * [-1 1];
            else
                ymid = (azim(end) + azim(1)) / 2;
                ywidth = (azim(end) - azim(1));
                ylimit = ymid + 1.12 * ywidth/2 * [-1 1];
            end
        else
            ymid = (1 + length(gcarc)) / 2;
            ywidth = (length(gcarc) - 1);
            ylimit = ymid + 1.12 * ywidth/2 * [-1 1];
        end
        
        xlabel('time since first picked arrival (s)');
        ylabel('epicentral distance (degrees)');
        
        % adjust y-axis for description labels
        set(ax2, 'Box', 'on', 'TickDir', 'out', 'XLim', window_plot, ...
            'YLim', ylimit, 'FontSize', 12);
        
        %% added plot label
        ax2b = boxedlabel(ax2, 'northwest', 0.03, 'norm', 'b', ...
            'FontSize', 12);
        ax2b.Position = [0.08 0.64 0.03 0.03];
        priority_values = [priority_values; 4];
        axes_collection = [axes_collection; ax2b];
        
        % positional parameters description labels
        is_label_left = true;
        prev_label_top = 0;
        
        for jj = 1:sum(whevent)
            % read the observed seismogram
            try
                obsfile = cindeks(ls2cell(sprintf('%s%d/*.%02d_*.sac', ...
                    obsmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
                [seis_o, hdr_o, ~, ~, tims_o] = readsac(obsfile);
            catch ME
                obsfile = cindeks(ls2cell(sprintf('%s%d/*.%04d_*.sac', ...
                    obsmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
                [seis_o, hdr_o, ~, ~, tims_o] = readsac(obsfile);
            end
            [~, ~, ~, fs_o] = gethdrinfo(hdr_o);
            tims_o = tims_o - hdr_o.T0;
            
            % remove instrument response
            pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], ...
                'sacpz', false);
            pres_o = real(pres_o);
            
            % read the synthetic vertical dispalcement at the ocean bottom
            try
                synfile = cindeks(ls2cell(sprintf('%s%d/*_%02d_0_*.sac', ...
                    synmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
                [~, hdr_s, ~, ~, tims_s] = readsac(synfile);
            catch ME
                synfile = cindeks(ls2cell(sprintf('%s%d/*_%04d_0_*.sac', ...
                    synmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
                [~, hdr_s, ~, ~, tims_s] = readsac(synfile);
            end
            
            % filter using the chosen corner frequencies
            if op4 == 2
                pres_o = bandpass(pres_o, fs_o, fcs(jj, 1), fcs(jj, 2), ...
                    4, 2, 'butter', 'linear');
                pres_o = pres_o .* shanning(length(pres_o), 0.05, 0);
            end
            
            % normalize the seismogram to 3.5% of the y-limit
            ep = 0.01/fs_o;
            t_min = window_plot(1);
            t_max = window_plot(2);
            pres_o2 = pres_o(and(geq(tims_o  + tsmul * t_shift(jj), t_min, ep), ...
                leq(tims_o + tsmul * t_shift(jj), t_max, ep)));
            o_norm = 0.035 * ywidth / max(abs(pres_o2));
            amp(jj) = max(abs(pres_o2));
            
            % plot together on a plot
            if op5 == 2 && CCmax(jj) <= 0.6
                color_obs = [0.6 0.8 1];
                color_tick = [0.6 0.6 0.6];
            else
                color_obs = [0 0.2 0.8];
                color_tick = [0 0 0];
            end
            % axes for plotting seismograms
            ax2ss = axes('Position', [0.08 0.08 0.84 0.56]);
            axes_collection = [axes_collection; ax2ss];
            priority_values = [priority_values; CCmax(jj)];
            
            if op2 == 1
                if op1 == 1
                    y_values = gcarc;
                else
                    y_values = azim;
                end
            else
                y_values = 1:length(gcarc);
            end
            
            % add ticks indicating a round trip between surface and bottom
            roundtrip_time = 2 * (-hdr_s.STEL) / 1500;
            ticks_x = 0:roundtrip_time:window_plot(2);
            plot(ax2ss, [ticks_x; ticks_x], repmat(y_values(jj) + ...
                0.035 * ywidth * [-1; 1], 1, length(ticks_x)), ...
                'Color', color_tick, 'LineWidth', 0.75);
            hold on
            
            % plot the seismograms
            signalplot(pres_o * o_norm + y_values(jj), fs_o, tims_o(1), ...
                ax2ss, '', [], color_obs, 'LineWidth', 1);
            
            ax2ss.Title.String = '';
            ax2ss.XAxis.Visible = 'off';
            ax2ss.YAxis.Visible = 'off';
            set(ax2ss, 'Box', 'off', 'TickDir', 'out', ...
                'XLim', window_plot, 'YLim', ylimit, 'FontSize', 12, ...
                'Color', 'none', 'XGrid', 'off', 'YGrid', 'off');
            
            % location of the labels
            if op3 == 1
                [x_norm, y_norm] = true2normposition(ax2ss, -39.0, ...
                    y_values(jj));
            else
                [x_norm, y_norm] = true2normposition(ax2ss, -19.5, ...
                    y_values(jj));
            end
            
            if is_label_left && prev_label_top > y_norm
                if op3 == 1
                    [x_norm, y_norm] = true2normposition(ax2ss, 58.5, ...
                        y_values(jj));
                else
                    [x_norm, y_norm] = true2normposition(ax2ss, 39.5, ...
                        y_values(jj));
                end
                is_label_left = false;
            else
                is_label_left = true;
            end
            prev_label_top = y_norm + 0.075;
            
            % add pressure text
            label_str = sprintf('$$ %.2f\\textnormal{~Pa} $$', ...
                amp(jj));
            
            % I set EdgeColor and FaceColor to 0.999 to avoid the patch
            % turning black when being saved.
            % [TODO] Figure out why it happens
            if is_label_left
                [x_text, y_text] = norm2trueposition(ax2ss, x_norm, ...
                    y_norm + 0.025);
                bhan = textbox(ax2ss, x_text, y_text, label_str, ...
                    {'FontSize', 11.5, 'Interpreter', 'latex'}, ...
                    {'EdgeColor', [1 1 1] * 0.999});
                bhan.FaceColor = [1 1 1] * 0.999;
            else
                [x_text, y_text] = norm2trueposition(ax2ss, x_norm, ...
                    y_norm + 0.025);
                bhan = textbox(ax2ss, x_text, y_text, label_str, ...
                    {'FontSize', 11.5, 'Interpreter', 'latex', ...
                     'HorizontalAlignment', 'right'}, ...
                    {'EdgeColor', [1 1 1] * 0.999});
                bhan.FaceColor = [1 1 1] * 0.999;
            end
            
            % add MERMAID number text (and corner frequencies)
            if op4 == 1
                number_str = sprintf('$$ \\textnormal{%s} $$', ...
                    stationname{jj});
            else
                number_str = sprintf('$$ \\textnormal{%s}, %4.2f-%4.2f\\textnormal{~Hz}$$', ...
                    stationname{jj}, fcs(jj,1), fcs(jj,2));
            end
            
            if is_label_left
                [x_text, y_text] = norm2trueposition(ax2ss, x_norm, ...
                    y_norm - 0.025);
                [bhan, than] = textbox(ax2ss, x_text + 1.5, y_text, number_str, ...
                    {'FontSize', 11.5, 'Interpreter', 'latex'}, ...
                    {'EdgeColor', [1 1 1] * 0.999});
                bhan.FaceColor = [1 1 1] * 0.999;
            else
                [x_text, y_text] = norm2trueposition(ax2ss, x_norm, ...
                    y_norm - 0.025);
                [bhan, than] = textbox(ax2ss, x_text, y_text, number_str, ...
                    {'FontSize', 11.5, 'Interpreter', 'latex', ...
                     'HorizontalAlignment', 'right'}, ...
                    {'EdgeColor', [1 1 1] * 0.999});
                bhan.FaceColor = [1 1 1] * 0.999;
            end
            
            % add MERMAID icon
            % get color icon
            color_icon = cmap(stationid(jj) == metadata.STNM(whevent),:);
            color_txt = [0 0 0];
            if op5 == 2 && CCmax(jj) <= 0.6
                % grey out if the CC is low
                color_icon = 0.5 * color_icon + 0.5;
                color_txt = [0.5 0.5 0.5];
            end
            % determine the icon's position
            icon_xposition = than.Extent(1) - 1;
            [~, icon_yposition] = norm2trueposition(ax2ss, x_norm, ...
                    y_norm - 0.018);
                
            scatter(icon_xposition, icon_yposition, 60, color_icon, 'v', ...
                'filled', 'MarkerEdgeColor', color_txt)
        end
        
        % add azimuth value to right y-axis
        axa = doubleaxes(ax2);
        axes_collection = [axes_collection; axa];
        priority_values = [priority_values; -2];
        
        if op2 == 1
            if op1 == 1
                axa.YTick = gcarc;
                axa.YTickLabel = num2str(round(azim));
                ax2.YGrid = 'on';
            else
                axa.YGrid = 'on';
                ax2.YTick = azim;
                ax2.YTickLabel = num2str(round(gcarc));
            end
        else
            ax2.YTick = 1:length(gcarc);
            ax2.YTickLabel = num2str(round(gcarc));
            axa.YTick = 1:length(gcarc);
            axa.YTickLabel = num2str(round(azim));
            if op1 == 1
                ax2.YGrid = 'on';
            else
                axa.YGrid = 'on';
            end
        end
        axa.XTickLabel = [];
        axa.YLabel.String = 'azimuth (degrees)';
        set(axa, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12, ...
            'Color', 'none');
        
        % rearrange the figure
        [~, i_sort] = sort(priority_values, 'ascend');
        for i_ax = i_sort'
            axes(axes_collection(i_ax));
        end
        
        % figure label
        if isempty(CMT)
            title(ax1, ...
                sprintf('Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', ...
                uniqevent(ii), hdr_o.MAG, hdr_o.EVDP));
        else
            title(ax1, ...
                sprintf('%s, Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', ...
                CMT.EventName, uniqevent(ii), hdr_o.MAG, hdr_o.EVDP));
        end
        % move the title up a little bit
        [ax1.Title.Position(1), ax1.Title.Position(2)] = ...
            norm2trueposition(ax1, 0.5, 1.06);
        
        % save figure
        set(gcf, 'Renderer', 'painters')
        fname = sprintf('%s_%d', mfilename, uniqevent(ii));
        figdisp(fname, [], [], 2, [], 'epstopdf');
    end
end
end

function addcontourlines(ax, lon, lat)
% distant lines
for distant = 10:10:180
    [latout, lonout] = reckon(lat, lon, distant, 0:360);
    lonout = mod(lonout, 360);
    
    % find if the track cross the cut-off longitude
    is_cross = (abs(lonout(2:end) - lonout(1:end-1)) > 90);
    where_cross = find(is_cross > 0);
    % add NaN points at the crossing
    latout = insert(latout, NaN(size(where_cross)), where_cross + 1);
    lonout = insert(lonout, NaN(size(where_cross)), where_cross + 1);
    
    plot(ax, lonout, latout, 'LineWidth', 0.5, ...
        'Color', [0.8 0.8 0.8]);
end

% azimuth lines
for azimuth = 0:30:330
    [latout, lonout] = reckon(lat, lon, 0:180, azimuth);
    lonout = mod(lonout, 360);
    
    % find if the track cross the cut-off longitude
    is_cross = (abs(lonout(2:end) - lonout(1:end-1)) > 90);
    where_cross = find(is_cross > 0);
    % add NaN points at the crossing
    latout = insert(latout, NaN(size(where_cross)), where_cross + 1);
    lonout = insert(lonout, NaN(size(where_cross)), where_cross + 1);
    
    plot(ax, lonout, latout, 'LineWidth', 0.5, ...
        'Color', [0.8 0.8 0.8]);
end
end

function r = leq(a, b, ep)
r = or(a < b, abs(a - b) < ep);
end

function r = geq(a, b, ep)
r = or(a > b, abs(a - b) < ep);
end