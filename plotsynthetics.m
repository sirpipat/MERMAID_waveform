function plotsynthetics(obsmasterdir, synmasterdir, specmasterdir, ...
    fcorners, CCmaxs, t_shifts, metadata)
% PLOTSYNTHETICS(obsmasterdir, synmasterdir, specmasterdir, fcorners, ...
%     CCmaxs, t_shifts, metadata)
%
% A cross-breed function between PLOTRECORDS and ARRAYCCSHIFTPLOT
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% specmasterdir     the master directory to the output folders from
%                   SPECFEM2D run for fluid-solid setting sorted into IRIS 
%                   event ID folders
% fcorners          corner frequencies used for comparing synthetic and
%                   observed acoustic pressures
% CCmaxs            Maximum correlation coefficient
% t_shifts          Best time shift where CC is maximum
% metadata          SAC header variables sorted by variable names
%
% SEE ALSO:
% PLOTRECORDS, ARRAYCCSHIFTPLOT
%
% Last modified by sirawich-at-princeton.edu, 07/14/2022

%% window lengths
window_waveform = [-5 5];
window_plot = [-40 60];

%% compute for extra metadata
% get station number
metadata.STNM = zeros(size(metadata.T0));
for ii = 1:length(metadata.STNM)
    metadata.STNM(ii) = str2double(indeks(metadata.KSTNM{ii},4:5));
end

% compute relative travel-time difference
metadata.DLNT = t_shifts ./ (metadata.T0 - metadata.USER8);

[uniqevent, ~, ~] = unique(metadata.USER7);

% iterate over all events
for ii = 1:length(uniqevent)
    whevent = (metadata.USER7 == uniqevent(ii));
    
    % only make a plot when there are more than one station
    if sum(whevent) >= 2
        %% list used metadata for making map
        % list of source and receriver locations for an event
        stlo = mod(metadata.STLO(whevent), 360);
        stla = metadata.STLA(whevent);
        evlo = indeks(unique(mod(metadata.EVLO(whevent), 360)), 1);
        evla = indeks(unique(metadata.EVLA(whevent)), 1);
        
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
        
        figure(3)
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
        colormap(gca, jet(n));
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
        
        % add event's focal mechanism
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
        
        %% list metadata for plotting traces
        % filter out the data from other events
        stationid = metadata.STNM(whevent);
        CCmax = CCmaxs(whevent);
        t_shift = t_shifts(whevent);
        fcs = fcorners(whevent, :);
        gcarc = metadata.GCARC(whevent);
        azim = metadata.AZ(whevent);
        dlnt = metadata.DLNT(whevent);
        
        % sort everything by epicentral distance
        [gcarc, i_gcarc] = sort(gcarc);
        stationid = stationid(i_gcarc);
        CCmax = CCmax(i_gcarc);
        t_shift = t_shift(i_gcarc);
        fcs = fcs(i_gcarc, :);
        azim = azim(i_gcarc);
        dlnt = dlnt(i_gcarc);
        
        % sort everything by CCmax to plot 
        
        ax2 = subplot('Position', [0.08 0.08 0.84 0.56]);
        cla
        for jj = 1:sum(whevent)
            % read the observed seismogram
            obsfile = cindeks(ls2cell(sprintf('%s%d/*.%02d_*.sac', ...
                obsmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
            [seis_o, hdr_o, ~, ~, tims_o] = readsac(obsfile);
            [dt_ref_o, dt_begin_o, ~, fs_o] = gethdrinfo(hdr_o);
            tims_o = tims_o - hdr_o.T0;
            pres_o = counts2pa(seis_o, fs_o, [0.05 0.1 10 20], [], 'sacpz', false);
            pres_o = real(pres_o);
            
            % filter
            pres_o = bandpass(pres_o, fs_o, fcs(jj, 1), fcs(jj, 2), ...
                4, 2, 'butter', 'linear');
            pres_o = pres_o .* shanning(length(pres_o), 0.05, 0);
            
            % read the synthetic vertical dispalcement at the ocean bottom
            synfile = cindeks(ls2cell(sprintf('%s%d/*_%02d_0_*.sac', ...
                synmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
            [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
            [~, dt_begin_s, ~, fs_s] = gethdrinfo(hdr_s);
            tims_s = tims_s + seconds(dt_begin_s - dt_begin_o) - hdr_o.T0;
            
            % obtain the response function
            ddir = sprintf('%sflat_%d_P%04d/', specmasterdir, ...
                uniqevent(ii), stationid(jj));
            [~, ~, t_r, seis_r, d] = cctransplot(ddir, ddir, [], ...
                {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, ...
                [], fs_o, false);
            
            % resample to MERMAID datetimes
            seis_s_interp = shannon(tims_s, seis_s, tims_o);
            t_r_interp = 0:(1/fs_o):t_r(end);
            seis_r = shannon(t_r, seis_r, t_r_interp);
            
            % convolve for synthetic pressure seismogram
            pres_s = conv(seis_s_interp, seis_r);
            pres_s = pres_s(1:length(seis_o), 1);
            pres_s = bandpass(pres_s, fs_o, fcs(jj, 1), fcs(jj, 2), ...
                4, 2, 'butter', 'linear');
            pres_s = pres_s .* shanning(length(pres_s), 0.05, 0);
            
            % determine the start and end of the window
            t_min = window_waveform(1); %max(tims_o(1), tims_o(1) + t_shift(jj));
            t_max = window_waveform(2); %min(tims_o(end), tims_o(end) + t_shift(jj));

            % determine scaling
            ep = 0.01/fs_o;
            pres_o1 = pres_o(and(geq(tims_o, t_min, ep), leq(tims_o, t_max, ep)));
            pres_s1 = pres_s(and(geq(tims_o  + t_shift(jj), t_min, ep), ...
                leq(tims_o + t_shift(jj), t_max, ep)));

            s = rms(pres_o1) / rms(pres_s1);
            
            % normalize the seismogram to 1
            t_min = window_plot(1);
            t_max = window_plot(2);
            pres_o2 = pres_o(and(geq(tims_o, t_min, ep), leq(tims_o, t_max, ep)));
            pres_s2 = pres_s(and(geq(tims_o  + t_shift(jj), t_min, ep), ...
                leq(tims_o + t_shift(jj), t_max, ep)));
            s_norm = 1 / max(max(abs(pres_o2)), max(abs(pres_s2 * s)));
            
            % plot together on a plot
            if CCmax(jj) > 0.6
                color_syn = [1 0 0];
                color_obs = [0 0 0];
            else
                color_syn = [1 0.6 0.6];
                color_obs = [0.6 0.6 0.6];
            end
            signalplot(pres_o * s_norm + gcarc(jj), fs_o, tims_o(1), ...
                ax2, '', [], color_obs, 'LineWidth', 1);
            hold on
            signalplot(pres_s * s_norm * s + gcarc(jj), fs_o, ...
                tims_o(1) + t_shift(jj), ax2, '', [], color_syn, ...
                'LineWidth', 1);
        end
        xlabel('time since first picked arrival (s)');
        ylabel('epicentral distance (degrees)');
        
        % adjust y-axis for description labels
        ymin = (min(gcarc) - 2) + (max(gcarc) - min(gcarc) + 2) * 2/93;
        ymax = max(max(gcarc) + 2, ...
            (min(gcarc) - 2) + (max(gcarc) - min(gcarc) + 2) * 98/93);
        ax2.Title.String = '';
        set(ax2, 'Box', 'off', 'TickDir', 'out', 'XLim', window_plot, ...
            'YLim', [ymin ymax], 'FontSize', 12, 'Color', 'none');
        
        % add azimuth value to right y-axis
        axa = doubleaxes(ax2);
        axa.XTickLabel = [];
        axa.YTick = gcarc;
        axa.YTickLabel = num2str(round(azim));
        axa.YLabel.String = 'azimuth (degrees)';
        set(axa, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12, ...
            'Color', 'none');
        
        % highlight the window for corrleation
        axh = doubleaxes(ax2);
        axes(axh)
        [xbox, ybox] = boxcorner(window_waveform, ax2.YLim);
        pgon = polyshape(xbox, ybox);
        bx = plot(axh, pgon, 'FaceColor', [1 0.9 0.4], 'FaceAlpha', 0.4, ...
                'EdgeColor', [0.7 0.6 0.1], 'EdgeAlpha', 1);
        hold on
        axh.XAxis.Visible = 'off';
        axh.YAxis.Visible = 'off';
        set(axh, 'Box', 'on', 'TickDir', 'both', 'XLim', ax2.XLim, ...
            'YLim', ax2.YLim, 'Position', ax2.Position);
        axes(ax2)
        
        % add description labels
        is_label_left = true;
        prev_label_top = 0;
        for jj = 1:sum(whevent)
            % read the observed seismogram
            obsfile = cindeks(ls2cell(sprintf('%s%d/*.%02d_*.sac', ...
                obsmasterdir, uniqevent(ii), stationid(jj)), 1), 1);
            [~, hdr_o] = readsac(obsfile);
            
            [x_norm, y_norm] = true2normposition(ax2, -39.5, gcarc(jj));
            
            if is_label_left && prev_label_top > y_norm
                [x_norm, y_norm] = true2normposition(ax2, 25.5, gcarc(jj));
                is_label_left = false;
            else
                is_label_left = true;
            end
            prev_label_top = y_norm + 0.035;
            
            axb = addbox(ax2, [max(x_norm,0) y_norm+0.01 0.34 0.035]);
            axes(axb)
            if CCmax(jj) > 0.6
                color_txt = [0 0 0];
            else
                color_txt = [0.5 0.5 0.5];
            end
            text(0.01, 0.4, sprintf(['$$ \\textnormal{P%04d,} X(%.2f\\ \\textnormal{s}) = ' ...
                '%.2f, \\Delta \\tau / \\tau = %.2f \\%% $$'], stationid(jj), t_shift(jj), ...
                CCmax(jj), dlnt(jj) * 100), ...
                'Interpreter', 'latex', 'FontSize', 10, 'Color', color_txt);
        end
        
        % move the description labels to the front
        ax2.Parent.Children = ax2.Parent.Children([1 3:(end-5) 2 end-3 ...
            end-4 (end-2):end]);
        
        % rearrange the traces by CCmax
        % traces with higher CCmax are in the front
        % repmat and [0 -1e-8] is for keeping synthetic traces in front of
        % their corresponding observed traces.
        CCmax_ord = reshape(repmat(CCmax,[1 2])', [length(CCmax)*2 1]);
        CCmax_ord = flip(CCmax_ord) + repmat([0 -1e-8]', [length(gcarc) 1]);
        [~,i_sort] = sort(CCmax_ord, 'descend');
        ax2.Children = ax2.Children(i_sort);
        
        % rearange the labels by CCmax
        [~, i_sort] = sort(CCmax, 'descend');
        ax2.Parent.Children(1:length(CCmax)) = ax2.Parent.Children(i_sort);
        
        % figure label
        [~, ~, CMT] = getfocalmech('PublicID', string(uniqevent(ii)));
        if isempty(CMT)
            title(ax1, ...
                sprintf('Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', ...
                uniqevent(ii), hdr_o.MAG, hdr_o.EVDP));
        else
            title(ax1, ...
                sprintf('%s, Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', ...
                CMT.EventName, uniqevent(ii), hdr_o.MAG, hdr_o.EVDP));
        end
        
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