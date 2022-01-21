function arrivalstats(sacfiles)
% ARRIVALSTATS(sacfiles)
%
% Plots histograms and scatter plots of epicentral distance, ray parameter,
% MERMAID dpeth, ocean depth, and azimuth of all earthquake arrivals.
%
% INPUT:
% sacfiles      cell array containing SAC files
%
% Last modified by sirawich-at-princeton.edu, 11/29/2021

BADVAL = -12345;

% stats to consider
gcarcs = zeros(size(sacfiles));
ps = zeros(size(sacfiles));
depths = zeros(size(sacfiles));
bottoms = zeros(size(sacfiles));
azs = zeros(size(sacfiles));
srcdepth = zeros(size(sacfiles));
r = presiduestat(sacfiles, false);

for ii = 1:length(sacfiles)
    [~, hdr] = readsac(sacfiles{ii});
    gcarcs(ii) = hdr.GCARC;
    depths(ii) = hdr.STDP;
    bottoms(ii) = -hdr.STEL;
    azs(ii) = hdr.AZ;
    srcdepth(ii) = hdr.EVDP;
    if hdr.USER9 == BADVAL
        % compute theoretical ray parameter at the ocean bottom below MERMAID.
        % [lat lon] of the receiver is slightly shifted if incident angle is not
        % close to zero.
        tt = taupPierce('ak135', HdrData.EVDP, ...
            'p,s,P,S,PP,SS,PKP,SKS,PKIKP,SKIKS', ...
            'sta', [HdrData.STLA HdrData.STLO], ...
            'evt', [HdrData.EVLA HdrData.EVLO], ...
            'pierce', -HdrData.STEL/1000, 'nodiscon');

        % remove all zero piercings
        for jj = 1:length(tt)
            index = length(tt(jj).pierce.p);
            while tt(jj).pierce.time(index) <= 0 && index > 1
                index = index - 1;
            end
            tt(jj).time = tt(jj).pierce.time(index);
            tt(jj).distance = tt(jj).pierce.distance(index);
        end

        % keep only one arrival for each phase
        ph = cell(size(tt));
        for jj = 1:length(ph)
            ph{jj} = tt(jj).phaseName;
        end
        [~, ia] = unique(ph);
        tt = tt(ia);

        % sort the arrivals by time
        tp = zeros(size(tt));
        for jj = 1:length(tp)
            tp(jj) = tt(jj).time;
        end
        [~, is] = sort(tp);
        tt = tt(is);

        ps(ii) = tt(1).rayParam;
    else
        ps(ii) = hdr.USER9;
    end
end

% making some plots
figure
set(gcf, 'Units', 'inches', 'Position', [0 0 16 16])
clf

% stores axes handles
axs = matlab.graphics.axis.Axes.empty;

for ii = 1:6
    for jj = 1:6
        index = 6 * (ii-1) + jj;
        ax = subplot('Position', subplotposition(6,6,index,0.05*[3 3 0 0]));
        axs(ii,jj) = ax;
        
        % defines what to plot for x-axis
        switch jj
            case 1
                xx = gcarcs;
                xxlabel = 'epicentral distance (degrees)';
            case 2
                xx = ps;
                xxlabel = 'ray parameter (rad s)';
            case 3
                xx = r;
                xxlabel = 'residual (s)';
%                 xx = depths;
%                 xxlabel = 'MERMAID depth (m)';
            case 4
                xx = srcdepth;
                xxlabel = 'source depth (km)';
%                 xx = bottoms;
%                 xxlabel = 'bathymetry (m)';
            case 5
                xx = azs;
                xxlabel = 'azimuth (degrees)';
            otherwise
                xx = zeros(size(ps));
                xxlabel = 'counts';
        end
        
        % defines what to plot for y-axis
        switch ii
            case 2
                yy = gcarcs;
                yylabel = 'epicentral distance (degrees)';
            case 3
                yy = ps;
                yylabel = 'ray parameter (rad s)';
            case 4
                yy = r;
                yylabel = 'residual (s)';
%                 yy = depths;
%                 yylabel = 'MERMAID depth (m)';
            case 5
                yy = srcdepth;
                yylabel = 'source depth (km)';
%                 yy = bottoms;
%                 yylabel = 'bathymetry (m)';
            case 6
                yy = azs;
                yylabel = 'azimuth (degrees)';
            otherwise
                yy = zeros(size(ps));
                yylabel = 'counts';
        end
        
        % filter any BADVAL depth
        if ii == 4
            wh = (yy ~= BADVAL);
            yy = yy(wh);
            xx = xx(wh);
        elseif jj == 3
            wh = (xx ~= BADVAL);
            xx = xx(wh);
            yy = yy(wh);
        end
        
        % vertical histograms
        if ii == 1
            switch jj
                case 1
                    histogram(xx)
                    title('Epicentral distance')
                case 2
                    histogram(xx)
                    title('Ray parameter')
                case 3
                    histogram(xx)
                    title('Residual')
                    %title('MERMAID depth')
                case 4
                    histogram(xx)
                    title('Source depth')
                    %title('Bathymetry')
                case 5
                    histogram(xx)
                    title('Azimuth')
                otherwise
                    delete(ax)
                    continue
            end
        % horizontal histograms
        elseif jj == 6
            switch ii
                case 2
                    histogram(yy, 'Orientation', 'horizontal')
                case 3
                    histogram(yy, 'Orientation', 'horizontal')
                case 4
                    histogram(yy, 'Orientation', 'horizontal')
                case 5
                    histogram(yy, 'Orientation', 'horizontal')
                case 6
                    histogram(yy, 'Orientation', 'horizontal')
                otherwise
                    delete(ax)
                    continue
            end
        % scater plot
        else
            scatter(xx, yy, 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
        end
        grid on
        ax.TickDir = 'both';
        ax.Box = 'on';
        
        % axis labels + ticks
        if ii == 6
            xlabel(xxlabel)
            ylim([0 360])
            yticks(0:60:360)
        end
        
        if jj == 1
            ylabel(yylabel)
        elseif jj == 5
            xlim([0 360])
            xticks(0:60:360)
        end
        
        % remove cluttering/redundant axis labels
        if ii < 6 && jj < 6
            nolabels(ax, 1);
        end
        if ii > 1 && jj > 1
            nolabels(ax, 2);
        end
    end
end

% fix the axis limit's alignment
for ii = 1:6
    for jj = 1:6
        if ii > 1 && jj < 6
            axs(ii,jj).XLim = axs(1,jj).XLim;
            axs(ii,jj).YLim = axs(ii,6).YLim;
        end
    end
end

% save figure
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s.eps', mfilename), [], [], 2, [], 'epstopdf');
end