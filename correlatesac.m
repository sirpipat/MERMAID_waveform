function [cc, tr, hdr] = correlatesac(allsacfiles, sortcriteria, arrcriteria, slide, win, plt, fmt)
% [ccs, tr, hdr] = CORRELATESAC(allsacfiles, sortcriteria, arrcriteria, slide, win, plt)
%
% Computes and plots the correlation coefficients between seismograms for
% the windows of all sac files. Use slide and win to adjust the beginning 
% and the end of the window. Slide can be either a scalar or vector.
% For the later case, the correlation coefficient output will be a 3D array
% with following indices: 
%       cc(slide_index, station1_index, station2_index).
% If plt is set to true, a PDF figure and a video are produced when slide
% is a scalar and a vector, respectively.
%
% INPUT
% allsacfiles       cell array to full path to all sac files
% sortcriteria      1 -- distance            [Default]
%                   2 -- azimuth
% arrcriteria       1 -- pick                [Default]
%                   2 -- predicted (ak135)
% slide             amount of window slide   [Default: 0]
% win               window length            [Default: 5]
% plt               whether to plot or not   [Default: true]
%
% OUTPUT
% cc                correlation coefficients
% tr                seismic traces, all sorted
% hdr               sac headers, all sorted
%
% Last modified by sirawich-at-princeton.edu, 09/29/2021

defval('sortcriteria', 1)
defval('arrcriteria', 1)
defval('slide', 0)
defval('win', 5)
defval('plt', true)
defval('fmt', 'eps')

if length(slide) > 1
    images = {};
    % run the first slide value to check the number of valid station
    try
        [cc_slide, tr, hdr] = correlatesac(allsacfiles, sortcriteria, ...
            arrcriteria, slide(1), plt, 'jpg');
        cc = zeros(length(slide), size(cc_slide, 1), size(cc_slide, 2));
        cc(1,:,:) = cc_slide;
        
        % limits of the time window of the seismograms
        window_left = 0 + slide(1);
        window_right = win + slide(1);
        
        eqid = hdr{1,1}.USER7;
        images{1} = sprintf('%s%s_%d%d_%d_%.1fs_%.1fs_ccplot.jpg', ...
            getenv('EPS'), mfilename, sortcriteria, arrcriteria, eqid, ...
            window_left, window_right);
    catch
        return
    end
    % run the remaining slides
    for ii = 2:length(slide)
        try
            [cc_slide, tr, hdr] = correlatesac(allsacfiles, sortcriteria, ...
                arrcriteria, slide(ii), plt, 'jpg');
            cc(ii,:,:) = cc_slide;
            
            % limits of the time window of the seismograms
            window_left = 0 + slide(ii);
            window_right = win + slide(ii);
            
            eqid = hdr{1,1}.USER7;
            images{ii} = sprintf('%s%s_%d%d_%d_%.1fs_%.1fs_ccplot.jpg', ...
                getenv('EPS'), mfilename, sortcriteria, arrcriteria, eqid, ...
                window_left, window_right);
        catch
            continue
        end
    end
    
    % convert saved images to a video
    videofilename = sprintf('%s_%d%d_%d_%.1fs_%.1fs_ccplot', ...
                mfilename, sortcriteria, arrcriteria, eqid, ...
                slide(1), slide(end));
    images2video(images, 'jpg', 10, videofilename);
    
    % remove the intermediate jpg files
    for ii = 1:length(images)
        system(sprintf('rm %s', images{ii}));
    end
    return
else
    % limits of the time window of the seismograms
    window_left = 0 + slide;
    window_right = win + slide;
end

n = length(allsacfiles);
% there is no point to correlate if the number of traces is below 2
if n < 2
    return
end

% store data from all SAC files
tr = cell(n, 1);
hdr = cell(n, 1);
tim = cell(n, 1);
dists = zeros(n, 1);
azs = zeros(n, 1);
sta = cell(n, 1);

% loop through all SAC files
for ii = 1:n 
    [SeisData, HdrData, ~, ~, tims] = readsac(allsacfiles{ii});
    fs = 1/(tims(2) - tims(1));
    
    if isnan(HdrData.T0) || HdrData.T0 == -12345
        tr{ii, 1} = zeros(floor((window_right-window_left) * fs), 1);
        hdr{ii, 1} = HdrData;
        tim{ii, 1} = zeros(floor((window_right-window_left) * fs), 1);
        dists(ii, 1) = NaN;
        azs(ii, 1) = NaN;
        sta{ii, 1} = 'NaN';
        continue
    end
    
    % convert digital counts to pressure
    x = real(counts2pa(SeisData, fs));
    % filter out the ambient noise below 1 Hz
    x = bandpass(x, fs, 1, 2, 2, 1, 'butter', 'linear');
    
    % cut sections between [window_left window_right] from picked arrival time
    if arrcriteria == 1
        wh = and(tims >= HdrData.T0 + window_left, ...
            tims <= HdrData.T0 + window_right);
        x = x(wh);
        tims = tims(wh) - HdrData.T0;
    else
        wh = and(tims >= HdrData.T0 - HdrData.USER4 + window_left, ...
            tims <= HdrData.T0 - HdrData.USER4 + window_right);
        x = x(wh);
        tims = tims(wh) - (HdrData.T0 - HdrData.USER4);
    end
    
    % Remove the last data point in case the number of data points is
    % rounded up. It is essential for correlation coefficient computation.
    x = x(1:floor((window_right-window_left) * fs));
    tims = tims(1:floor((window_right-window_left) * fs));
    
    tr{ii, 1} = x;
    hdr{ii, 1} = HdrData;
    tim{ii, 1} = tims;
    dists(ii, 1) = HdrData.GCARC;
    azs(ii, 1) = HdrData.AZ;
    sta{ii, 1} = HdrData.KSTNM(1:5);
end

% remove invalid data
notNaN = ~isnan(dists);
tr = tr(notNaN);
hdr = hdr(notNaN);
tim = tim(notNaN);
dists = dists(notNaN);
azs = azs(notNaN);
sta = sta(notNaN);
n = length(tr);

% remove redundant stations
[sta, IU, ~] = unique(sta);
tr = tr(IU);
hdr = hdr(IU);
tim = tim(IU);
dists = dists(IU);
azs = azs(IU);
n = length(tr);

% sort everything by distance
if sortcriteria == 1
    [~, I] = sort(dists);
else
    [~, I] = sort(azs);
end
tr = tr(I);
hdr = hdr(I);
tim = tim(I);
dists = dists(I);
azs = azs(I);
sta = sta(I);

% there is no point to correlate if the number of traces is below 2
if n < 2
    return
end

cc = zeros(n, n);
for ii = 1:n
    for jj = 1:n
        cc(ii,jj) = corr(tr{ii, 1}, tr{jj, 1});
    end
end

%% plot the correlation table
if plt
    figure(3)
    set(gcf, 'Units', 'inches', 'Position', [8 8 6 6])
    clf
    subplot('Position', [0.14 0.12 0.8 0.74]);
    cla
    imagesc(cc, [-1 1]);
    colorbar
    colormap('jet');
    set(gca, 'FontSize', 12)
    xticks(1:n)
    yticks(1:n)
    xticklabels(sta);
    xtickangle(45);
    yticklabels(sta);
    
    eqid = hdr{1,1}.USER7;
    mag = hdr{1,1}.MAG;
    evdp = hdr{1,1}.EVDP;
    if arrcriteria == 1
        arrstr = 'picked';
    else
        arrstr = 'ak135';
    end
    title(sprintf('Event ID: %d, window: %.2f to %.2f s around %s arrival', ...
        eqid, window_left, window_right, arrstr));
    
    ax2 = doubleaxes(gca);
    ax2.YAxis.Visible = 'off';
    if sortcriteria == 1
        ax2.XTickLabel = string(round(dists, 1));
    else
        ax2.XTickLabel = string(round(azs, 1));
    end
    ax2.XTickLabelRotation = 45;
    
    filename = sprintf('%s_%d%d_%d_%.1fs_%.1fs_ccplot', ...
        mfilename, sortcriteria, arrcriteria, eqid, window_left, ...
        window_right);
    if strcmp(fmt, 'jpg')
        filename = strcat(filename, '.jpg');
        print([getenv('EPS') filename], '-djpeg', '-r300');
    else
        filename = strcat(filename, '.eps');
        figdisp(filename, [], [], 2, [], 'epstopdf');
    end
end