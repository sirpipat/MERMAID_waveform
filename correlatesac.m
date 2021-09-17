function [cc, tr, hdr] = correlatesac(allsacfiles, sortcriteria, plt)
% [ccs, tr, hdr] = CORRELATESAC(allsacfiles, sortcriteria, plt)
%
% Computes and plots the correlation coefficients between seismograms for
% the first 5 seconds after the first picked P-wave arrivals of all
% sac files.
%
% INPUT
% allsacfiles       cell array to full path to all sac files
% sortcriteria      1 -- distance
%                   2 -- azimuth
% plt               whether to plot or not
%
% OUTPUT
% cc                correlation coefficients
% tr                seismic traces, all sorted
% hdr               sac headers, all sorted
%
% Last modified by sirawich-at-princeton.edu, 09/17/2021

% limits of the time window of the seismograms
window_left = 0;
window_right = 5;

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
    
    % cut sections between -10 and 5 s from picked arrival time
    wh = and(tims >= HdrData.T0 + window_left, ...
        tims <= HdrData.T0 + window_right);
    x = x(wh);
    tims = tims(wh) - HdrData.T0;
    
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

if plt
    figure(3)
    set(gcf, 'Units', 'inches', 'Position', [8 8 6 4.8])
    clf
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
    title(sprintf('Event ID: %d, Magnitude: %4.2f, Depth: %6.2f km', eqid, mag, evdp));
    
    filename = sprintf('%s_%d_ccplot.eps', mfilename, eqid);
    figdisp(filename, [], [], 2, [], 'epstopdf');
end