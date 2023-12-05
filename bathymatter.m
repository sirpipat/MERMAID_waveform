function [t_shifts, CCmaxs, fcorners, snr, s, depthstats, slopestats, ...
    peakstats, n, metadata] = ...
    bathymatter(obsmasterdir, synmasterdir, flatmasterdir, ...
    bathmasterdir, opt, plt)
% [t_shifts, CCmaxs, fcorners, snr, s, depthstats, slopestats, ...
%     peakstats, n, metadata] = ...
%     BATHYMATTER(obsmasterdir, synmasterdir, flatmasterdir, ...
%                 bathmasterdir, opt, plt)
%
% Compares 2 responses from flat ocean bottom and bathymetry from GEBCO in
% order to determine whether bathymatry matters and when.
%
% INPUT
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% flatmasterdir     the master directory to the output folders from
%                   SPECFEM2D run for fluid-solid setting with flat
%                   topography sorted into IRIS event ID folders
% bathmasterdir     the master directory to the output folders from
%                   SPECFEM2D run for fluid-solid setting with GEBCO
%                   bathymetry sorted into IRIS event ID folders
% opt               corner frequency options
%                       1             fixed at 1-2 Hz
%                       2             selected by FREQSELECT  [Defatult]
%                       [fc1 fc2]     user defined corner frequencies
% plt               whether to plot or not [Default: true]
%
% OUTPUT
% t_shifts          Best time shift where CC is maximum     [flat, bath]
% CCmaxs            Maximum correlation coefficient         [flat, bath]
% fcorners          corner frequencies used for comparing synthetic and
%                   observed acoustic pressures
% snr               best signal-to-noise ratio on the observed waveform
% s                 scaling to minimize the misfit          [flat, bath]
% depthstats        struct contatining the following fields
%       depth_mid           depth at the middle of the profile right below
%                           the hydrophone
%       depth_avg           mean of the depth of the bathymetry profile
%       depth_std           standard deviation
%       depth_range        difference between the highest and the lowest
% slopestats        struct containing the following fields every field has 
%                   a unit in degree and is positive when sloping upward 
%                   except SLOPE_EXTREME which is always positive
%       slope_full          general trend over the bathymetry profile. It
%                           is calculated from linear fitting.
%       slope_left          same as SLOPE_FULL but using the left side or
%                           the side where the wave is coming from
%       slope_right         same as SLOPE_FULL but using the right side or
%                           the side where the wave is going forward
%       slope_local         slope within 25th and 75th percentile of the
%                           distance along the profile
%       slope_extreme       maximum of absolute slope from moving average
%                           slope over 1000 m section
% peakstats         struct containing the following fields
%       N                   the number of observed peaks
%       prominence_max      maximum prominence of the peaks (0 if N == 0)
%       prominence_avg      average prominence of the peaks (0 if N == 0)
%       width_avg           avarage width of the peaks      (0 if N == 0)
% n                 the number of data points
% metadata          SAC header variables sorted by variable names
%
% SEE ALSO
% RUNFLATSIM, COMPARERESPONSEFUNCTIONS
%
% Last modified by sirawich-at-princeton.edu, 12/05/2023

defval('opt', 2)
defval('true', plt)

[allflatdirs, fndex] = allfile(flatmasterdir);
[allbathdirs, bndex] = allfile(bathmasterdir);

defval('sname', sprintf('%s_%s.mat', mfilename, ...
    hash([obsmasterdir synmasterdir flatmasterdir bathmasterdir ...
    double(plt) fndex bndex sum(opt) nargin nargout], 'SHA-1')))

pname = fullfile(getenv('IFILES'), 'HASHES', sname);

if plt || ~exist(pname, 'file')
    CCmaxs_flat = [];
    CCmaxs_bath = [];
    t_shifts_flat = [];
    t_shifts_bath = [];
    fcorners = [];
    snr = [];
    s = [];
    depth_mid = [];
    depth_avg = [];
    depth_std = [];
    depth_range = [];
    slope_full = [];
    slope_local = [];
    slope_left = [];
    slope_right = [];
    slope_extreme = [];
    N = [];
    prominence_max = [];
    prominence_avg = [];
    width_avg = [];
    fileused = {};
    n = 1;

    % loop through allflatdirs
    for ii = 1:fndex
        dir = removepath(allflatdirs{ii});
        eventid = cindeks(split(dir, '_'), 2);
        stationid = indeks(cindeks(split(dir, '_'), 3), 4:5);
        obsfile = cindeks(ls2cell(sprintf('%s%s/*.%s_*.sac', ...
            obsmasterdir, eventid, stationid), 1), 1);
        synfile = cindeks(ls2cell(sprintf('%s%s/*_%s_0_*.sac', ...
            synmasterdir, eventid, stationid), 1), 1);
        try
            if size(opt, 1) == 1
                [t_shift1, t_shift2, CCmax1, CCmax2, bath1, bath2, fcs, ...
                    snrr, s1, s2] = ...
                    compareresponsefunctions(obsfile, synfile, ...
                    [allflatdirs{ii} '/'], [allbathdirs{ii} '/'], opt, plt);
            else
                [t_shift1, t_shift2, CCmax1, CCmax2, bath1, bath2, fcs, ...
                    snrr, s1, s2] = ...
                    compareresponsefunctions(obsfile, synfile, ...
                    [allflatdirs{ii} '/'], [allbathdirs{ii} '/'], opt(n, :), plt);
            end
            bath1(:,2) = bath1(:,2) - 9600;
            bath2(:,2) = bath2(:,2) - 9600;
            t_shifts_flat(n,1) = t_shift1;
            t_shifts_bath(n,1) = t_shift2;
            fcorners(n,:) = fcs;
            snr(n) = snrr;
            s(n,1) = s1;
            s(n,2) = s2;
            CCmaxs_flat(n,1) = CCmax1;
            CCmaxs_bath(n,1) = CCmax2;
            
            % statistics of the bathymetry
            depth_mid(n,1) = bath2((size(bath2,1)+1)/2,2);
            depth_avg(n,1) = mean(bath2(:,2));
            depth_std(n,1) = std(bath2(:,2));
            depth_range(n,1) = range(bath2(:,2));
            
            % statistics of the slope
            dx = bath2(2,1) - bath2(1,1);
            x_mid   = (0.5 * bath2(1,1) + 0.5 * bath2(size(bath2,1),1));
            x_width = (-bath2(1,1) + bath2(size(bath2,1),1));
            wh_left = (bath2(:,1) <= x_mid);
            wh_right = (bath2(:,1) >= x_mid);
            wh_local = and(bath2(:,1) >= x_mid - 0.5 * x_width/2, ...
                bath2(:,1) <= x_mid + 0.5 * x_width/2);
            
            slope_full(n,1) = atan(indeks(polyfit(bath2(:,1), ...
                bath2(:,2), 1), 1)) * 180/pi;
            slope_left(n,1) = atan(indeks(polyfit(bath2(wh_left,1), ...
                bath2(wh_left,2), 1), 1)) * 180/pi;
            slope_right(n,1) = atan(indeks(polyfit(bath2(wh_right,1), ...
                bath2(wh_right,2), 1), 1)) * 180/pi;
            slope_local(n,1) = atan(indeks(polyfit(bath2(wh_local,1), ...
                bath2(wh_local,2), 1), 1)) * 180/pi;
            
            moving_slope = movmean((bath2(2:size(bath2,1),2) - ...
                bath2(1:size(bath2,1)-1,2)) / dx, 21);
            slope_extreme(n,1) = atan(max(abs(moving_slope))) * 180/pi;
            
            % peaks
            [pks,locs,w,p] = findpeaks(bath2(:,2),bath2(:,1));
            N(n,1) = length(w);
            if isempty(p)
                prominence_max(n,1) = 0;
                prominence_avg(n,1) = 0;
            else
                prominence_max(n,1) = max(p);
                prominence_avg(n,1) = mean(p);
            end
            if isempty(w)
                width_avg(n,1) = 0;
            else
                width_avg(n,1) = mean(w);
            end
            
            fileused{n,1} = synfile;
            n = n + 1;
        catch ME
            fprintf('%s\n', ME.getReport);
            continue
        end
    end

    % fix the number of files
    n = n - 1;

    % gather stats
    t_shifts = [t_shifts_flat, t_shifts_bath];
    CCmaxs = [CCmaxs_flat, CCmaxs_bath];
    depthstats = struct('depth_mid', depth_mid, 'depth_avg', depth_avg, ...
        'depth_std', depth_std, 'depth_range', depth_range);
    slopestats = struct('slope_full', slope_full, ...
        'slope_local', slope_local, 'slope_left', slope_left, ...
        'slope_right', slope_right, 'slope_extreme', slope_extreme);
    peakstats = struct('N', N, 'prominence_max', prominence_max, ...
        'prominence_avg', prominence_avg, 'width_avg', width_avg');

    % gather metadata
    metadata = getheaderarray(fileused);
    
    % save
    fprintf('save the output to a file to %s ...\n', pname);
    save(pname, 't_shifts', 'CCmaxs', 'fcorners', 'snr', 's', ...
        'depthstats', 'slopestats', 'peakstats', 'n', 'metadata');
else
    % load
    fprintf('found the save in a file in %s\n', pname);
    fprintf('load the variables ...\n');
    load(pname, 't_shifts', 'CCmaxs', 'fcorners', 'snr', 's', ...
        'depthstats', 'slopestats', 'peakstats', 'n', 'metadata');
end
end
