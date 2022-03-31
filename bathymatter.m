function [t_shifts, CCmaxs, depthstats, n, metadata] = ...
    bathymatter(obsmasterdir, synmasterdir, flatmasterdir, ...
    bathmasterdir, plt)
% [t_shifts, CCmaxs, depthstats, n, metadata] = ...
%     BATHYMATTER(obsmasterdir, synmasterdir, flatmasterdir, bathmasterdir)
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
% plt               whether to plot or not [Default: true]
%
% OUTPUT
% t_shifts          Best time shift where CC is maximum     [flat, bath]
% CCmaxs            Maximum correlation coefficient         [flat, bath]
% depthstats        struct contatining the following fields
%       - depth_mid         depth at the middle of the profile right below
%                           the hydrophone
%       - depth_avg         mean of the depth of the bathymetry profile
%       - depth_std         standard deviation
%       - depth_range       difference between the highest and the lowest
%       - slope             slope in degree, positive when sloping upward
% n                 the number of data points
% metadata          SAC header variables sorted by variable names
%
% SEE ALSO
% RUNFLATSIM, COMPARERESPONSEFUNCTIONS
%
% Last modified by sirawich-at-princeton.edu, 03/31/2022

defval('true', plt)

[allflatdirs, fndex] = allfile(flatmasterdir);
[allbathdirs, bndex] = allfile(bathmasterdir);

defval('sname', sprintf('%s_%s.mat', mfilename, ...
    hash([obsmasterdir synmasterdir flatmasterdir bathmasterdir ...
    double(plt) fndex bndex], 'SHA-1')))

pname = fullfile(getenv('IFILES'), 'HASHES', sname);

if ~exist(pname, 'file')
    CCmaxs_flat = [];
    CCmaxs_bath = [];
    t_shifts_flat = [];
    t_shifts_bath = [];
    depth_mid = [];
    depth_avg = [];
    depth_std = [];
    depth_range = [];
    slope = [];
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
            [t_shift1, t_shift2, CCmax1, CCmax2, bath1, bath2] = ...
                compareresponsefunctions(obsfile, synfile, ...
                [allflatdirs{ii} '/'], [allbathdirs{ii} '/'], plt);
            bath1(:,2) = bath1(:,2) - 9600;
            bath2(:,2) = bath2(:,2) - 9600;
            t_shifts_flat(n,1) = t_shift1;
            t_shifts_bath(n,1) = t_shift2;
            CCmaxs_flat(n,1) = CCmax1;
            CCmaxs_bath(n,1) = CCmax2;
            depth_mid(n,1) = bath2((size(bath2,1)+1)/2,2);
            depth_avg(n,1) = mean(bath2(:,2));
            depth_std(n,1) = std(bath2(:,2));
            depth_range(n,1) = range(bath2(:,2));
            slope(n,1) = atan(indeks(polyfit(bath2(:,1), ...
                bath2(:,2), 1), 1)) * 180/pi;
            fileused{n,1} = synfile;
            n = n + 1;
        catch ME
            continue
        end
    end

    % fix the number of files
    n = n - 1;

    % gather stats
    t_shifts = [t_shifts_flat, t_shifts_bath];
    CCmaxs = [CCmaxs_flat, CCmaxs_bath];
    depthstats = struct('depth_mid', depth_mid, 'depth_avg', depth_avg, ...
        'depth_std', depth_std, 'depth_range', depth_range, 'slope', slope);

    % gather metadata
    metadata = getheaderarray(fileused);
    
    % save
    fprintf('save the output to a file to %s ...\n', pname);
    save(pname, 't_shifts', 'CCmaxs', 'depthstats', 'n', 'metadata');
else
    % load
    fprintf('found the save in a file in %s\n', pname);
    frpintf('load the variables ...\n');
    load(pname, 't_shifts', 'CCmaxs', 'depthstats', 'n', 'metadata');
end
end