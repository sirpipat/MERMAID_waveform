function updateheader(ddirin, arrfile, evtfile, ddirout)
% UPDATEHEADER(ddirin, arrfile, evtfile, ddirout)
%
% Takes reported SAC files as an input. Then, applies firstarrival picks
% and identified earthquakes to the headers. Finally, saved SAC files with
% updated headers.
%
% INPUT:
% ddirin       directory for the original reported SAC files
% arrfile      full filename to the first arrival picks file
% evtfile      full filename to the identified earthquakes file
% ddirout      output directory for the updated SAC files
%
% It creates directories <ddirout>/<EventIDs>/ if they do not exist, and
% saves SAC files containing identifiable arrivals from an earthquake with
% EventID. It updates the following header fields:
% - T0          time of the first arrival pick 
% - KT0         name of the first arrival phase
% - T1          time of the first P-wave arrival at ocean bottom
% - KT1         name of the first P-wave arrival at ocean bottom
% - AZ          azimuth
% - BAZ         back azimuth
% - DIST        epicentral distance in km
% - EVDP        earthquake depth
% - EVLA        earthquake latitude
% - EVLO        earthquake longitude
% - GCARC       epicentral distance in degrees
% - MAG         magnitude
% - STEL        elevation of the seafloor at MERMAID
% - USER4       travel time residue: pick - synthetic
% - USER5       theoretical travel time of the first arrival phase using
%               taupTime.m and ak135 model w/o bathymetry
% - USER6       Time difference between reference model with bathymetry and
%               reference model w/o bathymetry -- Note that the theoretical
%               trevel time with ak135 reference model w/ bathymetry is
%               USER5 + USER6
% - USER7       IRIS event ID
% - USER8       event rupture time relatve to reference time in seconds
% - USER9       ray parameter of the first P-wave arrival phase (KT1)
%
% Last modified by sirawich-at-princeton.edu, 03/18/2024

%[allfiles, fndex] = gatherrecords(ddirin, [], [], 'sac', []);
[allfiles, fndex] = allfilen('/Users/sirawich/research/processed_data/MERMAID_reports_updated/', 2);

[sac, eqtime, eqlat, eqlon, eqregion, eqdepth, eqdist, eqmag, ...
    eqphase1, eqid, sacdate, eqdate] = readevt2txt(evtfile, [], ...
    [], 'ALL', true);

[s, ph, dat, tres, tptime, tadj, delay, twosd, maxc_y, SNR, ID, ...
    winflag, tapflag, zerflag] = readfirstarrival(arrfile);
    

for ii = 1:fndex
    % load the SAC file we want to update
    [SeisData, HdrData, ~, ~, ~] = readsac(allfiles{ii});
    [dt_ref, ~, ~, ~, ~, ~, ~] = gethdrinfo(HdrData);
    
    % determine which entries from the first arrival pick files and event
    % files for the update
    firstarrival_index = strcmp(s, removepath(allfiles{ii}));
    event_index = strcmp(sac, removepath(allfiles{ii}));
    
    %% update header
    if any(firstarrival_index)
        HdrData.T0 = dat(firstarrival_index);   % first arrival time
        HdrData.KT0 = ph{firstarrival_index};   % first arrival phase
        HdrData.USER4 = tres(firstarrival_index);
        HdrData.USER5 = tptime(firstarrival_index);
        HdrData.USER6 = tadj(firstarrival_index);
    end
    if any(event_index)
        evla = eqlat(event_index);
        evlo = eqlon(event_index);
        stla = HdrData.STLA;
        stlo = HdrData.STLO;
        publicid = eqid{event_index};
        % Sometimes publicid contains '*' in front. See evt2txt.m for
        % detail.
        if strcmp(publicid(1), '*')
            publicid = publicid(2:end);
        end
        HdrData.AZ = azimuth(evla, evlo, stla, stlo);
        HdrData.BAZ = azimuth(stla, stlo, evla, evlo);
        HdrData.DIST = deg2km(eqdist(event_index));
        HdrData.EVDP = eqdepth(event_index);
        HdrData.EVLA = eqlat(event_index);
        HdrData.EVLO = eqlon(event_index);
        HdrData.GCARC = eqdist(event_index);
        HdrData.MAG = eqmag(event_index);
        HdrData.USER7 = str2double(publicid);
        dt_event = datetime(eqdate(event_index));
        HdrData.USER8 = seconds(dt_event - dt_ref);
        foldername = publicid;
    else
        foldername = 'notevent';
    end
    [lons,lats,elev,~,~] = bathymetry([], [-0.1 0.1] + HdrData.STLO, ...
        [-0.1 0.1] + HdrData.STLA, false, []);
    HdrData.STEL = interp2(lats, lons, elev, HdrData.STLA, ...
        mod(HdrData.STLO, 360));
    
    %% determine the ray parameter of the first arrival phase
    % Try to use TauP 1.3.0 or newer where we can specify receiver depth
    if any(event_index)
        try
            tt = tauptime('mod', 'ak135', ...
                          'depth', HdrData.EVDP, ...
                          'gcarc', HdrData.GCARC, ...
                          'ph', 'p,P,PP,PKP,PKIKP', ...
                          'stadepth', -HdrData.STEL/1000);
            new_version = true;
        catch ME
        % OLD VERSION
        % compute theoretical travel times at the ocean bottom below MERMAID.
        % [lat lon] of the receiver is slightly shifted if incident angle is not
        % close to zero.
            tt = taupPierce('ak135', HdrData.EVDP, ...
                'p,P,PP,PKP,PKIKP', ...
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
            new_version = false;
        end

        % keep only one arrival for each phase
        phases = cell(size(tt));
        for jj = 1:length(phases)
            if new_version
                tt(jj).phaseName = tt(jj).phase;
            end
            phases{jj} = tt(jj).phaseName;
        end
        [~, ia] = unique(phases);
        tt = tt(ia);

        % sort the arrivals by time
        tp = zeros(size(tt));
        for jj = 1:length(tp)
            tp(jj) = tt(jj).time;
        end
        [~, is] = sort(tp);
        tt = tt(is);

        if new_version
            HdrData.USER9 = tt(1).rayparameter;
        else
            HdrData.USER9 = tt(1).rayParam;
        end

        HdrData.T1 = tt(1).time;
        HdrData.KT1 = tt(1).phaseName;
    end
    
    %% write the updated sac file
    system(sprintf('mkdir -p %s%s', ddirout, foldername));
    sacfile = strcat(ddirout, foldername, '/', removepath(allfiles{ii}));
    writesac(SeisData, HdrData, sacfile);
end
