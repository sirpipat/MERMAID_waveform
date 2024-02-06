function updateheader_from_tomocat(sacfiles, tomofile, ddirout)
% UPDATEHEADER_FROM_TOMOCAT(sacfiles, tomofile, ddirout)
%
% Takes reported SAC files as an input. Then, applies firstarrival picks
% and identified earthquakes to the headers. Finally, saved SAC files with
% updated headers.
%
% INPUT:
% sacfiles     cell array of the original reported SAC files
% tomofile     full filename to the identified earthquakes file
% ddirout      output directory for the updated SAC files
%
% It creates directories <ddirout>/<EventIDs>/ if they do not exist, and
% saves SAC files containing identifiable arrivals from an earthquake with
% EventID. It updates the following header fields:
% - T0          time of the first arrival pick 
% - KT0         name of the first arrival phase
% - AZ          azimuth
% - BAZ         back azimuth
% - DIST        epicentral distance in km
% - EVDP        earthquake depth
% - EVLA        earthquake latitude
% - EVLO        earthquake longitude
% - GCARC       epicentral distance in degrees
% - MAG         magnitude
% - STDP        elevation of the MERMAID
% - STEL        elevation of the seafloor at MERMAID
% - USER4       travel time residue: pick - synthetic (ak135 no bathymetry)
% - USER5       theoretical travel time of the first arrival phase using
%               taupTime.m and ak135 model
% - USER6       Time difference between reference model with bathymetry and
%               reference model w/o bathymetry (bathymetry - no bathymetry)
% - USER7       IRIS event ID
% - USER8       event rupture time relatve to reference time in seconds
%
% SEE ALSO:
% UPDATEHEADER
%
% Last modified by sirawich-at-princeton.edu, 02/06/2024

mermaid = read_tomocat1(tomofile);
for ii = 1:length(mermaid.seismogram_time)
    mermaid.seismogram_time{ii} = indeks(mermaid.seismogram_time{ii}, '1:end-3');
end

for ii = 1:length(sacfiles)
    % read the seismogram and metadata from the original SAC file
    [seis, hdr] = readsac(sacfiles{ii});
    [~, dt_B] = gethdrinfo(hdr);
    
    % determines which tomocat file entry used for updating the header
    % It uses begin time and staion ID to identify the entry
    dt_B.Format = 'uuuu-MM-dd''T''HH:mm:ss';
    %dt_B.Format = 'uuuu-MM-dd''T''HH:mm:ss.SS';
    dt_B_string = string(dt_B);
    begin_time_index = strcmp(mermaid.seismogram_time, dt_B_string);
    % remove the unprintable character and space when compare to the
    % tomocat database
    wh = ismember(hdr.KSTNM, 33:126);
    station_index = strcmp(mermaid.KSTNM, hdr.KSTNM(wh));
    mermaid_index = and(begin_time_index, station_index);
    
    if sum(mermaid_index) ~= 1
        continue
    end
    
    % update the header
    hdr.T0    = mermaid.obs_arvltime(mermaid_index);
    hdr.KT0   = mermaid.phase_name{mermaid_index};
    hdr.EVDP  = mermaid.evdp(mermaid_index) / 1000; % convert to km
    hdr.EVLA  = mermaid.evla(mermaid_index);
    hdr.EVLO  = mermaid.evlo(mermaid_index);
    hdr.GCARC = mermaid.gcarc_1D(mermaid_index);
    hdr.MAG   = mermaid.mag_val(mermaid_index);
    hdr.STDP  = mermaid.stdp(mermaid_index);
    hdr.STEL  = -mermaid.ocdp(mermaid_index);
    hdr.USER4 = mermaid.tres_1D(mermaid_index);
    hdr.USER5 = mermaid.travtime_1D(mermaid_index);
    hdr.USER6 = mermaid.travtime_1Dstar_adj(mermaid_index);
    hdr.USER7 = str2double(mermaid.IRIS_ID{mermaid_index});
    hdr.USER8 = mermaid.obs_arvltime(mermaid_index) - ...
        mermaid.obs_travtime(mermaid_index);
    hdr.AZ = azimuth(hdr.EVLA, hdr.EVLO, hdr.STLA, hdr.STLO);
    hdr.BAZ = azimuth(hdr.STLA, hdr.STLO, hdr.EVLA, hdr.EVLO);
    hdr.DIST = deg2km(hdr.GCARC);
    
    % write the updated sac file
    foldername = mermaid.IRIS_ID{mermaid_index};
    system(sprintf('mkdir -p %s%s', ddirout, foldername));
    writesac(seis, hdr, [ddirout foldername '/' ...
        mermaid.filename{mermaid_index}]);
end
end