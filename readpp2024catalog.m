function catalog = readpp2024catalog(fname)
% catalog = READPP2024CATALOG(fname)
%
% Read the catalog from a catalog file.
%
% INPUT:
% fname         catalog filename
%
% OUTPUT:
% catalog       a struct containing these variables
% - IRIS_ID         IRIS Event ID
% - CMT_ID          Global CMT Solution ID
% - EVLA            event latitude
% - EVLO            event longitude
% - EVDP            event depth (km)
% - ORIGIN_DATE     event origin date 'uuuu-MM-dd'
% - ORIGIN_TIME     event origin time 'HH:mm:ss'
% - MAG_VAL         event magnitude
% - MAG_TYPE        event magnitude type e.g. 'Mww'
% - STLA            station latitude
% - STLO            station longitude
% - STDP            station depth (m), below sea level
% - BATH            bathymetry depth at the station
% - KSTNM           station name
% - PHASE           phase name of the first arrival
% - RAYPARAM        ray parameter of the first arrival phase
% - GCARC           great-circle epicentral distance
% - AZ              azimuth (of the station from the source)
% - BAZ             backazimuth (azimuth of the source from the station)
% - FC_LOWER        chosen lower corner frequency from FREQSELECT
% - FC_UPPER        chosen upper corner frequency from FREQSELECT
% - SNR_ALL         signal-to-noise ratio of the record bandpass 0.4-10 Hz
% - SNR_PASS        signal-to-noise ratio of the record bandpass 
%                   FC_LOWER-FC_UPPER
% - SNR_STOP        signal-to-noise ratio of the record bandstop
%                   FC_LOWER-FC_UPPER
% - CC              optimal correlation coefficeint between the synthetic
%                   and the observed pressure waveform
% - TSHIFT          time shift applied to the synthetic waveform to get the
%                   optimal correlation coefficient
% - TSHIFT_REL      time shift divided by the ray-theoretical travel time
% - ADJUSTMENT      Instaseis ak135f_1s pick arrival minus ray-theoretical
%                   arrival time based on the same Earth model
%
% SEE ALSO:
% WRITEPP2024CATALOG, FREQSELECT, PRESIDUESTAT, COMPAREPRESSURE
%
% Last modified by sirawich-at-princeton.edu, 03/06/2023

% header names
fid = fopen(fname, 'r');
fgetl(fid);
line = fgetl(fid);
fclose(fid);

hdr = split(line);
hdr = hdr(2:end);

% specify formats                    ;
irisid_fmt        = '%8d      '      ;
cmtid_fmt         = '%14s      '     ;
evla_fmt          = '%9.4f      '    ;
evlo_fmt          = evla_fmt         ;
evdp_fmt          = evla_fmt         ;
origin_date_fmt   = '%10s      '     ;
origin_time_fmt   = '%8s      '      ;
mag_val_fmt       = '%6.4f      '    ;
mag_type_fmt      = '%3s      '      ;

stla_fmt          = evla_fmt         ;
stlo_fmt          = evla_fmt         ;
stdp_fmt          = '%6d      '      ;
bath_fmt          = '%5d      '      ;
kstnm_fmt         = '%5s      '      ;
phase_fmt         = '%5s      '      ;
rayparam_fmt      = '%9.4f      '    ;
gcarc_fmt         = '%8.4f      '    ;
az_fmt            = gcarc_fmt        ;
baz_fmt           = gcarc_fmt        ;

fcorner_lower_fmt = '%4.2f      '    ;
fcorner_upper_fmt = fcorner_lower_fmt;
snr_all_fmt       = '%11.4f      '   ;
snr_pass_fmt      = snr_all_fmt      ;
snr_stop_fmt      = snr_all_fmt      ;

cc_fmt            = '%6.4f      '    ;
tshift_fmt        = '%7.2f      '    ;
tshift_ratio_fmt  = '%9.4f      '    ;
correction_fmt    = '%6.2f'          ;

fmt =  [irisid_fmt ...
        cmtid_fmt ...
        evla_fmt ...
        evlo_fmt ...
        evdp_fmt ...
        origin_date_fmt ...
        origin_time_fmt ...
        mag_val_fmt ...
        mag_type_fmt ...
        ...
        stla_fmt ...
        stlo_fmt ...
        stdp_fmt ...
        bath_fmt ...
        kstnm_fmt ...
        phase_fmt ...
        rayparam_fmt ...
        gcarc_fmt ...
        az_fmt ...
        baz_fmt ...
        ...
        fcorner_lower_fmt ...
        fcorner_upper_fmt ...
        snr_all_fmt ...
        snr_pass_fmt ...
        snr_stop_fmt ...
        ...
        cc_fmt ...
        tshift_fmt ...
        tshift_ratio_fmt ...
        correction_fmt ...
        '\n'];
    
fid = fopen(fname, 'r');
C = textscan(fid, fmt,'HeaderLines', 3, ...
             'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fid);

for ii = 1:length(hdr)
    catalog.(hdr{ii}) = C{ii};
end
end