function indexlist = indexlistmaker(obsmasterdir, synmasterdir, metadata)
% indexlist = INDEXLISTMAKER(obsmasterdir, synmasterdir, metadata)
%
% Determines the indices for SAC files for observed waveforms (obs) and
% SAC files for synthetic waveforms (syn) that are used in MERMAID waveform
% modeling where metadata of the runs are stored.
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% metadata          the SAC header structure array like one from
%                   BATHYMATTER or GETHEADERARRAY
%
% OUTPUT:
% indexlist         a struct with following variables
%   obs                 indices for SAC files for observed waveforms
%   syn                 indices for SAC files for synthetic waveforms
%
% SEE ALSO:
% BATHYMATTER, GETHEADERARRAY
%
% Last modified by sirawich-at-princeton.edu, 03/06/2024

[allobsfiles, ondex] = allfilen(obsmasterdir, 2);
[allsynfiles, sndex] = allfilen(synmasterdir, 2);

allobsfiles = allobsfiles';
allsynfiles = allsynfiles';

% stores PublicID for each MERMAID report SAC file
events_obs = nan(size(allobsfiles));
% stores station number for each MERMAID report SAC file
stnms_obs = nan(size(allobsfiles));
for ii = 1:ondex
    words = split(allobsfiles{ii}, '/');
    events_obs(ii) = str2double(words{end-1});
    stnms_obs(ii) = str2double(cindeks(split(cindeks(split(words{end}, ...
        '.'), 2), '_'), 1));
end

% stores PublicID for each MERMAID synthetic SAC file
events_syn = nan(size(allsynfiles));
% stores station number for each MERMAID synthetic SAC file
stnms_syn = nan(size(allsynfiles));
for ii = 1:sndex
    words = split(allsynfiles{ii}, '/');
    events_syn(ii) = str2double(words{end-1});
    stnms_syn(ii) = str2double(cindeks(split(words{end}, '_'), 2));
end

indexlist.obs = nan(size(metadata.KSTNM));
indexlist.syn = nan(size(metadata.KSTNM));
% find the index
for ii = 1:length(metadata.KSTNM)
    publicid = metadata.USER7(ii);
    stnm = str2double(indeks(metadata.KSTNM{ii}, 2:5));
    
    wh_event_obs = (events_obs == publicid);
    wh_stnm_obs = (stnms_obs == stnm);
    [maxval, maxindex] = max(and(wh_event_obs, wh_stnm_obs));
    if maxval
        indexlist.obs(ii) = maxindex;
    end
    
    wh_event_syn = (events_syn == publicid);
    wh_stnm_syn = (stnms_syn == stnm);
    [maxval, maxindex] = max(and(wh_event_syn, wh_stnm_syn));
    if maxval
        indexlist.syn(ii) = maxindex;
    end
end