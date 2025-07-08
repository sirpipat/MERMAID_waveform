function computesyntheticspressure(obsmasterdir, synmasterdir, specmasterdir, presmasterdir)
% COMPUTESYNTHETICSPRSSURE(obsmasterdir, synmasterdir, specmasterdir, presmasterdir)
%
% Compute the synthetic pressure at the MERMAID depth by convolving the
% synthetic vertical displacement at the ocean bottom from Instaseis with
% the response function due to the ocean layer simulated by SPECFEM2D. It
% writes the synthetic pressure to SAC files. The time steps and metadata 
% are the same as the MERMAID record SAC file.
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% specmasterdir     the master directory to the output folders from
%                   SPECFEM2D run for fluid-solid setting sorted into IRIS 
%                   event ID folders
% presmasterdir     the master directory for the synthetic pressure
%                   seismograms will be created and sorted into IRIS eent
%                   ID folders
%
% Last modified by sirawich-at-princeton.edu, 07/08/2025

[allspecdirs, fndex] = allfile(specmasterdir);

% Initialize the output directory if it doesn't exist
if ~exist(presmasterdir, 'dir')
    mkdir(presmasterdir);
end

% loop over SPECFEM2D runs
for ii = 1016:fndex
    dir = removepath(allspecdirs{ii});
    eventid = cindeks(split(dir, '_'), 2);
    % keep all 4-digit station ID at first even if the files may
    % contain only the last 2 digits
    stationid = indeks(cindeks(split(dir, '_'), 3), 2:5);

    % Identify which MERMAID record file to read
    try
        obsfile = cindeks(ls2cell(sprintf('%s%s/*.%s_*.sac', ...
            obsmasterdir, eventid, stationid), 1), 1);
    catch ME
        if strcmp(ME.message, 'This directory or file does not exist')
            obsfile = cindeks(ls2cell(sprintf('%s%s/*.%s_*.sac', ...
                obsmasterdir, eventid, stationid(end-1:end)), 1), 1);
        end
    end

    % read the observed data
    [seis_o, hdr_o] = readsac(obsfile);
    [dt_ref_o, dt_begin_o, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);

    % Identify which synthetic file to read
    try
        synfile = cindeks(ls2cell(sprintf('%s%s/*_%s_0_*.sac', ...
            synmasterdir, eventid, stationid), 1), 1);
    catch ME
        if strcmp(ME.message, 'This directory or file does not exist')
            synfile = cindeks(ls2cell(sprintf('%s%s/*_%s_0_*.sac', ...
                synmasterdir, eventid, stationid(end-1:end)), 1), 1);
        end
    end

    % read the synthetic data
    [seis_s, hdr_s] = readsac(synfile);
    [~, dt_begin_s, ~, fs_s, ~, dts_s] = gethdrinfo(hdr_s);

    % obtain the response function
    try
        [~, ~, t_r, seis_r, d] = cctransplot([allspecdirs{ii} '/'], ...
            [allspecdirs{ii} '/'], [], ...
            {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, ...
            [], fs_o, false);
    catch ME
        ME.getReport
        continue
    end

    % resample to MERMAID datetimes
    seis_s_interp = shannon(dts_s, seis_s, dts_o);
    t_r_interp = 0:(1/fs_o):t_r(end);
    seis_r = shannon(t_r, seis_r, t_r_interp);

    % convolve
    pres_s = conv(seis_s_interp, seis_r);
    pres_s = pres_s(1:length(seis_o), 1);

    % Update the ray parameter of the first arrival phase (USER9)
    % See UPDATESYNTHETICS
    hdr_o.USER9 = hdr_s.USER9;

    % Write the synthetic pressure to SAC files
    % Create a subdirectory for each event ID
    eventDir = fullfile(presmasterdir, eventid);
    if ~exist(eventDir, 'dir')
        mkdir(eventDir);
    end
    % Prepare the output file path
    filename = replace(removepath(synfile), 'SYNTHETIC', ...
        'SYNTHETIC_PRESSURE');
    outputFile = sprintf('%s/%s', eventDir, filename);
    writesac(pres_s, hdr_o, outputFile);
end
end