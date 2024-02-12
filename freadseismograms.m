function seis = freadseismograms(fname, nrec, npts, dtype)
% seis = freadseismograms(fname, nrec, npts, dtype)
% seis = freadseismograms(fname, par_file)
%
% Reads seismogram(s) from a seismogram binary file output from SPECFEM2D.
%
% INPUT:
% fname         full filename to the binary file output
% nrec          number of receivers
% npts          number of samples in the seismograms
% dtype         datatype in the binary file ('single' (default) or 'double')
% par_file      full filename of the parameter file (use this input
%               parameter if you want the function to infer NREC, NPTS, 
%               and DTYPE from the parameter file instead)
%
% OUTPUT:
% seis          seismograms with size(SEIS) = [NPTS NREC]
%
% Last modified by sirawich-at-princeton.edu, 02/12/2024


% determine number of receivers and number of samples from Par_file
if ischar(nrec)
    params = loadparfile(nrec);
    nrec = 0;
    for ii = 1:length(params.RECEIVERS)
        nrec = nrec + params.RECEIVERS{ii}.nrec;
    end
    npts = params.NSTEP / params.subsamp_seismos;
end

% try to determine data type from filename first
if ~exist('dtype', 'var') || ~ischar(dtype)
    [~, filename] = fileparts(fname);
    if contains(filename, 'double')
        dtype = 'double';
    elseif contains(filename, 'single')
        dtype = 'single';
% then try to determine data type from parameter file if possible
    elseif exist('params', 'var')
        if params.save_binary_seismograms_single
            dtype = 'single';
        elseif params.save_binary_seismograms_double
            dtype = 'double';
        end
    else
% otherwise assume the datatype to single
        dtype = 'single';
    end
end

% read the seismograms from a binary file
fid = fopen(fname);
seis = fread(fid, [npts nrec], dtype);
fclose(fid);
end
