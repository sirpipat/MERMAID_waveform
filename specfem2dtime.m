function t = specfem2dtime(par_file, opt, source_file)
% t = SPECFEM2DTIME(par_file, opt, source_file)
%
% Figures out the time steps in the seismograms from input files. This
% feature is useful when the seismogram outputs do not contain the time,
% and reading the source-time function file may not be practical due to the
% file size and how the subsampling affect the time in the seismogram.
%
% INPUT:
% par_file      full filename of the parameter file
% opt           definition of the time
%               1 - starting time is at zero [default]
%               2 - peak source-time function of the earliest source is zero
% source_file   full filename of the soruce file
%               [default: 'directory-to-par_file/SOURCE']
%
% OUTPUT:
% t             time steps in the seismogram
%
% SEE ALSO:
% READ_SEISMOGRAM, FREADSEISMOGRAMS
%
% Last modified by sirawich-at-princeton.com, 02/14/2024

% figures out number of samples and time step between two samples
params = loadparfile(par_file);
dt = params.DT;
nsteps = params.NSTEP;
subsample = params.subsamp_seismos;
t = (0:(subsample * dt):((nsteps - 1) * dt))';

% shifts the time so that the zero is the peak of the source-time function
if opt == 2
    defval('source_file', replace(par_file, removepath(par_file), ...
        'SOURCE'));
    sources = loadsource(source_file);
    % figures out the peak of the earliest source
    mintshift = Inf;
    for ii = 1:length(sources)
        tshift = sources{ii}.tshift;
        % adjusts timeshift based on the source-time function
        if any(sources{ii}.source_type == [5 11])
            tshift = tshift - 2.0 / sources{ii}.f0;
        else
            tshift = tshift - 1.2 / sources{ii}.f0;
        end
        if tshift < mintshift
            mintshift = tshift;
        end
    end
    % applies the timeshift
    t = t + mintshift;
end
end