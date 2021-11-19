function outputdirs = runflatsim(sacfile, ddir, specfembin)
% outputdirs = RUNFLATSIM(sacfile, ddir, specfembin)
%
% Sets up fluid-solid simulations in SPECFEM2D with a flat bathymetry. 
% Then, runs the simulations and then computes the correlation coefficients
% and transfer functions between the ocean bottom seismometer and the 
% hydrophone right above (optional). For more general purpose fluid-solid
% simulations in SPECFEM2D, please use SPECFEM2D_INPUT_SETUP instead.
%
% INPUT:
% sacfile       SAC file containing source/reciver info. The file should be
%               generated by UPDATEHEADER and UPDATESYNTHETICS in order to 
%               work properly. The station elevation (STEL) is assumed to 
%               be the depth of the hydrophone, and the station depth 
%               (STDP) is assumed to be the ocean bottom elevation.
% ddir          directory to store simulations' files
% specfembin    directory to specfem2d binaries
%
% OUTPUT:
% outputdirs    directories to the two simulations
%       outputdirs{1} -- simulation#1 : pressure receivers
%       outputdirs{2} -- simulation#2 : X/Z displacement receivers
%
% SEE ALSO:
% SPECFEM2D_INPUT_SETUP, RUNTHISEXAMPLE, UPDATEHEADER, UPDATESYNTHETICS
%
% Last modified by sirawich-at-princeton.edu, 11/18/2021

% specify where you want to keep the simulations input/output files
defval('ddir', getenv('REMOTE2D'))
% specify where you keep SPECFEM2D software
defval('specfembin', strcat(getenv('SPECFEM2D'), 'bin/'))

% bad value in SAC files
badval = -12345;

% read the sac file
[~, HdrData] = readsac(sacfile);
fs = 1 / HdrData.DELTA;

%% compute the incident angle
if HdrData.USER9 ~= badval
    theta = atan(HdrData.USER9 * 3.4 / (6371-8.88)) * 180 / pi;
else
    % compute theoretical travel times at the ocean bottom below MERMAID.
    % [lat lon] of the receiver is slightly shifted if incident angle is not
    % close to zero.
    tt = taupPierce('ak135', HdrData.EVDP, ...
        'p,s,P,S,PP,SS,PKP,SKS,PKIKP,SKIKS', ...
        'sta', [HdrData.STLA HdrData.STLO], ...
        'evt', [HdrData.EVLA HdrData.EVLO], ...
        'pierce', -HdrData.STEL/1000, 'nodiscon');

    % remove all zero piercings
    for ii = 1:length(tt)
        index = length(tt(ii).pierce.p);
        while tt(ii).pierce.time(index) <= 0 && index > 1
            index = index - 1;
        end
        tt(ii).time = tt(ii).pierce.time(index);
        tt(ii).distance = tt(ii).pierce.distance(index);
    end
    
    % keep only one arrival for each phase
    ph = cell(size(tt));
    for ii = 1:length(ph)
        ph{ii} = tt(ii).phaseName;
    end
    [~, ia] = unique(ph);
    tt = tt(ia);

    % sort the arrivals by time
    tp = zeros(size(tt));
    for ii = 1:length(tp)
        tp(ii) = tt(ii).time;
    end
    [~, is] = sort(tp);
    tt = tt(is);
    
    theta = atan(tt(1).rayParam * 3.4 / (6371-8.88)) * 180 / pi;
end

% create a name for the output directories
example = sprintf('flat_%d_%s', HdrData.USER7, ...
    replace(HdrData.KSTNM, ' ', ''));
outputdir = sprintf('%s%s/', ddir, example);

if HdrData.STDP == badval
    depth = 1500;
else
    depth = HdrData.STDP;
end

if HdrData.STEL == badval
    [lons,lats,elev,~,~] = bathymetry([], [-0.1 0.1] + HdrData.STLO, ...
        [-0.1 0.1] + HdrData.STLA, false, []);
    bottom = interp2(lats, lons, elev, HdrData.STLA, ...
        mod(HdrData.STLO, 360));
else
    bottom = HdrData.STEL;
end

%% create the input files and directories for the SPECFEM2D runs
outputdirs = specfem2d_input_setup_flat(example, -bottom, ...
    depth, 'homogeneous', 1, theta, [], outputdir);

%% run the simulation
poolobj = parpool('local', 2);
parfor ii = 1:2
    runthisexample(example, outputdirs{ii}, specfembin);
end
delete(poolobj)
%% analyze the data
% for SYNTHETIC output in dislacement vs. OBSERVED MERMAID pressure
cctransplot(outputdirs{1}, outputdirs{2}, example, ...
    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, fs, true);

% for response function at the ocean bottom
cctransplot(outputdirs{1}, outputdirs{2}, example, ...
    {'bottom', 'displacement'}, {'bottom', 'pressure'}, fs, true);

% displacment to pressure at the hydrophone
cctransplot(outputdirs{1}, outputdirs{2}, example, ...
    {'hydrophone', 'displacement'}, {'hydrophone', 'pressure'}, fs, true);

% for pressure propagation from the bottom to the hydrophone
cctransplot(outputdirs{1}, outputdirs{2}, example, ...
    {'bottom', 'pressure'}, {'hydrophone', 'pressure'}, fs, true);

% for reflection pattern
cctransplot(outputdirs{1}, outputdirs{2}, example, ...
    {'hydrophone', 'pressure'}, {'hydrophone', 'pressure'}, fs, true);
end