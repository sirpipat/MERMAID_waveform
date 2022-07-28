function outputdirs = runflatsim_routine(obsmasterdir, synmasterdir, i_begin, i_end, is_create, is_run, is_plt, branch, gpu_mode)
% outputdirs = RUNFLATSIM_ROUTINE(obsmasterdir, synmasterdir, i_begin, i_end, is_create, is_run, is_plt, branch, gpu_mode)
%
% A script for run fluid-solid simulation to find the response function
% between z-displacement at the ocean bottom and the pressure at the
% hydrophone in the water column. Then, it convolves with the synthetic
% seismograms to compare with the observed pressure record from MERMAID.
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% i_begin           first index for IRIS event ID folders
% i_end             last index for IRIS event ID folders
% is_create         whether to create input files instead of using existing
%                   dirctories [default: true]
% is_run            whether to (re)run runflatsim or just plot the result
%                   [default: true]
% is_plt            whether to plot or not [default: true]
% branch            SPECFEM2D branch [default: 'master']
%                   'master' (commit: e937ac2f74f23622f6ebbc8901d30fb33c1a2c38)
%                   'devel'  (commit: cf89366717d9435985ba852ef1d41a10cee97884)
% gpu_mode          whether to enable GPU MODE [Default: false]
%
%   i_begin and i_end must satisfy the following condition
%   1 <= i_begin <= i_end <=dndex    where dndex is the number of IRIS
%   event ID folders in synmasterdir
%
% OUTPUT:
% outputdirs        directories to the two simulations
%       --- 'master' branch ---
%       outputdirs{1} -- simulation#1 : pressure receivers
%       outputdirs{2} -- simulation#2 : X/Z displacement receivers
%       --- 'devel' branch ---
%       outputdirs    -- both pressures and displacements at OBS and
%                        hydrophone
%                   from the last runflatsim call
%
% SEE ALSO:
% RUNFLATSIM, CCTRANSPLOT, COMPAREPRESSURE
%
% Last modified by sirawich-at-princeton.edu, 07/28/2022

defval('obsmasterdir', '/home/sirawich/research/processed_data/MERMAID_reports_updated/')
defval('synmasterdir', '/home/sirawich/research/SYNTHETICS/')
defval('is_create', true)
defval('is_run', true)
defval('is_plt', true)
defval('branch', 'master')
defval('gpu_mode', false)

badval = -12345;

[allsyndirs, dndex] = allfile(synmasterdir);

% input validation
if i_begin < 1 || i_end > dndex || i_begin > i_end
    error('i_begin and i_end must satisfy the condition: 1 <= i_begin <= i_end <= dndex')
end

for ii = i_begin:i_end
    evid = removepath(allsyndirs{ii});
    [allobsfiles, ondex] = allfile([obsmasterdir evid '/']);
    [allsynfiles, sndex] = allfile([allsyndirs{ii} '/']);
    receiverid = cell(1,ondex);
    for jj = 1:ondex
        receiverid{jj} = indeks(cindeks(split(cindeks(split(...
            allobsfiles{jj}, '/'), 'end'), '.'), 2), 1:2);
    end
    [~, ia, ~] = unique(receiverid);
    allobsfiles = allobsfiles(ia);
    for jj = 1:sndex
        [seis_o, hdr_o] = readsac(allobsfiles{jj});
        [seis_s, hdr_s] = readsac(allsynfiles{jj});
        example = sprintf('flat_%d_%s', hdr_o.USER7, ...
            replace(hdr_o.KSTNM, ' ', ''));
        % create directories for files I/O for SPECFEM2D
        if is_create
            use_bathymetry = true;
            if ~use_bathymetry
                outputdirs = runflatsim(allobsfiles{jj}, [], [], is_run, false, branch, gpu_mode);
            else
                % use GEBCO bathymetry instead
                [x, z] = bathymetryprofile(20000, 401, ...
                    [hdr_s.STLO hdr_s.STLA], mod(180 + hdr_s.BAZ, 360));
                tparams.X = x;
                tparams.Z = z + 9600;
                
                % compute the incident angle
                if hdr_s.USER9 ~= badval
                    theta = atan(hdr_s.USER9 * 3.4 / (6371-8.88)) * 180 / pi;
                else
                    % compute theoretical travel times at the ocean bottom below MERMAID.
                    % [lat lon] of the receiver is slightly shifted if incident angle is not
                    % close to zero.
                    tt = taupPierce('ak135', hdr_s.EVDP, ...
                        'p,s,P,S,PP,SS,PKP,SKS,PKIKP,SKIKS', ...
                        'sta', [hdr_s.STLA hdr_s.STLO], ...
                        'evt', [hdr_s.EVLA hdr_s.EVLO], ...
                        'pierce', -hdr_s.STEL/1000, 'nodiscon');

                    % remove all zero piercings
                    for kk = 1:length(tt)
                        index = length(tt(kk).pierce.p);
                        while tt(kk).pierce.time(index) <= 0 && index > 1
                            index = index - 1;
                        end
                        tt(kk).time = tt(kk).pierce.time(index);
                        tt(kk).distance = tt(kk).pierce.distance(index);
                    end

                    % keep only one arrival for each phase
                    ph = cell(size(tt));
                    for kk = 1:length(ph)
                        ph{kk} = tt(kk).phaseName;
                    end
                    [~, ia] = unique(ph);
                    tt = tt(ia);

                    % sort the arrivals by time
                    tp = zeros(size(tt));
                    for kk = 1:length(tp)
                        tp(kk) = tt(kk).time;
                    end
                    [~, is] = sort(tp);
                    tt = tt(is);

                    theta = atan(tt(1).rayParam * 3.4 / (6371-8.88)) * 180 / pi;
                end
                
                if hdr_s.STDP == badval
                    depth = 1500;
                else
                    depth = hdr_s.STDP;
                end
                
                outputdir = sprintf('%sAK135_RUNS_BATHYMETRY_MUNK/%s/', ...
                    getenv('REMOTE2D'), example);
                outputdirs = specfem2d_input_setup_response(example, ...
                    'custom', tparams, depth, 'munk', ...
                    'homogeneous', 1, theta, [], outputdir, false, branch, ...
                    gpu_mode);
                
                % plot the bathymetry
                try
                    if strcmp(branch, 'master')
                        for kk = 1:length(outputdirs)
                            plotoutput(outputdirs{kk}, example);
                            plotoutput(outputdirs{kk}, example);
                        end
                    else
                        plotoutput(outputdir, example);
                    end
                catch ME
                    delete(gcf)
                end
                if is_run
                    if strcmpi(branch, 'master')
                        runthisexample(example, outputdirs{1}, specfembin);
                        runthisexample(example, outputdirs{2}, specfembin);
                    else
                        runthisexample(example, outputdirs, specfembin);
                    end
                end
            end
        % use the existing directories
        else
            if strcmpi(branch, 'master')
                outputdirs = cell(2,1);
                outputdirs{1} = sprintf('%s%s_1/', getenv('REMOTE2D'), example);
                outputdirs{2} = sprintf('%s%s_2/', getenv('REMOTE2D'), example);
            else
                outputdirs = sprintf('%s%s/', getenv('REMOTE2D'), example);
            end
            
            % (re)run the existing directories
            if is_run 
                if strcmpi(branch, 'master')
                    runthisexample(example, outputdirs{1}, specfembin);
                    runthisexample(example, outputdirs{2}, specfembin);
                else
                    runthisexample(example, outputdirs, specfembin);
                end
            end
        end
        if is_plt
            if strcmpi(branch, 'master')
                [~, ~, t_rf, rf] = cctransplot(outputdirs{1}, outputdirs{2}, example, ...
                    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, [], 1/hdr_o.DELTA, false);
            else
                [~, ~, t_rf, rf] = cctransplot(outputdirs, outputdirs, example, ...
                    {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, [], 1/hdr_o.DELTA, false);
            end
            try
                comparepressure(seis_s, hdr_s, seis_o, hdr_o, rf, t_rf, ...
                    [-10 20], [-10 5], [0.4 1], true, 5, false);
            catch ME
                fprintf('%s\n', ME.getReport);
                fprintf('Error occured. Move on to the next iteration.\n');
            end
        end
    end
end


end