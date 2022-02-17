function outputdirs = runflatsim_routine(obsmasterdir, synmasterdir, i_begin, i_end, is_create, is_run, is_plt, branch)
% outputdirs = RUNFLATSIM_ROUTINE(obsmasterdir, synmasterdir, i_begin, i_end, is_create, is_run, is_plt, branch)
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
%
%   i_begin and i_end must satisfy the following condition
%   1 <= i_begin <= i_end <=dndex    where dndex is the number of IRIS
%   event ID folders in synmasterdir
%
% OUTPUT:
% outputdirs        directories to the two simulations
%       outputdirs{1} -- simulation#1 : pressure receivers
%       outputdirs{2} -- simulation#2 : X/Z displacement receivers
%                   from the last runflatsim call
%
% SEE ALSO:
% RUNFLATSIM, CCTRANSPLOT, COMPAREPRESSURE
%
% Last modified by sirawich-at-princeton.edu, 02/17/2022

defval('obsmasterdir', '/home/sirawich/research/processed_data/MERMAID_reports_updated/')
defval('synmasterdir', '/home/sirawich/research/SYNTHETICS/')
defval('is_create', true)
defval('is_run', true)
defval('is_plt', true)
defval('branch', 'master')

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
            outputdirs = runflatsim(allobsfiles{jj}, [], [], is_run, false, branch);
        % use the existing directories
        else
            outputdirs = cell(2,1);
            outputdirs{1} = sprintf('%s%s_1/', getenv('REMOTE2D'), example);
            outputdirs{2} = sprintf('%s%s_2/', getenv('REMOTE2D'), example);
            
            if is_run
                runthisexample(example, outputdirs{1}, specfembin);
                runthisexample(example, outputdirs{2}, specfembin);
            end
        end
        if is_plt
            [~, ~, t_rf, rf] = cctransplot(outputdirs{1}, outputdirs{2}, example, ...
                {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, 1/hdr_o.DELTA, false);
            try
                comparepressure(seis_s, hdr_s, seis_o, hdr_o, rf, t_rf, ...
                    [-10 20], [-10 5], true, 5, false);
            catch ME
                fprintf('%s\n', ME.getReport);
                fprintf('Error occured. Move on to the next iteration.\n');
            end
        end
    end
end


end