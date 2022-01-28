function [t_shifts, CCmaxs, scales, n, metadata] = ...
    comparepressure_routine(obsmasterdir, synmasterdir, i_begin, i_end, plt)
% [t_shifts, CCmaxs, scales, n] = ...
%     COMPAREPRESSURE_ROUTINE(obsmasterdir, synmasterdir, i_begin, i_end, plt)
%
% A script for run COMPAREPRESSURE over computed the response functions
% between z-displacement at the ocean bottom and the pressure at the
% hydrophone in the water column generated by RUNFLATSIM.
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% i_begin           first index for IRIS event ID folders
% i_end             last index for IRIS event ID folders
% plt               whether to plot the figure from COMPAREPRESSURE or not
%                   [Default: false]
%
% OUTPUT:
% t_shifts          Best time shift where CC is maximum
% CCmaxs            Maximum correlation coefficient
% scales            Scaling to minimize the misfit
% n                 the number of data points
% metadata          SAC header variables sorted by variable names
%
% SEE ALSO:
% COMPAREPRESSURE
%
% Last modified by sirawich-at-princeton.edu, 01/25/2022

defval('obsmasterdir', '/home/sirawich/research/processed_data/MERMAID_reports_updated/')
defval('synmasterdir', '/home/sirawich/research/SYNTHETICS/')
defval('plt', false)

[allsyndirs, dndex] = allfile(synmasterdir);

% input validation
if i_begin < 1 || i_end > dndex || i_begin > i_end
    error('i_begin and i_end must satisfy the condition: 1 <= i_begin <= i_end <= dndex')
end

CCmaxs = [];
t_shifts = [];
scales = [];
fileused = {};
n = 1;

% loop over IRIS event ID folders
for ii = i_begin:i_end
    evid = removepath(allsyndirs{ii});
    [allobsfiles, ondex] = allfile([obsmasterdir evid '/']);
    [allsynfiles, sndex] = allfile([allsyndirs{ii} '/']);
    % identify available receivers
    receiverid = cell(1,ondex);
    for jj = 1:ondex
        receiverid{jj} = indeks(cindeks(split(cindeks(split(...
            allobsfiles{jj}, '/'), 'end'), '.'), 2), 1:2);
    end
    % keep only unique receivers
    [~, ia, ~] = unique(receiverid);
    allobsfiles = allobsfiles(ia);
    
    % loop over receivers
    for jj = 1:sndex
        % read observed and synthetic SAC files
        [seis_o, hdr_o] = readsac(allobsfiles{jj});
        [seis_s, hdr_s] = readsac(allsynfiles{jj});
        
        % locate the output folders from RUNFLATSIM
        example = sprintf('flat_%d_%s', hdr_o.USER7, ...
            replace(hdr_o.KSTNM, ' ', ''));
        outputdirs = cell(2,1);
        outputdirs{1} = sprintf('%s%s_1/', getenv('REMOTE2D'), example);
        outputdirs{2} = sprintf('%s%s_2/', getenv('REMOTE2D'), example);

        % get best timeshift, CC, and scaling if the output folders exist
        try
            [~, ~, t_rf, rf] = cctransplot(outputdirs{1}, ...
                outputdirs{2}, example, {'bottom', 'displacement'}, ...
                {'hydrophone', 'pressure'}, 1/hdr_o.DELTA, false);
            [t_shifts(n,1), CCmaxs(n,1), ~, ~, scales(n,1)] = ...
                comparepressure(seis_s, hdr_s, seis_o, hdr_o, rf, t_rf, ...
                plt);
            fileused{n,1} = allsynfiles{jj};
            n = n + 1;
        catch ME
            fprintf('%s\n', ME.getReport);
            fprintf('Error occured. Move on to the next iteration.\n');
        end
    end
end
% fix the number of files
n = n - 1;
% gather metadata
metadata = getheaderarray(fileused);
end