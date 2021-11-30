function outputdirs = runflatsim_routine(obsmasterdir, synmasterdir, i_begin, i_end)
% outputdirs = RUNFLATSIM_ROUTINE(obsmasterdir, synmasterdir, i_begin, i_end)
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
% Last modified by sirawich-at-princeton.edu, 11/18/2021

defval('obsmasterdir', '/home/sirawich/research/processed_data/MERMAID_reports_updated/')
defval('synmasterdir', '/home/sirawich/research/SYNTHETICS/')

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
        outputdirs = runflatsim(allobsfiles{jj});
        example = sprintf('flat_%d_%s', hdr_o.USER7, ...
            replace(hdr_o.KSTNM, ' ', ''));
        [~, ~, t_rf, rf] = cctransplot(outputdirs{1}, outputdirs{2}, example, ...
            {'bottom', 'displacement'}, {'hydrophone', 'pressure'}, 1/hdr_o.DELTA, true);
        comparepressure(seis_s, hdr_s, seis_o, hdr_o, rf, t_rf);
    end
end


end