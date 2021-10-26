function obssynshift_routine(obsmasterdir, synmasterdir)
% OBSSYNSHIFT_ROUTINE(obsmasterdir, synmasterdir)
%
% A script to plot the best timeshift for all observed-synthetic
% seismogram pairs.
%
% INPUT:
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
%
% SEE ALSO:
% OBSSYNSHIFT
%
% Last modified by sirawich-at-princeton.edu, 10/26/2021

[allsyndirs, dndex] = allfile(synmasterdir);

for ii = 1:dndex      
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
        obssynshift(allobsfiles{jj}, allsynfiles{jj});
    end
end