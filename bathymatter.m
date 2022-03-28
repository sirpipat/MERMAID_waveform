function bathymatter(obsmasterdir, synmasterdir, flatmasterdir, bathmasterdir)
% BATHYMATTER(obsmasterdir, synmasterdir, flatmasterdir, bathmasterdir)
%
% Compares 2 responses from flat ocean bottom and bathymetry from GEBCO in
% order to determine whether bathymatry matters and when.
%
% INPUT
% obsmasterdir      the master directory to the observed files sorted into
%                   IRIS event ID folders
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% flatmasterdir     the master directory to the output folders from
%                   SPECFEM2D run for fluid-solid setting with flat
%                   topography sorted into IRIS event ID folders
% bathmasterdir     the master directory to the output folders from
%                   SPECFEM2D run for fluid-solid setting with GEBCO
%                   bathymetry sorted into IRIS event ID folders
%
% SEE ALSO
% RUNFLATSIM, COMPARERESPONSEFUNCTIONS
%
% Last modified by sirawich-at-princeton.edu, 03/28/2022

[allflatdirs, fndex] = allfile(flatmasterdir);
[allbathdirs, bndex] = allfile(bathmasterdir);

% loop through allflatdirs
for ii = 1:fndex
    dir = removepath(allflatdirs{ii});
    eventid = cindeks(split(dir, '_'), 2);
    stationid = indeks(cindeks(split(dir, '_'), 3), 4:5);
    obsfile = cindeks(ls2cell(sprintf('%s%s/*.%s_*.sac', ...
        obsmasterdir, eventid, stationid), 1), 1);
    synfile = cindeks(ls2cell(sprintf('%s%s/*_%s_0_*.sac', ...
        synmasterdir, eventid, stationid), 1), 1);
    try
        compareresponsefunctions(obsfile, synfile, [allflatdirs{ii} '/'], ...
            [allbathdirs{ii} '/']);
    catch
        continue
    end
end


end