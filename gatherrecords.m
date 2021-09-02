function [allfiles, fndex] = gatherrecords(ddir, dt_begin, dt_end, format, tipe)
% [allfiles, fndex] = GATHERRECORDS(ddir, dt_begin, dt_end, format, tipe)
%
% Collects all files of MERMAID reported seismograms from a local data
% storage. (TODO: pull from the online server)
%
% INPUT
% ddir          where reported seismograms are stored (backslash needed)
% dt_begin      beginning datetime
% dt_end        ending datetime
% format        file format to look for
%               'sac'   (default)
%               'seed'
% tipe          type of report
%               'DET'   determined
%               'REQ'   requested
%               If you want both, leave it empty.
%
% OUTPUT
% allfiles      bottom-file list with complete file names
% fndex         the total number of elements in the list
%
% Last modified by sirawich-at-princeton.edu, 08/25/2021

defval('ddir', '/Users/sirawich/research/MERMAID_REPORT/forPete/processed/')
defval('format', 'sac')

fndex = 0;
% loop through the directory for instrument subdirectories
[mdirs, mndex] = allfile(ddir);
for ii = 1:mndex
    % loop through the instrument subdirctory
    [sdirs, sndex] = allfile([mdirs{ii} '/']);
    for jj = 1:sndex
        s = removepath(sdirs{jj});
        % check if a element is the subdirectory for reported sesimograms
        % and not the metadata file.
        if strlength(s) < 18
            continue
        end
        if all([strcmp(s(9),'-'), strcmp(s(12),'h'), strcmp(s(15),'m'), ...
                strcmp(s(18),'s')])
            dt_dir = datetime(s(1:18), 'InputFormat', ...
                'uuuuMMdd-HH''h''mm''m''ss''s', 'Format', ...
                'uuuu-MM-dd''T''HH:mm:ss.SSSSSS','TimeZone','UTC');
            % skip reading the subdirectory if the report is not between
            % dt_begin and dt_end
            if ~isempty(dt_end) && dt_dir > dt_end
                continue
            end
            % loop through the surfacing subdirctory
            [files, lndex] = allfile([sdirs{jj} '/']);
            for kk = 1:lndex
                % check whether the file is the reported seismogram
                % The title begins with the datetime of the record.
                f = removepath(files{kk});
                s = split(f, '.');
                try
                    dt = datetime(f(1:15), 'Format', ...
                        'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', 'TimeZone', 'UTC');
                catch
                    continue
                end
                % only gather records within [dt_begin, dt_end]
                if all([strcmp(s{end}, format), ...
                        isempty(dt_begin) || dt >= dt_begin, ...
                        isempty(dt_end) || dt <= dt_end])
                    % check for tipe
                    if isempty(tipe)
                        fndex = fndex + 1;
                        allfiles{fndex} = files{kk};
                    else
                        if strcmp(s{end-2}, tipe)
                            fndex = fndex + 1;
                            allfiles{fndex} = files{kk};
                        end
                    end
                end
            end
        end
    end
end
end