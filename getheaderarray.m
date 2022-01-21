function HdrData = getheaderarray(sacfiles)
% HdrData = getheaderarray(sacfiles)
%
% Retrieves SAC headers from multiple SAC files and structures into a
% struct of arrays grouping the value of the same variable together.
%
% INPUT:
% sacfiles      cell array containing full path SAC files
%
% OUTPUT:
% HdrData       the header structure array
% 
% SEE ALSO:
% READSAC
%
% Last modified by sirawich-at-princeton.edu, 01/20/2022

N = length(sacfiles);

% get the header fields to construct the header array
[~, hdr] = readsac(sacfiles{1});
hdrfields = fieldnames(hdr);

% construct the header array
cmd = ['HdrData = struct(''' hdrfields{1} ''', []'];
for ii = 2:length(hdrfields)
    cmd = [cmd ', ''' hdrfields{ii} ''', []'];
end
cmd = [cmd ');'];

% evaluate 'HdrData = struct(hdrfields{1}, [], ..., hdrfields{end}, []);'
eval(cmd);

% read SAC files and sort into struct of arrays
for ii = 1:N
    [~, hdr] = readsac(sacfiles{ii});
    hdrfields = fieldnames(hdr);
    for jj = 1:length(hdrfields)
        eval(sprintf('fieldvar = hdr.%s;', hdrfields{jj}));
        if ischar(fieldvar)
            eval(sprintf('HdrData.%s{%d,1} = ''%s'';', hdrfields{jj}, ...
                ii, fieldvar));
        else
            eval(sprintf('HdrData.%s(%d,1) = %f;', hdrfields{jj}, ii, ...
                fieldvar));
        end
    end
end
end