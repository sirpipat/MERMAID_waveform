function [allfiles, fndex] = allfilen(ddir, n)
% [allfiles,fndex] = ALLFILEN(ddir, n)
%
% Makes cell array with the names in a directory and subdirectories up to N
% levels. For N=1, the output is the same as ALLFILE(ddir).
%
% INPUT:
%
% ddir     Where you keep the data (trailing slash needed!)
% n        the number of level in the directory to look for
%          if n=1, it will look for only contents in the directory but 
%                  any content in the subdirectories
%          if n=2, it will look for contents in the directory and the
%                  subdirectories but ignore any content in the
%                  sub-subdirectories.
%
% OUTPUT:
%
% allfiles   Bottom-file list with complete file names
% fndex      The total number of elements in the list
%            -1 if dir does not exist as a directory
%             0 if dir is an empty directory
%
% SEE ALSO:
% ALLFILE, LS2CELL, DIR
%
% Last modified by Sirawich Pipatprathanporn, 11/01/2021

% make a table of contents
allfiles = {};
fndex = 0;

% check input argument
if ~isint(n) || n < 1
    fprintf('n must be a positive integer\n');
    return
end

% add a trailing slash to the directory if there is not one
if ~strcmp(ddir(end), '/')
    ddir = strcat(ddir, '/');
end

% look for the content in the directory
[nodes, nndex] = allfile(ddir);

% exit if ddir is not a non-empty directory
if nndex <= 0
    fndex = nndex;
    return
end

% base case: n = 1
if n == 1
   allfiles = nodes;
   fndex = nndex;
% recursion
else
    for ii = 1:nndex
        [cls, cndex] = allfilen(nodes{ii}, n-1);
        % it is a file: append the file to the list
        if cndex == -1
            allfiles{fndex + 1} = nodes{ii};
            fndex = fndex + 1;
        % it is an empty directory: skip
        elseif cndex == 0
            continue
        % it is not empty directory
        else
            for jj = 1:cndex
                allfiles{fndex + jj} = cls{jj};
            end
            fndex = fndex + cndex;
        end
    end
end
end