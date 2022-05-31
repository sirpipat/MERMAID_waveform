function loadbackup(n)
% LOADBACKUP(n)
%
% Loads the n-th backup files from $MFILES/backup/ to the current working
% directory. WARNING: It overwrites the existing variables in the working
% directory.
%
% INPUT:
% n         1 - most recent         [default]
%           2 - second most recent
%           3 - least recent
%
% SEE ALSO:
% BACKUP
% 
% Last modified by sirawich-at-princeton.edu, 05/31/2022

defval('n', 1)
fname = sprintf('%sbackup/backup_%02d.mat', getenv('MFILES'), n);

% display the loaded variable tothe screen
showloadedvars(fname);

% actual loads to the working space
evalin('base', sprintf('load(''%s'')', fname));
end

% load all variables from a MAT file and display on the screen
function showloadedvars(fname)
    fprintf('Loading %s\n\n', fname);
    fprintf('List of loaded variables\n')
    load(fname);
    whos
end