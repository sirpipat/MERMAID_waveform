function s = listvars(fname)
% s = LISTVARS(fname)
%
% Lists variables in a MAT file. Best used for learning file content before
% loading the variables to workspace which may overwrite the existing
% variables.
%
% LISTVARS(fname) will output the list of variables to the console (command
% window) unless the output variable s is specified.
%
% INPUT:
% fname         full filename of the MAT file
%
% OUTPUT:
% s             struct listing the variables in the file
%
% SEE ALSO:
% WHOS, MATFILE, LOAD
%
% Last modified by sirawich-at-princeton.edu, 09/21/2023

fprintf('Loading %s\n\n', fname);
fprintf('List of variables\n')
load(fname);
clear fname       % do not show variable fname on the list
if nargout == 0
    whos
else
    s = whos;
end
end