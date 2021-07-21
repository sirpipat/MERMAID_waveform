function s = read_sources(source_file)
% s = READ_SOURCES(source_file)
%
% Reads sources from source file.
% OBSOLETED: use LOADSOURCE instead which can handle multiple sources,
% more source variables, and flexible comments.
%
% INPUT
% source_file   full filename of a source file
%
% OUTPUT
% s             sources, struct with a following fields
%               - x
%               - z
%               - st        source type
%               - tt        time function type
%               - f0        dominant frequency
%               - tshift    time shift
%               - angle     angle source
%               - M         moment tensor
%               - A         amplification factor
%
% SEE ALSO
% LOADSOURCE, WRITESOURCE
%
% Last modified by Siraich Pipatprathanporn, 07/13/2021

fid = fopen(source_file);

% skip header and source_surf
for ii = 1:2
    fgetl(fid);
end

% read positions
line = fgetl(fid);
x = sscanf(line, 'xs = %f');
line = fgetl(fid);
z = sscanf(line, 'zs = %f');

% read source type
line = fgetl(fid);
st = sscanf(line, 'source_type = %d');

% read time function type
line = fgetl(fid);
tt = sscanf(line, 'time_function_type = %d');

for ii = 1:5
    fgetl(fid);
end

% read dominant source frequency
line = fgetl(fid);
f0 = sscanf(line, 'f0 = %f');

% read time shift
line = fgetl(fid);
tshift = sscanf(line, 'tshift = %f');

% read angle source
line = fgetl(fid);
angle = sscanf(line, 'anglesource = %f');

% read moment tensor
line = fgetl(fid);
Mxx = sscanf(line, 'Mxx = %f');
line = fgetl(fid);
Mzz = sscanf(line, 'Mzz = %f');
line = fgetl(fid);
Mxz = sscanf(line, 'Mxz = %f');
M = [Mxx Mxz; Mxz Mzz];

% read amplification factor
line = fgetl(fid);
factor = sscanf(line, 'factor = %f');

s = struct('x', x, 'z', z, 'st', st, 'tt', tt, 'f0', f0, ...
    'tshift', tshift, 'angle', angle, 'M', M, 'factor', factor);

fclose(fid);
end