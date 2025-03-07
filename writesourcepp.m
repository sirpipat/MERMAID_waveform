function writesourcepp(sources, fname, branch)
% WRITESOURCEPP(sources, fname, branch)
%
% Writes a source file following SPECFEMPP yaml format.
%
% INPUT
% sources       sources
% fname         name of the source file
% branch        SPECFEM2D branch [Default: 'devel']
%               'devel'  (commit: cf89366717d9435985ba852ef1d41a10cee97884)
%
% Last modified by Sirawich Pipatprathanporn, 11/11/2024

defval('branch', 'devel')

num_sources = length(sources);

if isempty(fname)
    fid = 1;
else
    fid = fopen(fname, 'w');
end

writeint(fid, 'number-of-sources', num_sources, 0);
writestring(fid, 'sources', [], 0);
for ii = 1:num_sources
    source = sources{ii};
    switch source.source_type
        case 1
            writestring(fid, '- force', [], 1);
            writefloat(fid, 'x', source.xs, 2);
            writefloat(fid, 'z', source.zs, 2);
            writefloat(fid, 'angle', source.anglesource, 2);
            writestf(fid, source, 3);
        case 2
            writestring(fid, '- moment-tensor', [], 1);
            writefloat(fid, 'x', source.xs, 2);
            writefloat(fid, 'z', source.zs, 2);
            writefloat(fid, 'Mxx', source.Mxx, 2);
            writefloat(fid, 'Mzz', source.Mzz, 2);
            writefloat(fid, 'Mxz', source.Mxz, 2);
            writefloat(fid, 'angle', source.anglesource, 2);
            writestf(fid, source, 2);
        otherwise
            writestring(fid, '- adjoint-source', [], 1);
            writestring(fid, 'station_name', source.station_name, 2);
            writestring(fid, 'network_name', source.network_name, 2);
            writefloat(fid, 'x', source.xs, 2);
            writefloat(fid, 'z', source.zs, 2);
            writestf(fid, source, 2);
    end
end

% close the file
if fid >= 3
    fclose(fid);
end
end

function writeblank(fid)
fprintf(fid, '\n');
end

function writecomment(fid, indentlevel, comment)
defval('indentlevel', 0)
% white space for indentation is 2 (no tab please)
indent = repmat('    ', 1, indentlevel);
if isempty(comment)
    writeblank(fid);
else
    comment = char(comment);
    if ~strcmp(comment(1), '#')
        fprintf(fid, '%s# %s\n', indent, comment);
    else
        fprintf(fid, '%s%s\n', indent, comment);
    end
end
end

function writebool(fid, name, value, indentlevel, comment)
defval('indentlevel', 0)
defval('comment', [])
if value == 0
    var_string = 'False';
else
    var_string = 'True';
end
writestring(fid, name, var_string, indentlevel, comment);
end


function writefloat(fid, name, value, indentlevel, comment)
defval('indentlevel', 0)
defval('comment', [])
if value == 0
    var_string = '0';
elseif and(abs(value) >= 1, abs(value) < 100)
    var_string = sprintf('%.3f', value);
else
    var_string = sprintf('%.3e', value);
end
writestring(fid, name, var_string, indentlevel, comment);
end

function writeint(fid, name, value, indentlevel, comment)
defval('indentlevel', 0)
defval('comment', [])
% white space for indentation is 2 (no tab please)
indent = repmat('    ', 1, indentlevel);
if isempty(comment)
    fprintf(fid, '%s%s: %d\n', indent, name, value);
else
    fprintf(fid, '%s%s: %s # %d\n', indent, name, value, comment);
end
end

function writestring(fid, name, value, indentlevel, comment)
defval('value', [])
defval('indentlevel', 0)
defval('comment', [])
% white space for indentation is 2 (no tab please)
indent = repmat('    ', 1, indentlevel);
if isempty(value)
    printstr = sprintf('%s%s:', indent, name);
else
    printstr = sprintf('%s%s: %s', indent, name, value);
end
if isempty(comment)
    fprintf(fid, '%s\n', printstr);
else
    fprintf(fid, '%s # %s\n', printstr, comment);
end
end

function writestf(fid, source, indentlevel)
defval('indentlevel', 0)
switch source.time_function_type
    case 1
        writestring(fid, 'Ricker', [], indentlevel);
        writefloat(fid, 'factor', source.factor, indentlevel + 1);
        writefloat(fid, 'tshift', source.tshift, indentlevel + 1); 
        writefloat(fid, 'f0', source.f0, indentlevel + 1);
    case 4
        writestring(fid, 'Dirac', [], indentlevel);
        writefloat(fid, 'factor', source.factor, indentlevel + 1);
        writefloat(fid, 'tshift', source.tshift, indentlevel + 1);  
    case 7
        writestring(fid, 'External', [], indentlevel);
        writestring(fid, 'format', source.stfformat, indentlevel + 1);
        writestring(fid, 'stf', [], indentlevel + 1);
        writestring(fid, 'X-component', source.stffile.x, indentlevel + 2);
        writestring(fid, 'Y-component', source.stffile.y, indentlevel + 2);
        writestring(fid, 'Z-component', source.stffile.z, indentlevel + 2);
    otherwise
        error('Only implemented for Ricker (1), Dirac (4), or External (7) source time function.')
end
end