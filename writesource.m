function writesource(sources, fname)
% WRITESOURCE(sources, fname)
%
% Writes a source file.
%
% INPUT
% sources       sources
% fname         name of the source file
%
% Last modified by Sirawich Pipatprathanporn, 06/08/2021

num_sources = length(sources);

if isempty(fname)
    fid = 1;
else
    fid = fopen(fname, 'w');
end

for ii = 1:num_sources
    source = sources{ii};
    writecomment(fid, sprintf(['# source %d.  The components of a ' ...
        'moment tensor source must be given in N.m, not in dyne.cm as ' ...
        'in the DATA/CMTSOLUTION source file of the 3D version of the ' ...
        'code.'], ii));
    writebool(fid, 'source_surf', source.source_surf, ...
        'source inside the medium or at the surface');
    writefloat(fid, 'xs', source.xs, 'source location x in meters');
    writefloat(fid, 'zs', source.zs, 'source location z in meters');
    writeint(fid, 'source_type', source.source_type, ...
        'elastic force or acoustic pressure = 1 or moment tensor = 2');
    writeint(fid, 'time_function_type', source.time_function_type, ...
        ['Ricker = 1, first derivative = 2, Gaussian = 3, Dirac = 4, ' ...
        'Heaviside = 5']);
    writecomment(fid, '# time function_type == 8 source read from file, if time function_type == 9 : burst');
    writecomment(fid, '# If time_function_type == 8, enter below the custom source file to read (two columns file with time and amplitude) :');
    writecomment(fid, "# (For the moment dt must be equal to the dt of the simulation. File name can't exceed 150 characters)");
    writestring(fid, 'name_of_source_file', source.name_of_source_file, ...
        'Only for option 8 : file containing the source wavelet');
    writefloat(fid, 'burst_band_width', source.burst_band_width, ...
        'Only for option 9 : band width of the burst');
    writefloat(fid, 'f0', source.f0, ...
        'dominant source frequency (Hz) if not Dirac or Heaviside');
    writefloat(fid, 'tshift', source.tshift, ...
        'time shift when multi sources (if one source, must be zero)');
    writefloat(fid, 'anglesource', source.anglesource, ...
        'angle of the source (for a force only)');
    writefloat(fid, 'Mxx', source.Mxx, ...
        'Mxx component (for a moment tensor source only)');
    writefloat(fid, 'Mzz', source.Mzz, ...
        'Mzz component (for a moment tensor source only)');
    writefloat(fid, 'Mxz', source.Mxz, ...
        'Mxz component (for a moment tensor source only)');
    writefloat(fid, 'factor', source.factor, ...
        'amplification factor');
    writeblank(fid);
end

% close the file
if fid >= 3
    fclose(fid);
end
end

function writeblank(fid)
fprintf(fid, '\n');
end

function writecomment(fid, comment)
if isempty(comment)
    writeblank(fid);
else
    comment = char(comment);
    if ~strcmp(comment(1), '#')
        fprintf(fid, '# %s\n', comment);
    else
        fprintf(fid, '%s\n', comment);
    end
end
end

function writebool(fid, name, value, comment)
if value == 0
    var_string = '.false.';
else
    var_string = '.true.';
end
writestring(fid, name, var_string, comment);
end


function writefloat(fid, name, value, comment)
if value == 0
    var_string = '0';
elseif and(abs(value) >= 1, abs(value) < 100)
    var_string = sprintf('%.3f', value);
else
    var_string = sprintf('%.3e', value);
    var_string = replace(var_string, 'e', 'd');
end
writestring(fid, name, var_string, comment);
end

function writeint(fid, name, value, comment)
if isempty(comment)
    fprintf(fid, '%-31s = %d\n', name, value);
else
    fprintf(fid, '%-31s = %-14d # %s\n', name, value, comment);
end
end

function writestring(fid, name, value, comment)
if isempty(comment)
    fprintf(fid, '%-31s = %s\n', name, value);
else
    fprintf(fid, '%-31s = %-14s # %s\n', name, value, comment);
end
end