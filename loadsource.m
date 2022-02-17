function sources = loadsource(fname)
% sources = LOADSOURCE(fname)
%
% Loads a source file. 
%
% INPUT
% fname         name of the source file
%
% OUTPUT
% sources       sources, array of struct(s) with following fields
%               - source_surf
%               - xs
%               - zs
%               - source_type
%               - time_function_type
%               - name_of_source_file
%               - burst_band_width
%               - f0
%               - tshift
%               - anglesource
%               - Mxx
%               - Mzz
%               - Mxz
%               - factor
%               - vx         (set to zero for 'master' branch)
%               - vz         (set to zero for 'master' branch)
%
% Last modified by Sirawich Pipatprathanporn, 02/16/2022

sources = {};
num_sources = 0;

fid = fopen(fname, 'r');

line = fgetl(fid);
while ischar(line)
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.source_surf = readbool(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.xs = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.zs = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.source_type = readint(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.time_function_type = readint(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.name_of_source_file = readstring(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.burst_band_width = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.f0 = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.tshift = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.anglesource = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.Mxx =  readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.Mzz = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.Mxz  = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    source.factor = readfloat(line);
    
    % skip comments / headers
    line = fgetl(fid);
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    % try to read vx and vz. If they are not there, move on to the next.
    if strcmp(line(1:2), 'vx')
        source.vx = readfloat(line);
        
        % skip comments / headers
        line = fgetl(fid);
        while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
            line = fgetl(fid);
            if ~ischar(line)
                break;
            end
        end
        if ~ischar(line)
            break;
        end
        
        source.vz = readfloat(line);
        
        line = fgetl(fid);
    else
        source.vx = 0.0;
        source.vz = 0.0;
    end
    
    num_sources = num_sources + 1;
    sources{num_sources} = source;
    
    
end

fclose(fid);
end

function value = readstring(line)
% find equal sign
where_start = strfind(line, '=');

% find # where the comment starts
where_end = strfind(line, '#');
if isempty(where_end)
    where_end = length(line) + 1;
end

% read the value
value = strip(sscanf(line((where_start+1):(where_end-1)), '%c'));
end

function value = readbool(line)
value = readstring(line);
if strcmp(value, '.true.')
    value = true;
elseif strcmp(value, '.false.')
    value = false;
else
    % do not know what to do
    error(strcat('ValueError: cannot read a boolean\n', ...
                 sprintf('line >> %s\n', line)));
end
end

function value = readint(line)
% find equal sign
where = strfind(line, '=');
% read the value
value = sscanf(line((where+1):end), '%d', 1);
end

function value =  readfloat(line)
% find equal sign
where = strfind(line, '=');
% change the exponent notation syntax from 'd' to 'e'
line = replace(line, 'd', 'e');
% read the value
value = sscanf(line((where+1):end), '%f', 1);
end