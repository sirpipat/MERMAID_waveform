function writeinterfacefile(itfs, layers, fname)
% WRITEINTERFACEFILE(itfs, layers, fname)
%
% Writes interfaces and layers to an interface file.
%
% INPUT:
% itfs          interfaces, an array of struct with following fields
%                   npts    the number of points
%                   pts     list of [x, z] coordinates in meters
% layers        number of vertical spectral elements for each layer
% fname         name of the interface file you want to create
%
% SEE ALSO:
% LOADINTERFACEFILE
%
% Last modified by Sirawich Pipatprathanporn, 07/22/2021

n = length(itfs);

if isempty(fname)
    fid = 1;
else
    fid = fopen(fname, 'w');
end

fprintf(fid, '#\n');
fprintf(fid, '# number of interfaces\n');
fprintf(fid, '#\n');
fprintf(fid, ' %d\n', n);
fprintf(fid, strcat('# for each interface below, we give the number of', ...
    ' points and then x,z for each point\n'));
fprintf(fid, '#\n');

for ii = 1:n
    if ii == 1
        str = ' (bottom of the mesh)';
    elseif ii == n
        str = ' (top of the mesh)';
    else
        str = '';
    end
    fprintf(fid, '#\n');
    fprintf(fid, '# interface number %d%s\n', ii, str);
    fprintf(fid, '#\n');
    fprintf(fid, ' %d\n', itfs{ii}.npts);
    fprintf(fid, ' %6d %6d\n', round(itfs{ii}.pts'));
end

fprintf(fid, '#\n');
fprintf(fid, strcat('# for each layer, we give the number of', ...
    ' spectral elements in the vertical direction\n'));
fprintf(fid, '#\n');

for ii = 1:length(layers)
    if ii == 1
        str = ' (bottom layer)';
    elseif ii == length(layers)
        str = ' (top layer)';
    else
        str = '';
    end
    fprintf(fid, '#\n');
    fprintf(fid, '# layer number %d%s\n', ii, str);
    fprintf(fid, '#\n');
    fprintf(fid, ' %d\n', layers(ii));
end

if fid >= 3
    fclose(fid);
end
end