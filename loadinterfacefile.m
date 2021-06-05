function [itfs, layers] = loadinterfacefile(fname)
% [itfs, layers] = LOADINTERFACEFILE(fname)
%
% Reads interfaces and layer thicknesses an interface file.
%
% INPUT:
% fname         name of the interface file
%
% OUTPUT:
% itfs          interfaces, an array of struct with following fields
%                   npts    the number of points
%                   pts     list of [x, z] coordinates in meters
% layers        number of vertical spectral elements for each layer
%
% SEE ALSO:
% writeinterfacefile
%
% Last modified by Sirawich Pipatprathanporn, 06/05/2021

% open the file
fid = fopen(fname, 'r');

line = fgetl(fid);
% skip the comments
while strcmp(line(1), '#')
    line = fgetl(fid);
end

% the number of interfaces
n = sscanf(line, '%d', 1);

% storage
itfs = cell(n, 1);
layers = zeros(n-1, 1);

for ii = 1:n
    line = fgetl(fid);
    % skip the comments
    while strcmp(line(1), '#')
        line = fgetl(fid);
    end
    % the number of interface points
    itf.npts = sscanf(line, '%d', 1);
    % load the interface points
    itf.pts = fscanf(fid, '%f', [2 itf.npts])';
    
    itfs{ii} = itf;
    % skip the newline character at the end
    fgetl(fid);
end

for ii = 1:(n-1)
    line = fgetl(fid);
    % skip the comments
    while strcmp(line(1), '#')
        line = fgetl(fid);
    end
    layers(ii) = sscanf(line, '%d', 1);
end

% close the file
fclose(fid);
end