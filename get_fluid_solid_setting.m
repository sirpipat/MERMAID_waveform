function [d_fluid, d_solid] = get_fluid_solid_setting(interface_file)
% [d_fluid, d_solid] = GET_FLUID_SOLID_SETTING(interface_file)
%
% Reads an interface file to get the thicknesses of fluid and solid layers.
%
% INPUT
% interface_file        full filename of the interface file
%
% OUTPUT
% d_fluid               thickness of the fluid layer
% d_solid               thickness of the solid layer
%
% Last modified by Sirawich Pipatprathanporn, 02/25/2021

fid = fopen(interface_file, 'r');
% skip to interface number 2: fluid-solid
for ii = 1:16
    fgetl(fid);
end

% read interface number 2
npts2 = fscanf(fid, '%d', 1);
sizeData = [2 npts2];
interface2 = fscanf(fid, '%f %f', sizeData);
seafloor = interface2';

% skip to interface number 3: sea level
for ii = 1:4
    fgetl(fid);
end

% read interface number 2
npts3 = fscanf(fid, '%d', 1);
sizeData = [2 npts3];
interface3 = fscanf(fid, '%f %f', sizeData);
seasurface = interface3';

fclose(fid);

% calculates thickness of the two layers
d_fluid = mean(seasurface(:,2)) - mean(seafloor(:,2));
d_solid = mean(seafloor(:,2));
end