function [x, z] = bathymetryprofile(length, npts, lonlat, az)
% [x, z] = bathymetryprofile(width, npts, lonlat, az)
%
% Generate a bathymetry profile centered at [LON LAT] along a straight path
% in a specified azimuthal direction
%
% INPUT:
% length        length of the path
% npts          number of points on the path
% [lon lat]     longitude and latitude of the center of the path
% az            azimuthal direction of the path
%
% OUTPUT:
% x             x-coordinate along the path from 0 to LENGTH
% z             elevation
% 
% Last modified by sirawich-at-princeton.edu, 03/17/2022

% Earth's radius in meter
R = 6371000;

% convert meters to degrees
m2deg = 180 / pi / R;

% x-coordinate of the bathymetry
x = linspace(-length/2, length/2, npts)';

% latitude and longitude of the bathymetry
[lats, lons] = reckon(lonlat(2), lonlat(1), x * m2deg, az);

% request grided elevation
[longrid, latgrid, elev] = bathymetry([], ...
    [min(lons) max(lons)] + [-0.1 0.1], ...
    [min(lats) max(lats)] + [-0.1 0.1], ...
    false);

% convert longitude to [-180 180]
longrid = mod(longrid + 180, 360) - 180;

% get the elevation along the straight path
z = interp2(latgrid, longrid, double(elev), lats, lons);

% shift x-cooridnate to [0 width]
x = x - x(1);
end