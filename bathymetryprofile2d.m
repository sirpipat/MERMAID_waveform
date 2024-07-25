function [ll, tt, zz, llons, llats] = bathymetryprofile2d(size, npts, lonlat, az)
% [ll, tt, zz, llons, llats] = ...
%     BATHYMETRYPROFILE2D([length width], [npts_l npts_t], [lon lat], az)
%
% Generate a 2D bathymetry profile centered at [LON LAT] along a straight 
% path in a specified azimuthal direction (longitudinal). The azimuth will 
% be the direction from the left to the right of the profile on the 
% longitudinal direction. The transverse direction is AZ - PI/2.
%
% INPUT
% [length width]        size on longitudinal and transverse direction
% [npts_l npts_t]       number of grid points on longitudinal and 
%                       transverse direction
% [lon lat]             longitude and latitude of the center of the profile
% az                    azimuthal direction of the path
% 
% OUTPUT
% ll                    longitudinal (azimuthal) coordinate of the grid
%                       locations
% tt                    transverse coordinate of the grid locations
% zz                    elevation at the grid locations
% llons                 longitude at the grid locations
% llats                 latitude  at the grid locations
%
% SEE ALSO
% BATHYMETRY, BATHYMETRYPROFILE
%
% Last modified by sirawich-at-princeton.edu, 07/01/2024

% Earth's radius in meter
R = 6371000;

% convert meters to degrees
m2deg = 180 / pi / R;

%% computation
% STEP 1: Compute the central points of the 1D profiles along the
% longitudinal direction
% coordinates of the bathymetry
y = linspace(-size(2)/2, size(2)/2, npts(2))';
x = linspace(-size(1)/2, size(1)/2, npts(1))';
[ll, tt] = meshgrid(x, y);
tt = flip(tt, 1);

% latitude and longitude of the bathymetry
[lats_center, lons_center] = reckon(lonlat(2), lonlat(1), y * m2deg, az-90);

llons = nan(npts(2), npts(1));
llats = nan(npts(2), npts(1));

for ii = 1:npts(2)
    [lats, lons] = reckon(lats_center(ii), lons_center(ii), x * m2deg, az);
    llons(ii, :) = lons';
    llats(ii, :) = lats';
end
llons = flip(llons);
llats = flip(llats);

% This handles cases where lons span across 180 E/W longitude
% Determine if there is any jump in lons
difflons = llons(:,2:end) - llons(:,1:(end-1));
if ~or(all(difflons < 0, 'all'), all(difflons > 0, 'all'))
    llons = mod(llons, 360);
    is_across180 = true;
else
    is_across180 = false;
end

% request grided elevation
[longrid, latgrid, elev] = bathymetry([], ...
    [min(llons, [], 'all') max(llons, [], 'all')] + [-0.1 0.1], ...
    [min(llats, [], 'all') max(llats, [], 'all')] + [-0.1 0.1], ...
    false);

% convert longitude to [-180 180] when lons do not span across 180 E/W
if ~is_across180
    longrid = mod(longrid + 180, 360) - 180;
end

% handle edge cases where longrid span across 180 E/W longitude
if longrid(1) > longrid(end)
    longrid = mod(longrid, 360);
    llons = mod(llons, 360);
end

% get the elevation along the straight path
zz = interp2(latgrid, longrid, double(elev), llats, llons);
end