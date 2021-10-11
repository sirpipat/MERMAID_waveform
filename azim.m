function [azDeg, bazDeg] = azim(lon1lat1, lon2lat2)
% [azDeg, bazDeg] = AZIM([lon1 lat1], [lon2 lat2])
% Calculates the azimuth and back azimuth of two points.
%
% INPUT:
% [lon1 lat1]       longitude and latitude of the starting point (degrees)
% [lon2 lat2]       longitude and latitude of the ending point (degrees)
%
% OUTPUT:
% azDeg             azimuth (in degrees) of the ending point direction from
%                   the starting point
% bazDeg            (back) azimuth (in degrees) of the starting point 
%                   direction from the ending point
%
% SEE ALSO:
% GRCDIST, AZIMDIST
%
% Last modified by sirawich-at-princeton.edu, 10/11/2021

% Conversion to radians
lon1lat1=lon1lat1 * pi / 180;
lon2lat2=lon2lat2 * pi / 180;

% longitude difference between the two points
londiff = mod(lon2lat2(1) - lon1lat1(1), 2*pi);

% check whether the ending point is east of the starting point
% if so swap the positions
if londiff > pi
    temp = lon1lat1;
    lon1lat1 = lon2lat2;
    lon2lat2 = temp;
    londiff = 2 * pi - londiff;
    swap = true;
else
    swap = false;
end

% distance between the two points using Rule of Cosine
dist = acos(sin(lon1lat1(2)) * sin(lon2lat2(2)) + ...
    cos(lon1lat1(2)) * cos(lon2lat2(2)) * cos(londiff));

% compute azimuth and back azimuth using Rule of Cosine
az = acos((sin(lon2lat2(2)) - sin(lon1lat1(2)) * cos(dist)) / ...
    (cos(lon1lat1(2)) * sin(dist)));
baz = acos((sin(lon1lat1(2)) - sin(lon2lat2(2)) * cos(dist)) / ...
    (cos(lon2lat2(2)) * sin(dist)));

% remove any imaginary components caused by az,baz = acos[+-(1 + \epsilon)] 
% where \epsilon is from catastrophic cancellation.
az = real(az);
baz = real(baz);

% check whether the ending point is east of the starting point
% if so, swap back
if swap
    temp = az;
    az  = mod(2*pi - baz, 2*pi);
    baz = temp;
else
    baz = mod(2*pi - baz, 2*pi);
end

% convert back to degrees
azDeg   = az * 180 / pi;
bazDeg  = baz * 180 / pi;
end