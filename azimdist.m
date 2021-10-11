function [azDeg, bazDeg, distDeg] = azimdist(lon1lat1, lon2lat2)
% [azDeg, bazDeg, distDeg] = AZIMDIST([lon1 lat1], [lon2 lat2])
%
% Calculates the azimuth, the back azimuth, and the distance between two 
% points on a great circle.
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
% distDeg           distance in degrees
%
% SEE ALSO:
% AZIM, GRCDIST
%
% Last modified by sirawich-at-princeton.edu, 10/11/2021

[~, distDeg] = grcdist(lon1lat1, lon2lat2);
[azDeg, bazDeg] = azim(lon1lat1, lon2lat2);
end

% EOF