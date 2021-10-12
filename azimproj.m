function [azDeg, distDeg] = azimproj(center, pts)
% [azDeg, distDeg] = AZIMPROJ(center, pts)
% 
% Computes the azimuths and distances of multiple points from the center of
% an azimuthal map projection.
%
% INPUT:
% center        [lon lat] of the center of the azimuthal map projection
% pts           [lons lats] of the points to be projected
%               lons and lats may be column vectors of equal length
%
% OUTPUT:
% azDeg         azimuths in degrees of the projected points
% distDeg       distances in degrees of the projected points
%
% SEE ALSO:
% AZIM, GRCDIST, AZIMDIST
%
% Last modified by sirawich-at-princeton.edu, 10/12/2021

azDeg = zeros(size(pts,1), 1);
distDeg = zeros(size(pts,1), 1);

% find azimuth and distance from the center of map projection
for ii = 1:size(pts, 1)
    [azDeg(ii, 1), ~, distDeg(ii, 1)] = azimdist(center, pts(ii, :));
end
end