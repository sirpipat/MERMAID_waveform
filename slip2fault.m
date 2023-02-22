function [strike, dip, rake] = slip2fault(slip, normal, coord)
% [strike, dip, rake] = SLIP2FAULT(slip, normal, coord)
%
% Computes strike, dip, and rake from the slip vector and the normal vector
% to the fault plane given the normal vector is pointing upward (i.e. away
% from Earth's center).
%
% INPUT:
% slip          the slip vector
% normal        the normal vector to the fault plane
%               counterclockwise from the strike
% coord         coordinate systems of the output vectors
%       'xyz' : north-east-down
%       'rtp' : up-south-east  (radius-theta-phi as in spherical coord)
%
% OUTPUT:
% strike        azimuth of the strike in degrees from north
% dip           dip angle in degrees from the horizontal plane
% rake          angle of the slip on the fault plane in degrees
%
% SEE ALSO:
% FAULT2SLIP
%
% Last modifed by sirawich-at-princeton.edu: 02/22/2023


% normalized the input
slip = slip / norm(slip, 2);
normal = normal / norm(normal, 2);

% convert to rtp coordinate
if strcmpi(coord, 'xyz')
    slip = [-slip(3); -slip(1); slip(2)];
    normal = [-normal(3); -normal(1); normal(2)];
elseif strcmpi(coord, 'rtp')
else
    error("coord must be either 'xyz' or 'rtp'.");
end

% flip the slip and normal if the normal points downward
if normal(1) < 0
    normal = -normal;
    slip = -slip;
end

% compute strike and dip from normal vector to the fault plane
phi = pi/2 - atan2(normal(3), normal(2));
delta = acos(normal(1));

% construct two orthonormal basis vectors on the fault plane
% == > ex - along the strike
% == > ey - up the fault plane
% == > ez - normal to the plane (not used here)
ex = [0; -cos(phi); sin(phi)];
ey = cross(normal, ex);

% decompose the slip vector to the two orthonormal basis
% slip = x * ex + y * ey where x,y are scalars
x = dot(slip, ex);
y = dot(slip, ey);

% compute the rake
lambda = atan2(y, x);

% convert radian to degrees
rad2deg = 180 / pi;
strike = mod(phi * rad2deg, 360);
dip = delta * rad2deg;
rake = mod(lambda * rad2deg, 360);
end