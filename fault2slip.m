function [slip, normal] = fault2slip(strike, dip, rake, coord)
% [slip, normal] = FAULT2SLIP(strike, dip, rake, coord)
%
% Computes the unit vector of the slip direction and the normal vector to
% the fault plane
%
% see Aki&Richards, Quantitative Seismology, 2nd ed. page 108-109
%
% INPUT:
% strike        azimuth of the strike in degrees from north
% dip           dip angle in degrees from the horizontal plane
% rake          angle of the slip on the fault plane in degrees
%               counterclockwise from the strike
% coord         coordinate systems of the output vectors
%       'xyz' : north-east-down
%       'rtp' : up-south-east  (radius-theta-phi as in spherical coord)
%
% OUTPUT:
% slip          the slip vector
% normal        the normal vector to the fault plane
%
% SEE ALSO:
% SLIP2FAULT
%
% Last modified by sirawich-at-princeton.edu, 02/21/2023

% convert to radians
deg2rad = pi / 180;
lambda = rake * deg2rad;
delta = dip * deg2rad;
phi = strike * deg2rad;

% Aki&Richards (1980) Equation 4.88 and Figure 4.20
if strcmpi(coord, 'xyz')
    slip = [cos(lambda)*cos(phi) + sin(lambda)*cos(delta)*sin(phi); ...
            cos(lambda)*sin(phi) - sin(lambda)*cos(delta)*cos(phi); ...
            -sin(lambda)*sin(delta)];
    normal = [-sin(delta)*sin(phi); sin(delta)*cos(phi); -cos(delta)];
elseif strcmpi(coord, 'rtp')
    slip = [sin(lambda)*sin(delta); ...
            -cos(lambda)*cos(phi) - sin(lambda)*cos(delta)*sin(phi); ...
            cos(lambda)*sin(phi) - sin(lambda)*cos(delta)*cos(phi)];
    normal = [cos(delta); sin(delta)*sin(phi); sin(delta)*cos(phi)];
else
    error("coord must be either 'xyz' or 'rtp'.");
end
end