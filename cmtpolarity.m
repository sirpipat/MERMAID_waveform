function [fp,fsv,fsh,v] = cmtpolarity(M, d, az, p, mod, down)
% [fp,fsv,fsh,v] = CMTPOLARITY(M, d, az, p, mod, is_P)
%
% Compute the polarity (first motion) of the P-wave from the moment tensor
% and the take-off direction.
%
% see Dahlen and Tromp, Theoretical Global Seismology, 1998 page 529
%
% INPUT:
% M         moment tensor M=[Mrr Mtt Mpp Mrt Mrp Mtp]
% d         centroid depth in km
% az        azimuth of the station from the epicenter-\
% p         ray parameter
% mod       1D-Earth model [default: 'ak135']
% down      whether the path going downward from the hypocenter
%
% OUTPUT:
% fp        polarization of P-wave
%           +1 for tension axis
%           -1 for compressional axis
%            0 for nodal planes
% fsv       polarization of SV-wave
% fsh       polarization of SH-wave
% v         take-off direction
%
% Last modified by sirawich-at-princeton.edu: 03/01/2023

defval('mod', 'ak135')

% convert to 3x3 symmetric matrix
MM = [M(1) M(4) M(5); M(4) M(2) M(6); M(5) M(6) M(3)];

% determine P-wave speed at the event focus
if strcmpi(mod, 'ak135')
    fname = fullfile(getenv('IFILES'), 'EARTHMODELS', 'PHYSICAL', ...
        'ak135.mat');
else
    fprintf('To do: including other models than ak135.\n');
    return
end
load(fname, 'model');
% add a very small value to duplicated depth to have all depths unique
for ii = 1:(length(model.d) - 1)
    if model.d(ii) >= model.d(ii+1)
        model.d(ii+1) = model.d(ii) + 1e-8;
    end
end
vp = interp1(model.d , model.vp, d);

% determine "take-off" angle
R = 6371;       % Earth radius in km
r = R - d;
sintheta = p * vp / r;
if abs(sintheta) > 1
    fprintf('Error: take-off angle is not physical.\n')
    fp = nan;
    return
end
theta = asin(p * vp / r);

% determine "take-off" unit vector
az_rad = az * pi / 180;
if down
    v = [-cos(theta); -sin(theta) * cos(az_rad); sin(theta) * sin(az_rad)];
else
    v = [ cos(theta); -sin(theta) * cos(az_rad); sin(theta) * sin(az_rad)];
end

% polarization vectors
np = v;
nsh = cross(v, [1 0 0]'); nsh = nsh / norm(nsh, 2);
nsv = cross(nsh, np);

% compute radiation pattern
% Dahlen and Tromp, Theoretical Global Seismology, 1998 page 529 Eq. 12.262
M0 = sqrt(sum(sum(MM .* MM)) / 2);
fp = v' * (MM / M0) * np;
fsv = 1/2 * (v' * (MM / M0) * nsv + nsv' * (MM / M0) * v);
fsh = 1/2 * (v' * (MM / M0) * nsh + nsh' * (MM / M0) * v);
end