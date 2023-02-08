function [cmtp,T,P,v] = cmtpolarity(M, d, az, p, mod, is_P)
% [cmtp,T,P,v] = CMTPOLARITY(M, d, az, p, mod, is_P)
%
% Compute the polarity (first motion) of the P-wave from the moment tensor
% and the take-off direction.
%
% INPUT:
% M         moment tensor M=[Mrr Mtt Mpp Mrt Mrp Mtp]
% d         centroid depth in km
% az        azimuth of the station from the epicenter
% p         ray parameter
% mod       1D-Earth model [default: 'ak135']
% is_P      whether the path going downward from the hypocenter
%
% OUTPUT:
% cmtp      polarity
%           +1 for tension axis
%           -1 for compressional axis
%            0 for nodal planes
% T         tension axis
% P         compressional axis
% v         take-off direction
%
% Last modified by sirawich-at-princeton.edu: 02/08/2023

defval('mod', 'ak135')

% convert to 3x3 symmetric matrix
MM = [M(1) M(4) M(5); M(4) M(2) M(6); M(5) M(6) M(3)];

% compute the principal axes
[V, D] = eig(MM);

% T-axis
[~,iT] = max(diag(D));
T = V(:,iT);

% P-axis
[~,iP] = min(diag(D));
P = V(:,iP);

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
phi = asin(p * vp / r);

% determine "take-off" unit vector
az_rad = az * pi / 180;
if is_P
    v = [-cos(phi); -sin(phi) * cos(az_rad); sin(phi) * sin(az_rad)];
else
    v = [ cos(phi); -sin(phi) * cos(az_rad); sin(phi) * sin(az_rad)];
end

% calculate the polarity
cmtp = (T' * v) ^ 2 - (P' * v) ^ 2;
end