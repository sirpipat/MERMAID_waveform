function [fp,fsv,fsh,v] = pPpolarity(M, evla, evlo, evdp, stla, stlo, model)
% [fp,fsv,fsh,v] = PPPOLARITY(M, evla, evlo, evdp, stla, stlo, model)
%
% Compute the polarity (first motion) of the pP-wave from the moment tensor
% and the take-off direction.
%
% see Dahlen and Tromp, Theoretical Global Seismology, 1998 page 529
%
% INPUT:
% M         moment tensor either 6 elements or 3x3 matrix
%           M=[Mrr Mtt Mpp Mrt Mrp Mtp]
%              / Mrr Mrt Mrp \
%           M=|  Mrt Mtt Mtp  |
%              \ Mrp Mtp Mpp /
% evla      event latitude
% evlo      event longitude
% evdp      event depth in km
% stla      station latitude
% stlo      station longitude
% model     Earth model [Default: 'ak135']
%
% OUTPUT:
% fp        polarization of pP-wave
%           +1 for tension axis
%           -1 for compressional axis
%            0 for nodal planes
% fsv       polarization of SV-wave
% fsh       polarization of SH-wave
% v         take-off direction
%
% Returns NaN if pP does not arrive at the station.
%
% Known issues: when the wave departs vertically, SV is not defined since
% all S-wave motions are horizontal. For that case, fsv = fsh = NaN.
%
% Last modified by sirawich-at-princeton.edu: 05/16/2024

defval('mod', 'ak135')

%% Part I: verify the inputs
% event is single
if or(all(size(evla) ~= [1 1]), all(size(evlo) ~= [1 1]))
    error('This function works for a single event only\n')
end

% source lat and lon are the same size
if or(ndims(stla) ~= ndims(stlo), size(stla) ~= size(stlo))
    error('STLA and STLO must have the same size.\n')
end

%% Part II: preprocessing
% convert to 3x3 symmetric matrix
if length(M) == 6
    MM = [M(1) M(4) M(5); M(4) M(2) M(6); M(5) M(6) M(3)];
else
    MM = M;
end

% convert longitude to [0 360]
stlo = mod(stlo, 360);
evlo = mod(evlo, 360);

%% Part III: determine paths from source-receivers
sname = sprintf('%s_%s.mat', mfilename, ...
    hash([evla evlo evdp reshape(stla, 1, []) reshape(stlo, 1, []) ...
    double(char(model))], 'SHA-1'));
pname = fullfile(getenv('IFILES'), 'HASHES', sname);

if ~exist(pname, 'file')
    [azim, dist, theta] = ...
        expensivefunction(evla, evlo, evdp, stla, stlo, model);
    
    % save
    fprintf('save the output to a file to %s ...\n', pname);
    save(pname, 'azim', 'dist', 'theta');
else
    % load
    fprintf('found the save in a file in %s\n', pname);
    fprintf('load the variables ...\n');
    load(pname, 'azim', 'dist', 'theta');
end

%% Part IV: calculate the polarity
v = nan([size(stlo) 3]);
fp = nan(size(stlo));
fsv = nan(size(stlo));
fsh = nan(size(stlo));
for ii = 1:size(stlo, 1)
    for jj = 1:size(stlo, 2)
        % determine "take-off" unit vector
        if isnan(theta(ii,jj))
            continue
        end
        azim_rad = azim(ii,jj) * pi / 180;
        theta_rad = theta(ii,jj) * pi / 180;
        % take off vector is always upward
        k = [ cos(theta_rad); -sin(theta_rad) * cos(azim_rad); ...
            sin(theta_rad) * sin(azim_rad)];
        % polarization unit vectors
        np = k;
        nsh = cross(k, [1 0 0]') / norm(cross(k, [1 0 0]'), 2);
        nsv = cross(nsh, np);

        % compute radiation pattern
        % Dahlen and Tromp, Theoretical Global Seismology, 1998 page 529 Eq. 12.262
        M0 = sqrt(sum(sum(MM .* MM)) / 2);
        v(ii,jj,:) = k;
        fp(ii,jj) = k' * (MM / M0) * np;
        fsv(ii,jj) = 1/2 * (k' * (MM / M0) * nsv + nsv' * (MM / M0) * k);
        fsh(ii,jj) = 1/2 * (k' * (MM / M0) * nsh + nsh' * (MM / M0) * k);
    end
end
end

% Gathers the expensive steps (especially TAUPTIME) to clearly show the
% expensive computations' inputs and outputs for saving/loading variables
% to speed up recurring function calls. This function is meant to be called
% by CMTPOLARITY only (for now). If this function is also called by other
% functions, please consider making it a stand-alone function.
%
% INPUT:
% evla          event latitude
% evlo          event longitude
% evdp          event depth (in km)
% stla          station latitude(s)
% stlo          station longitude(s)
% model         1D Earth model name either 'prem', 'ak135', or 'iasp91'
%
% OUTPUT:
% azim          azimuth at the event heading to the station
% dist          epicentral distance from the event
% theta         take-off angle from the vertical [0 <= theta <= 90 degrees]
function [azim, dist, theta] = ...
    expensivefunction(evla, evlo, evdp, stla, stlo, model)
% CONSTANTS
R_EARTH = 6371; % km

% compute the azimuth and epicentral distance
azim = nan(size(stlo));
dist = nan(size(stlo));
theta = nan(size(stlo));
for ii = 1:numel(stlo)
    [azim(ii), ~, dist(ii)] = azimdist([evlo evla], [stlo(ii) stla(ii)]);
    
    if isnan(azim(ii))
        azim(ii) = 0;
    end

    % compte the ray parameter of pP wave
    tt = tauptime('mod', model, 'depth', evdp, ...
        'phases', 'pP', 'evt', [evla evlo], ...
        'sta', [stla(ii) stlo(ii)]);
    if isempty(tt)
        theta(ii) = NaN;
        continue
    end
    if strcmpi(model, 'ak135')
        vp = ak135('depths', evdp, 'dcbelow', false).vp;
    elseif strcmpi(model, 'iasp91')
        vp = iasp91('depths', evdp, 'dcbelow', false).vp;
    elseif strcmpi(model, 'prem')
        vp = prem('depths', evdp, 'dcbelow', false).vp;
    else
        warning(['This model (%s) is not implemented for this function' ...
            'yet. PREM is used instead.\n'], model)
        vp = prem('depths', evdp, 'dcbelow', false).vp;
    end
    
    % prevent unphysical values
    sintheta = max(min(tt(1).rayparameter * vp / (R_EARTH - evdp), 1), -1);
    theta(ii) = asin(sintheta) * 180 / pi;
end
end