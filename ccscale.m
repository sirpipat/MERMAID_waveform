function [t_shift, CCmax, lag, CC, s] = ...
    ccscale(x1, x2, dt_begin1, dt_begin2, fs, maxmargin, windowtype, cc_env)
% [t_shift, CCmax, lag, CC, s] = ...
%   CCSCALE(x1,x2,dt_begin1,dt_begin2,fs,maxmargin, windowtype, cc_env)
%
% Computes correlation coefficients for all lags in [-maxmargin, maxmargin]
% between two signals. It also finds the scaling to minimize the misfit of
% the two signals after the best timeshift is applied.
%
% The function expects the relation between two signals to be
% x1(t) == s * x2(t - t_shift)
%
% INPUT:
% x1            A signal containing x2
% x2            A signal contained in x1
% dt_begin1     Begin datetime of x1
% dt_begin2     Begin datetime of x2
% fs            Sampling rate of both signals
% maxmargin     Maximum time shift as a duration [default: seconds(inf)]
% windowtype    How to apply window to CC outside [-maxmargin maxmargin]
%               options are the following
%               'hard' -- boxcar window, zero outside  [default]
%               'soft' -- Gaussian curve outside the window with the
%                         variance of (maxmargin / 2)^2
% cc_env        Whether to compare envelope instead of waveform 
%               [default: true]
%
% OUTPUT:
% t_shift       Best time shift where CC is maximum
% CCmax         Maximum correlation coefficient
% lag           Vector of all time shifts
% CC            Vector of CC for every time shift in lag
% s             Scaling to minimize the misfit
%
% SEE ALSO:
% CCSHIFT, XCORR
%
% Last modified by Sirawich Pipatprathanporn: 12/05/2023

defval('maxmargin', seconds(inf))
defval('windowtype', 'hard')
defval('cc_env', true)

% convert x1 and x2 to column vectors
if size(x1, 1) == 1
    x1 = x1';
end
if size(x2, 1) == 1
    x2 = x2';
end

if cc_env
    x1 = envelope(x1);
    x2 = envelope(x2);
end

%% find best CC and timeshift
if size(x1, 1) == size(x2, 1)
    [CC, lag] = xcorr(x1, x2, 'coeff');
    lag = lag / fs + seconds(dt_begin1 - dt_begin2);
    
    switch lower(windowtype)
        case 'hard'
            CC(abs(lag) > seconds(maxmargin)) = 0;
        case 'soft'
            % soft window (flat within maxmargin, normal outside) to put less
            % weight on peaks outside the maxmargin window
            mm = seconds(maxmargin);
            w = exp(-(max(abs(lag), mm) - mm).^2 / (2 * (mm/2).^2))';
            CC = CC .* w;
        otherwise
            fprintf('Invalid option. Hard window is applied\n');
            CC(abs(lag) > seconds(maxmargin)) = 0;
    end
    [CCmax, IImax] = max(CC);
    t_shift = lag(IImax);
else
    [t_shift, CCmax, lag, CC] = ccshift(x1, x2, dt_begin1, dt_begin2, ...
        fs, maxmargin, windowtype);
end

%% optimal scaling (rms ratio)
% define the time for each sample
dt1 = dt_begin1 + seconds((0:length(x1)-1)' / fs);
dt2 = dt_begin2 + seconds((0:length(x2)-1)' / fs);

% determine the start and end of the window
dt_min = max(dt1(1), dt2(1) + seconds(t_shift));
dt_max = min(dt1(end), dt2(end) + seconds(t_shift));

% get the window to determine the scaling
ep = seconds(0.01/fs);
xs1 = x1(and(geq(dt1, dt_min, ep), leq(dt1, dt_max, ep)));
xs2 = x2(and(geq(dt2  + seconds(t_shift), dt_min, ep), ...
    leq(dt2 + seconds(t_shift), dt_max, ep)));

if cc_env
    s = rms(xs1) / rms(xs2);
else
    s = std(xs1) / std(xs2);
end
end

function r = leq(a, b, ep)
r = or(a < b, abs(a - b) < ep);
end

function r = geq(a, b, ep)
r = or(a > b, abs(a - b) < ep);
end