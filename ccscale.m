function [t_shift, CCmax, lag, CC, s] = ccscale(x1, x2, dt_begin1, dt_begin2, fs, maxmargin)
% [t_shift, CCmax, lag, CC, s] = ...
%   CCSCALE(x1,x2,dt_begin1,dt_begin2,fs,maxmargin)
%
% Computes correlation coefficients for all lags in [-maxmargin, maxmargin]
% between two signals. It also finds the scaling to minimize the misfit of
% the two signals after the best timeshift is applied.
%
% The function expects the relation between two signals to be
% x1(t) == s * x2(t + t_shift)
%
% INPUT:
% x1            A signal containing x2
% x2            A signal contained in x1
% dt_begin1     Begin datetime of x1
% dt_begin2     Begin datetime of x2
% fs            Sampling rate of both signals
% maxmargin     Maximum time shift as a duration [default: seconds(inf)]
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
% Last modified by Sirawich Pipatprathanporn: 11/28/2021

defval('maxmargin', seconds(inf))

% convert x1 and x2 to column vectors
if size(x1, 1) == 1
    x1 = x1';
end
if size(x2, 1) == 1
    x2 = x2';
end

%% find best CC and timeshift
if size(x1, 1) == size(x2, 1)
    [CC, lag] = xcorr(x1, x2, 'coeff');
    lag = lag / fs + seconds(dt_begin2 - dt_begin1);
    
    CC(abs(lag) > seconds(maxmargin)) = 0;

    [CCmax, IImax] = max(CC);
    t_shift = lag(IImax);
else
    [t_shift, CCmax, lag, CC] = ccshift(x1, x2, dt_begin1, dt_begin2, ...
        fs, maxmargin);
end

%% optimal scaling (least-squared method)
% define the time for each sample
dt1 = dt_begin1 + seconds((0:length(x1)-1)' / fs);
dt2 = dt_begin1 + seconds((0:length(x2)-1)' / fs + t_shift);

% determine the start and end of the window
dt_min = max(dt1(1), dt2(1));
dt_max = min(dt1(end), dt2(end));

% get the window to determine the scaling
xs1 = x1(and(dt1 > dt_min, dt1 < dt_max));
xs2 = x2(and(dt2 > dt_min, dt2 < dt_max));

s = (xs2' * xs1) / (xs2' * xs2);
end