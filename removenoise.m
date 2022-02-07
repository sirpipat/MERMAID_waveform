function xs = removenoise(t, x, fs, arrival, n)
% xs = REMOVENOISE(x, fs, t0, arrival, nfreq)
%
% Determines the Fourier coefficients of the noise time-series before 
% ARRIVAL which separate the noise-only section from the other section and 
% subtracts the time-series generated from those coefficients from the
% signal.
%
% INPUT:
% t             time
% x             input signal
% fs            rate of sampling
% arrival       when the signal arrives
% n             number of frequencies to fit
%
% OUTPUT:
% xs            output signal with noise removed
%
% Last modified by sirawich-at-princeton.edu. 02/07/2022

% convert all inputs to row vectors
if size(t, 1) > 1
    t = t';
end
if size(x, 1) > 1
    x = x';
    converted = true;
else
    converted = false;
end

% slice for section to determine the noise
xb = x(t <= arrival);
tb = t(t <= arrival);

% determine target frequencies
max_f = fs/2;
min_f = 1/(tb(end)-tb(1));

f = min_f:min_f:(max_f-min_f);
max_n = length(f);
n = min(n, max_n);

% check whether sinefit throws a warning about singular matrix
lastwarn('', '');
[A, B, C, ~, F, P] = sinefit(tb, xb, [], f);
% now if a warning was raised, warnmsg and warnid will not be empty
[warnmsg, warnid] = lastwarn();
if ~isempty(warnid)
    error(warnmsg, warnid);
end

A = [0 sqrt(2)*A];
B = [C sqrt(2)*B];

% find k most significant frequency
[~, I] = maxk(P, n);

A = A(I);
B = B(I);
F = F(I);

% constrct the noise signal
xn = A * sin(2 * pi * F' * t) + B * cos(2 * pi * F' * t);

% remove the noise
xs = x - xn;

% do not forget to convert back
if converted
    xs = xs';
end
end