function bandstopdemo(x,Fs,colo,cohi,npol,npas,tipe,trending)
% BANDSTOPDEMO(x,Fs,colo,cohi,npol,npas,tipe,trending)
%
% Compare the MATLAB built-in verions of BANDSTOP with the FJS
% implementation of BANDSTOP as well as the complementary signal (i.e. the
% original signal subtracted by FJS implementation of BANDPASS).
%
% INPUT:
%
% x         The signal
% Fs        Its sampling frequency [Hz]
% colo      The lower corner frequency [Hz]
% cohi      The higher corner frequency [Hz]
% npol      The number of poles [default: 2]
% npas      The number of passes [default: 1]
% tipe      The filter name [default: 'butter']
% trending  'linear' or 'constant' [default: 'linear']
%
% NOTE: 
%
% MATLAB-own BANDSTOP exists, but I write this function to follow the
% function signature of Frederik J Simons's BANDPASS
%
% Returns the npas frequency response and the effective stop band for
% one or two passes (3 dB level)
%
% You'll see that plot(F,decibel(HABS2)) (this is what FREQZ plots)
% shows how' you concentrate between cohi and colo at the 
% 3 dB-level
%
% SEE ALSO:
% BANDSTOP, BANDPASS
%
% Last modified by sirawich@princeton.edu, 10/24/2022

defval('npol',2)
defval('npas',1)
defval('colo',0.05)
defval('cohi',0.50)
defval('Fs',110)
defval('tipe','butter')
defval('trending','linear')

% bandpass signal
xf = bandpass(x, Fs, colo, cohi, npol, npas, tipe, trending);

% complementary signal (original - bandpass)
xc = x - xf;

% path to all bandstop.m
paths = which('bandstop', '-all');

% MATLAB's BANDSTOP
currpath = pwd;
cd(fileparts(paths{2}));
xstop_MATLAB = bandstop(x, [colo cohi], Fs);

% BANDSTOP version of FJS's BANDPASS
cd(fileparts(paths{1}));
xstop_FJS = bandstop(x, Fs, colo, cohi, npol, npas, tipe, trending);
cd(currpath);

% comparison
fig = figure;
set(fig, 'Units', 'inches', 'Position', [0 1 6 8]);
ax1 = subplot('Position', subplotposition(4, 1, 1, ...
    [0.08 0.08 0.08 0.08], [0.03 0.05 0.03 0.03]));
plot(x, 'LineWidth', 1)
grid on
title('original')
set(ax1, 'TickDir', 'out')
nolabels(ax1, 1);

ax2 = subplot('Position', subplotposition(4, 1, 2, ...
    [0.08 0.08 0.08 0.08], [0.03 0.05 0.03 0.03]));
plot(xf, 'LineWidth', 1)
grid on
title(sprintf('FJS bandpass (bp): %.2f-%.2f Hz, %d poles, %d pass', colo, ...
    cohi, npol, npas))
set(ax2, 'TickDir', 'out')
nolabels(ax2, 1);

ax3 = subplot('Position', subplotposition(4, 1, 3, ...
    [0.08 0.08 0.08 0.08], [0.03 0.05 0.03 0.03]));
plot(xstop_FJS, 'LineWidth', 1.75)
hold on
plot(xstop_MATLAB, 'LineWidth', 1)
plot(xc, 'LineWidth', 0.75)
grid on
title('bandstop (bs)')
set(ax3, 'TickDir', 'out')
nolabels(ax3, 1);
legend('FJS bs', 'MATLAB bs', ...
    'complement: original - FJS bp', 'Location', 'northwest')

ax4 = subplot('Position', subplotposition(4, 1, 4, ...
    [0.08 0.08 0.08 0.08], [0.03 0.05 0.03 0.03]));
plot(xstop_FJS - xstop_MATLAB, 'LineWidth', 1)
hold on
plot(xstop_FJS - xc, 'LineWidth', 1)
grid on
xlabel('sample number')
title('differences')
set(ax4, 'TickDir', 'out')
legend('FJS bs - MATLAB bs', 'FJS bs - complement', ...
    'Location', 'northwest')

% save the figure
set(fig, 'Renderer', 'painters')
figdisp(sprintf('%s.eps', mfilename), [], [], 2, [], 'epstopdf')
end

