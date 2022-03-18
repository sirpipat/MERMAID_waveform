function [t_cc, cc, t_rf, rf] = cctransplot(ddir1, ddir2, example, channel1, channel2, fs, plt)
% [t_cc, cc, t_rf, rf] = CCTRANSPLOT(ddir1, ddir2, example, channel1, channel2, fs, plt)
%
% Plot the transfer function and cross correlation coefficient between the
% z-component displacement at the ocean bottom seismometer and the pressure
% at the hydrophone right above the ocean bottom seismometer. These 
% directories have to be created by SPECFEM2D_INPUT_SETUP_FLAT.
%
% INPUT:
% ddir1         directory to the simulation with the hydrophone output
% ddir2         directory to the simulation with the 3-comp seismograms output
%               For 'devel' branch where the simulation can handle multiple
%               type of receivers, use the same directory for both DDIR1
%               and DDIR2.
% example       name of the simulation
% channel1      input channel       {RECIEVER TYPE} 
%                   [default: {'bottom' 'displacement'}]
% channel2      output channel      {RECEIVER TYPE}
%                   [default: {'hydrophone' 'pressure'}]
%
%       RECEIVER could be the following:
%           'bottom'
%           'hydrophone'
%       TYPE could be the following:
%           'displacement'
%           'pressure'
%
% fs            sampling rate of the output [Default: the sampling rate of
%               the two channels]
% plt           whether to plot or not [Default: true]
%
% OUTPUT:
% t_cc      time for cross correlation coefficient
% cc        cross correlation coefficient
% t_rf      time for response function
% rf        response funtion
%
% SEE ALSO:
% SPECFEM2D_INPUT_SETUP_FLAT, RUNFLATSIM
% 
% Last modified by sirawich-at-princeton.edu, 03/18/2022

defval('channel1', {'bottom' 'displacement'})
defval('channel2', {'hydrophone' 'pressure'})
defval('plt', true)

% read the first arrival of the input signal
if strcmpi(channel1{2}, 'displacement')
    [tims_i, seisdata_i] = getarrivaltemplate(ddir2, example, channel1{1});
else
    [tims_i, seisdata_i] = getarrivaltemplate(ddir1, example, channel1{1});
end

% read the output signal
if strcmpi(channel2{2}, 'displacement')
    ddir = ddir2;
    chan = 'BXZ.semd';
else
    ddir = ddir1;
    chan = 'PRE.semp';
end
if strcmpi(channel2{1}, 'bottom')
    network = 'AB';
else
    network = 'AA';
end
[tims_o, seisdata_o] = read_seismogram(sprintf('%sOUTPUT_FILES/%s.%s.%s', ...
    ddir, network, 'S0001', chan));

% lowpass the signals
dt = tims_i(2) - tims_i(1);
seisdata_i = lowpass(seisdata_i, 1/dt, 5, 2, 2, 'butter', 'linear');
seisdata_o = lowpass(seisdata_o, 1/dt, 5, 2, 2, 'butter', 'linear');

% resample the data to the given sampling rate
if ~isempty(fs)
    dt = 1/fs;
    tt = (tims_i(1):dt:tims_i(end))';
    seisdata_i = shannon(tims_i, seisdata_i, tt);
    seisdata_o = shannon(tims_o, seisdata_o, tt);
    tims_i = tt;
    tims_o = tt;
end

% number of frequencies
% N = length(tims_i);
N = length(tims_i);
nf = 2 ^ nextpow2(N);

% tapering window
w = shanning(length(tims_i), 0.05);
% w = ones(size(tims_i));
% w(1:10) = 0;
% w(end-9:end) = 0;

% compute correlation coefficients
[cc, lags] = xcorr(seisdata_o .* w, seisdata_i .* w, 'coeff');
t_cc = lags' * dt;

asq = max(abs(seisdata_i)).^2;

[rf, x_h_tran, n, d] = spectraldivision(seisdata_o, seisdata_i, w, 'damp', '');
t_rf = (0:(N-1))' * dt;

if plt
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 6 6 9])

    ax1 = subplot('Position', [0.08 0.86 0.86 0.12]);
    plot(tims_i, seisdata_i, 'Color', 'k')
    xlim([tims_i(1) tims_i(end)])
    grid on
    xlabel('time (s)')
    title(sprintf('First P-wave arrival, %s, %s', channel1{2}, channel1{1}))
    set(ax1, 'Box', 'on', 'TickDir', 'both');

    ax2 = subplot('Position', [0.08 0.66 0.86 0.12]);
    plot(tims_o, seisdata_o, 'Color', 'k')
    xlim([tims_i(1) tims_i(end)])
    grid on
    hold on
    plot(tims_o, x_h_tran, 'Color', [0.9 0.3 0.1])
    xlabel('time (s)')
    title(sprintf('%s record, %s', channel2{2}, channel2{1}))
    legend('observed', 'from response', 'Location', 'northeast')
    set(ax2, 'Box', 'on', 'TickDir', 'both');

    ax3 = subplot('Position', [0.08 0.46 0.86 0.12]);
    plot(tims_o, x_h_tran-seisdata_o, 'Color', 'k')
    grid on
    xlim([tims_o(1), tims_o(end)])
    xlabel('time (s)')
    title(sprintf('from response - observed, norm = %0.5g', n))
    set(ax3, 'Box', 'on', 'TickDir', 'both');

    ax4 = subplot('Position', [0.08 0.26 0.86 0.12]);
    plot(t_cc, cc, 'Color', 'k')
    xlim([tims_i(1) tims_i(end)])
    grid on
    xlabel('time (s)')
    title('Correlation coefficient')
    set(ax4, 'Box', 'on', 'TickDir', 'both');

    ax5 = subplot('Position', [0.08 0.06 0.86 0.12]);
    plot(t_rf, rf, 'Color', 'k')
    xlim([0 tims_i(end)-tims_i(1)])
    grid on
    xlabel('time (s)')
    title(sprintf('response, damping factor = %0.5g x input amplitude^2', d/asq))
    set(ax5, 'Box', 'on', 'TickDir', 'both');

    % save figure
    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%s_in_%s-%s_out_%s-%s.eps', mfilename, example, ...
        channel1{1}, channel1{2}, channel2{1}, channel2{2});
    figdisp(savename, [], [], 2, [], 'epstopdf');
end
end