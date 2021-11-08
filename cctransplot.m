function cctransplot(ddir1, ddir2, example)
% CCTRANSPLOT(ddir1, ddir2, example)
%
% Plot the transfer function and cross correlation coefficient between the
% z-component displacement at the ocean bottom seismometer and the pressure
% at the hydrophone right above the ocean bottom seismometer. These 
% directories have to be created by SPECFEM2D_INPUT_SETUP_FLAT.
%
% INPUT:
% ddir1         directory to the simulation with the hydrophone output
% ddir2         directory to the simulation with the 3-comp seismograms output
% example       name of the simulation
%
% SEE ALSO:
% SPECFEM2D_INPUT_SETUP_FLAT, RUNFLATSIM
% 
% Last modified by sirawich-at-princeton.edu, 11/05/2021

% read the first arrival of at the ocean bottom
[tims_o, seisdata_o] = getarrivaltemplate(ddir2, example);

% read the hydrophone pessure
[tims_h, seisdata_h] = read_seismogram(sprintf('%sOUTPUT_FILES/%s.%s.PRE.semp', ...
    ddir1, 'AA', 'S0001'));

% tapering window
%w = hann(length(tims_o))';
w = ones(size(tims_o));
% w(1:10) = 0;
% w(end-9:end) = 0;

% compute correlation coefficients
[cc, lags] = xcorr(seisdata_h .* w, seisdata_o .* w, 'coeff');
t_cc = lags * (tims_o(2)-tims_o(1));

% compute the transfer function
tf = ifft((fft(seisdata_h .* w) .* conj(fft(seisdata_o .* w))) ./ ...
    (fft(seisdata_o .* w) .* conj(fft(seisdata_o .* w)) + 1e-42));
t_tf = (0:(length(tf)-1)) * (tims_o(2)-tims_o(1));

% apply transfer function to obtain hydrophone pressure
x_h_tran = conv(seisdata_o, tf);
x_h_tran = x_h_tran(1:length(seisdata_h));

figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [0 6 6 8])

ax1 = subplot('Position', [0.08 0.83 0.86 0.13], 'Box', 'on');
plot(tims_o, seisdata_o, 'Color', 'k')
xlim([tims_o(1) tims_o(end)])
grid on
xlabel('time (s)')
title('First P-wave arrival, z-displacement, OBS')

ax2 = subplot('Position', [0.08 0.57 0.86 0.13], 'Box', 'on');
plot(tims_h, seisdata_h, 'Color', 'k')
xlim([tims_o(1) tims_o(end)])
grid on
hold on
plot(tims_h, x_h_tran, 'Color', 'r')
xlabel('time (s)')
title('Pressure record, MERMAID hydrophone')
legend('observed', 'from response', 'Location', 'northwest')

ax3 = subplot('Position', [0.08 0.32 0.86 0.13], 'Box', 'on');
plot(t_cc, cc, 'Color', 'k')
xlim([tims_o(1) tims_o(end)])
grid on
xlabel('time (s)')
title('Correlation coefficient')

ax4 = subplot('Position', [0.08 0.07 0.86 0.13], 'Box', 'on');
plot(t_tf, tf, 'Color', 'k')
xlim([0 tims_o(end)-tims_o(1)])
grid on
xlabel('time (s)')
title('response')

figure(2)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 6 4])
ax5 = subplot(1,1,1, 'Box', 'on');
plot(tims_h, x_h_tran-seisdata_h, 'Color', 'k')
grid on
xlim([tims_h(1), tims_h(end)])
xlabel('time (s)')
title('Predicted - observed')
end