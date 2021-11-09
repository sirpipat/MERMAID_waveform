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
% Last modified by sirawich-at-princeton.edu, 11/09/2021

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
SEISDATA_O = fft(seisdata_o .* w);
SEISDATA_H = fft(seisdata_h .* w);
% Use damping factor
d = abs(max(seisdata_o)).^2 * 1000;
% d = 10 .^ (-45:0.1:-35)';
% % minimize General Cross Validation function (GCV)
% X = (abs(SEISDATA_O).^2) ./ (abs(SEISDATA_O).^2 + d);
% GCV = sum((SEISDATA_H .* (1 - X)).^2, 2) ./ ...
%     (length(SEISDATA_O) - sum(X, 2)) .^2;
% d = d(GCV == min(GCV));

tf = ifft(SEISDATA_H .* conj(SEISDATA_O) ./ ...
    (SEISDATA_O .* conj(SEISDATA_O) + d));
t_tf = (0:(length(tf)-1)) * (tims_o(2)-tims_o(1));

% apply transfer function to obtain hydrophone pressure
x_h_tran = conv(seisdata_o, tf);
x_h_tran = x_h_tran(1:length(seisdata_h));

figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [0 6 6 9])

ax1 = subplot('Position', [0.08 0.86 0.86 0.12], 'Box', 'on');
plot(tims_o, seisdata_o, 'Color', 'k')
xlim([tims_o(1) tims_o(end)])
grid on
xlabel('time (s)')
title('First P-wave arrival, z-displacement, OBS')

ax2 = subplot('Position', [0.08 0.66 0.86 0.12], 'Box', 'on');
plot(tims_h, seisdata_h, 'Color', 'k')
xlim([tims_o(1) tims_o(end)])
grid on
hold on
plot(tims_h, x_h_tran, 'Color', 'r')
xlabel('time (s)')
title('Pressure record, MERMAID hydrophone')
legend('observed', 'from response', 'Location', 'northwest')

ax3 = subplot('Position', [0.08 0.46 0.86 0.12], 'Box', 'on');
plot(tims_h, x_h_tran-seisdata_h, 'Color', 'k')
grid on
xlim([tims_h(1), tims_h(end)])
xlabel('time (s)')
title('Predicted - observed')

ax4 = subplot('Position', [0.08 0.26 0.86 0.12], 'Box', 'on');
plot(t_cc, cc, 'Color', 'k')
xlim([tims_o(1) tims_o(end)])
grid on
xlabel('time (s)')
title('Correlation coefficient')

ax5 = subplot('Position', [0.08 0.06 0.86 0.12], 'Box', 'on');
plot(t_tf, tf, 'Color', 'k')
xlim([0 tims_o(end)-tims_o(1)])
grid on
xlabel('time (s)')
title(sprintf('response, dampling factor = %0.5g', d))

% save figure
set(gcf, 'Renderer', 'painters')
savename = sprintf('%s_%s.eps', mfilename, example);
figdisp(savename, [], [], 2, [], 'epstopdf');
end