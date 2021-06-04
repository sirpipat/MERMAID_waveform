function plot_seismogram(filename, ax)
% PLOT_SEISMOGRAM(filename,ax)
% plot a seismogram
%
% INPUT
% filename      Filename of the seismogram
% ax            Axes to plot

% read seismograms
sizeData = [2 Inf];
% convert index to string
fid = fopen(filename,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);
% sampling rate
fs = 1 / (data(1,2) - data(1,1));

% plot the signal
signalplot(data(2,:), fs, data(1,1), ax, removepath(filename), 'left',[]);
xlabel('Time [s]');
end