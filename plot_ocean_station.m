function plot_ocean_station(ddir, network, station)
% PLOT_OCEAN_STATION(ddir, network, station)
% plot seismograms of a station
% 
% INPUT
% ddir      Directory of the seismogram files
% network   Network code
% station   Station number

filename = strcat(ddir, network, sprintf('.S%s.PRE.semp', ...
                  num2str(station,'%4.4i')));
% read seismograms
sizeData = [2 Inf];
fid = fopen(filename,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);
% sampling rate
fs = 1 / (data(1,2) - data(1,1));
% raw signal
x = data(2,:);

% plot the signal
figure;
ax = gca;
signalplot(x, fs, data(1,1), ax, removepath(filename), 'left',[]);
end