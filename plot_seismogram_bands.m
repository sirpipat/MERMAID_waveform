function plot_seismogram_bands(ddir, network, station,tipe,elem)
% PLOT_SEISMOGRAM_BANDS(ddir, network, station)
% plot a seismogram
%
% INPUT
% ddir      Directory of the seismogram files
% network   Network code
% station   Station number
% tipe      Type of measurement
%           1 = displacement
%           2 = velocity (Not implemented)
%           3 = acceleration (Not implemented)
%           4 = pressure
%           5 = curl of displacement (Not implemented)
%           6 = fluid potential (Not implemented)

types = {'semd','n/a','n/a','semp','n/a','n/a'};

if tipe == 1
    if elem == 'X'
        elem_tag = 'BXX';
    elseif elem == 'Z'
        elem_tag = 'BXZ';
    else
        fprintf('Not valid elem. Exit.');
        return
    end
elseif tipe == 4
    elem_tag = 'PRE';
else
    fprintf('Not valid type. Exit.');
    return
end

filename = strcat(ddir, network, sprintf('.S%s.%s.%s', ...
                  num2str(station,'%4.4i'),elem_tag,types{tipe}));

% read seismograms
sizeData = [2 Inf];
% convert index to string
fid = fopen(filename,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);
% sampling rate
fs = 1 / (data(1,2) - data(1,1));
% raw signal
x = data(2,:);
x_d2 = detrend(decimate(x,2),1);
x_d4 = detrend(decimate(x,4),1);

% filter to various bands
xf1 = bandpass(x,fs,6,10,2,2,'butter','linear');
xf2 = bandpass(x,fs,4,6,2,2,'butter','linear');
xf3 = bandpass(x,fs,2,4,2,2,'butter','linear');
xf4 = bandpass(x,fs,1,2,2,2,'butter','linear');
xf5 = bandpass(x_d2,fs/2,0.5,1,2,2,'butter','linear');
xf6 = bandpass(x_d4,fs/4,0.1,0.5,2,2,'butter','linear');

figure;
% plot the signal
ax1 = subplot(7,1,1);
signalplot(x, fs, data(1,1), ax1, removepath(filename), 'left',[]);
% remove xlabel
nolabels(ax1, 1);
ax1.XAxis.Label.Visible = 'off';
ylims = ax1.YLim;

ax2 = subplot(7,1,2);
signalplot(xf1, fs, data(1,1), ax2, '6-10 Hz', 'left',[]);
% remove xlabel
nolabels(ax2, 1);
ax2.XAxis.Label.Visible = 'off';
ax2.YLim = ylims;

ax3 = subplot(7,1,3);
signalplot(xf2, fs, data(1,1), ax3, '4-6 Hz', 'left',[]);
% remove xlabel
nolabels(ax3, 1);
ax3.XAxis.Label.Visible = 'off';
ax3.YLim = ylims;

ax4 = subplot(7,1,4);
signalplot(xf3, fs, data(1,1), ax4, '2-4 Hz', 'left',[]);
% remove xlabel
nolabels(ax4, 1);
ax4.XAxis.Label.Visible = 'off';
ax4.YLim = ylims;

ax5 = subplot(7,1,5);
signalplot(xf4, fs, data(1,1), ax5, '1-2 Hz', 'left',[]);
% remove xlabel
nolabels(ax5, 1);
ax5.XAxis.Label.Visible = 'off';
ax5.YLim = ylims;

ax6 = subplot(7,1,6);
signalplot(xf5, fs/2, data(1,1), ax6, '0.5-1 Hz', 'left',[]);
% remove xlabel
nolabels(ax6, 1);
ax6.XAxis.Label.Visible = 'off';
ax6.YLim = ylims;

ax7 = subplot(7,1,7);
signalplot(xf6, fs/4, data(1,1), ax7, '0.1-0.5 Hz', 'left',[]);
ax7.YLim = ylims;
end