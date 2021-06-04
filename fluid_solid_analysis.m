function fluid_solid_analysis(ddir_fluid, ddir_solid)
% FLUID_SOLID_ANALYSIS(ddir_fluid, ddir_solid)
%
% Performs the following (more to come)
% 1. make videos from simulation snapshots photos
% 2. plot seismograms (pressure for hydrophone in the water column) and (X
% and Z components for seismometers burried under the ocean bottom)
% 3. compare a ocean-bottom hydrophone to a ocean-bottom seismometer
%
% INPUT
% ddir_fluid        directory to a simulation with hydrophones output
% ddir_solid        directory to the same simulation but with burried
%                   seismometers output

if true
    % make videos from simulation snapshots
    ddirs = {ddir_fluid, ddir_solid};
    for ii = 1:2
        ddir = ddirs{ii};
        [allfiles, fndex] = allfile(ddir);
        first_filename = strcat(ddir, 'forward_image0000005.jpg');
        last_filename = strcat(ddir, 'forward_image0050000.jpg');

        % determine the begin and end indices of the snapshots
        first_index = strcmp(allfiles, first_filename) * (1:fndex)';
        last_index = strcmp(allfiles, last_filename) * (1:fndex)';
        images = allfiles(first_index:last_index);

        % frame rate = # of snapshots / (# of steps * stepping size)
        fr = 500/(100000*0.00085);
        savename = strcat('wave_propagation_', removepath(ddir(1:end-1)));
        images2video(images, 'jpg', fr, savename);
    end

    % plot seismograms
    plot_all_ocean_stations(ddir_fluid, 'AA', 32, ...
        strcat('ocean_seismograms_', removepath(ddir_fluid(1:end-1))));
    plot_many_stations(ddir_solid, 'AA', 1:11, 'X', ...
        strcat('seismograms_X_', removepath(ddir_solid(1:end-1))))
    plot_many_stations(ddir_solid, 'AA', 1:11, 'Z', ...
        strcat('seismograms_Z_', removepath(ddir_solid(1:end-1))))
end
% compare a ocean-bottom hydrophone to a ocean-bottom seismometer
figure
ax1 = subplot(2,1,1);
filename = strcat(ddir_fluid, 'AA.S0001.PRE.semp');
% read seismograms
sizeData = [2 Inf];
fid = fopen(filename,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);

if false
    % filter the signal
    fs = 1 / (data(1,2) - data(1,1));
    d_factor = 1;
    xd = detrend(decimate(detrend(data(2,:), 1), d_factor), 1);
    td = decimate(data(1,:), d_factor);
    xf = bandpass(xd, fs/d_factor, 0.05, 1, 2, 2, 'butter', 'linear');
else
    td = data(1,:);
    xf = data(2,:);
end

plot(td, xf, 'Color', 'k');
xlabel('Time [s]');
xlim([min(data(1,:)) max(data(1,:))]);
title('pressure records at ocean station');
grid on
ax1.GridAlpha = 0.3;

ax2 = subplot(2,1,2);
filename = strcat(ddir_solid, 'AA.S0011.BXZ.semd');
% read seismograms
sizeData = [2 Inf];
fid = fopen(filename,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);

if false
    % filter the signal
    fs = 1 / (data(1,2) - data(1,1));
    d_factor = 1;
    xd = detrend(decimate(detrend(data(2,:), 1), d_factor), 1);
    td = decimate(data(1,:), d_factor);
    xf = bandpass(xd, fs/d_factor, 0.05, 1, 2, 2, 'butter', 'linear');
else
    td = data(1,:);
    xf = data(2,:);
end

plot(td, xf, 'Color', 'k');
xlabel('Time [s]');
xlim([min(data(1,:)) max(data(1,:))]);
title('Displacement in Z-direction');
grid on
ax2.GridAlpha = 0.3;
end