function plot_all_ocean_stations(ddir,network,num_station,savename)

defval('savename', 'ocean_seismograms')

% read stations
station_fname = strcat(ddir, 'DATA/STATIONS');
[n, name, network, ~, z] = read_stations(station_fname);
% depth spacing
dz = abs(z(2)-z(1));

% setting parameters
% TODO: figure out how to read the parameters from SOURCE, Par_file etc.
interface_fname = strcat(ddir,'DATA/interfaces_fluid_flat.dat');
source_depth = 4080;
[water_depth, solid_depth] = get_fluid_solid_setting(interface_fname);
vp = 3400;
vw = 1500;

depth = water_depth + solid_depth - z;

figure;
ax = gca;
axis ij
hold on
for jj = 1:n
    filename = strcat(ddir, network{jj}, '.', name{jj}, '.PRE.semp');
    labels = removepath(filename);
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
    % normalize the amplitude by twice of the maximum/minimum of the first
    % seismograms (in order to keep the relative amplitudes) and then
    % scale to the depth spacing
    if jj == 1
        scale = 2 * max(abs(xf));
    end
    xf = xf / scale * dz;
    
    % plot the signal
    plot(td, xf + depth(jj), 'Color', 'k');
    hline(gca, depth(jj), '-', 0.5, [0.6 0.6 0.6]);
    
    % plot the peaks
    if false
        [peaks,locs] = findpeaks(xf,td','MinPeakHeight',scale*0.25);
        [mpeaks,mlocs] = findpeaks(-xf,td','MinPeakHeight',scale*0.25);
%         [peaks,locs] = findpeaks(xf,td','MinPeakWidth',0.25);
%         [mpeaks,mlocs] = findpeaks(-xf,td','MinPeakWidth',0.25);
        scatter(locs, peaks + depth(jj) * scale, 10, 'o', ...
            'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        scatter(mlocs, -mpeaks + depth(jj) * scale, 10, 'o', ...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    end
    
    % plot(data(1,:), data(2,:) + (ii-1) * scale, 'Color', 'k');
    % text(150, (jj-0.7)*scale, labels);
end
% add expected arrivals of each phase

t1 = source_depth / vp + (water_depth - depth) / vw;
t2 = source_depth / vp + (water_depth + depth) / vw;
t3 = t1 + 2 * water_depth / vw;
t4 = t2 + 2 * water_depth / vw;
t5 = t3 + 2 * water_depth / vw;
t6 = t4 + 2 * water_depth / vw;
t7 = t5 + 2 * water_depth / vw;
t8 = t6 + 2 * water_depth / vw;
plot(t1, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t2, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t3, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t4, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t5, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t6, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t7, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
plot(t8, depth, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)

xlabel('Time [s]');
ylabel('Depth [m]');
xlim([min(data(1,:)) max(data(1,:))]);
xlim([0 15]);
title(sprintf('pressure records at ocean stations (scale = %0.3e)', scale));
ax.XGrid = 'on';
ax.GridAlpha = 0.3;
ax.TickDir = 'both';
% adjust YAxisTickLabel
% ax.YTickLabel = string(num2cell(4800 - 150 * ax.YTick));

% print the figure
% save the figure
eps_filename = strcat(savename, '.epsc');
pdf_filename = strcat(savename, '.pdf');
print(eps_filename, '-depsc');
system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
system(sprintf('rm %s', eps_filename));
end