function plot_many_stations(ddir,network,station_number,direction,savename)

defval('savename', sprintf('seismograms_%s', direction))

if ~or(direction == 'X', direction == 'Z')
    fprintf('Invalid direction. Exit.\n');
end

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

% index for station
i_n = 1:n;
% keep track for receivers in solid or fluid layer
i_solid = (depth >= water_depth);
i_fluid = (depth < water_depth);

depth_solid = depth(i_solid) - water_depth;
depth_fluid = depth(i_fluid);

figure;
ax = gca;
axis ij
hold on
for ii = 1:length(depth_solid)
    filename = strcat(ddir, network{i_n(ii)}, '.', name{i_n(ii)}, '.BX', direction, '.semd');
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
    % seismograms (in order to keep the relative amplitudes)
    if ii == 1
        scale = 1 * max(abs(xf));
    end
    
    xf = xf / scale * dz;
    
    % plot the signal
    plot(data(1,:), xf + depth_solid(ii), 'Color', 'k');
    hline(gca, depth_solid(ii), '-', 0.5, [0.6 0.6 0.6]);
    % text(150, (ii-0.7), num2str(station_number(ii),'%2.2i'));
end
% add expected arrivals of each phase
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ w 
%  |\  |\  |\  |
%  | \ | \ | \ |          water
%  |  \|  \|  \|
% -------------------------------- b
%  |\   \   \   \         rock
%  o o   o   o   o <- P, PbP, PAwAP, PAwAbAwAP, PAwAbAwAbAwAP respectively
%  |  \   \   \   \
%  +
% direct P-wave (P)
plot((source_depth - depth_solid) / vp, depth_solid, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
% PbP
tPbP = (source_depth + depth_solid) / vp;
plot(tPbP, depth_solid, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
% PAwAP
tPAwAP = tPbP + 2 * water_depth / vw;
plot(tPAwAP, depth_solid, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
% PAwAbAwAP
tPAwAbAwAP = tPAwAP + 2 * water_depth / vw;
plot(tPAwAbAwAP, depth_solid, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)
% PAwAbAwAbAwAP
tPAwAbAwAbAwAP = tPAwAbAwAP + 2 * water_depth / vw;
plot(tPAwAbAwAbAwAP, depth_solid, 'color', [0.7 0.7 0.7], 'LineWidth', 0.5)

xlabel('Time [s]');
ylabel('Depth [m]');
xlim([min(data(1,:)) max(data(1,:))]);
xlim([0 15]);
title(sprintf('Displacement in %s-direction  (scale = %0.3e)', ...
    direction, scale));
ax.XGrid = 'on';
ax.GridAlpha = 0.3;

% print the figure
% save the figure
eps_filename = strcat(savename, '.epsc');
pdf_filename = strcat(savename, '.pdf');
print(eps_filename, '-depsc');
system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
system(sprintf('rm %s', eps_filename));
end