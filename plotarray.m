function fig = plotarray(ddir, network, savedir, savename, sv)
% fig = PLOTARRAY(ddir, network, savedir, savename, sv)
%
% Plots seismograms from all stations in a network.
%
% INPUT:
% ddir          main directory of the simulation
% network       network name
% savedir       directory for the saved file    [Default: $EPS]
% savename      name of the saved file
% sv            whether to save or not          [Default: true]
%
% OUTPUT:
% fig           figure handle for the plot
%
% Last modified by sirawich@princeton.edu, 08/30/2021

defval('savedir', getenv('EPS'))
defval('savename', network)
defval('sv', true)

% identify all seismogram files
seisfiles = ls2cell([ddir 'OUTPUT_FILES/' network '*'], 1);

words = split(removepath(seisfiles{1}),'.');
tipe = words{3};

switch tipe
    case 'BXX'
        title_string = 'Displacement in X-direction';
    case 'BXZ'
        title_string = 'Displacement in Z-direction';
    case 'PRE'
        title_string = 'Pressure';
    otherwise
        title_string = 'Record';
end

% plot
fig = figure;
set(fig, 'Units', 'inches', 'Position', [2 2 6 8], 'Renderer', 'painters');
clf
ax1 = subplot('Position', [0.10 0.48 0.88 0.44], 'Box', 'on', ...
    'TickDir', 'both', 'FontSize', 12);
hold on
% find absolute maximum of all time-series
xmax_all = 0;
for ii = 1:length(seisfiles)
    [~, x] = read_seismogram(seisfiles{ii});
    xmax = 1 * max(abs(x), [], 'all');
    if xmax > xmax_all
        xmax_all = xmax;
    end
end

% plot seismograms
for ii = 1:length(seisfiles)
    [t, x] = read_seismogram(seisfiles{ii});
    x = x / xmax_all;
    plot(ax1, t, (0*x) + ii, 'Color', [0.75 0.75 0.75]);
    plot(ax1, t, x + ii, 'k');
    xlim([0 t(end)]);
    xlabel('time (s)');
    ylabel('station number');
    title(sprintf('Network %s : %s', network, title_string));
end
ax1.YLim = [0 length(seisfiles)+1];

% draw the setting
ax2 = subplot('Position', [0.02 0.02 0.96 0.4], 'Box', 'on', 'FontSize', 12);
example = removepath(ddir(1:end-1));
interfacefile = [ddir 'DATA/interfaces_' example '.dat'];
ax2 = drawbackground(interfacefile, ax2);
% locate stations in the array
[~, ~, networks, x, z] = read_stations([ddir 'DATA/STATIONS']);
x_network = x(strcmp(networks, network));
z_network = z(strcmp(networks, network));
scatter(ax2, x_network, z_network, 5, 's', 'MarkerEdgeColor', ...
    [0.1 0.6 0.1], 'MarkerFaceColor', [0.1 0.6 0.1]);
scatter(ax2, x_network(1), z_network(1), 15, 's', 'MarkerEdgeColor', ...
    [0.1 0.1 0.1], 'MarkerFaceColor', [0.8 0.1 0.1], ...
    'LineWidth', 0.25);
scatter(ax2, x_network(end), z_network(end), 15, 's', 'MarkerEdgeColor', ...
    [0.1 0.1 0.1], 'MarkerFaceColor', [0.8 0.8 0.1], ...
    'LineWidth', 0.25);

% label first and last station
if any(abs(x_network(2:end) - x_network(1:end-1)) < 1e-5 * (ax2.XLim(2) - ax2.XLim(1)))
    x_shift =  0.015 * (ax2.XLim(2) - ax2.XLim(1));
else
    x_shift = -0.01 * (ax2.XLim(2) - ax2.XLim(1));
end
if any(abs(z_network(2:end) - z_network(1:end-1)) < 1e-5 * (ax2.YLim(2) - ax2.YLim(1)))
    z_shift = 0.05 * (ax2.YLim(2) - ax2.YLim(1));
else
    z_shift = 0;
end
width = 0.03 * (ax2.XLim(2) - ax2.XLim(1));
height = 0.065 * (ax2.YLim(2) - ax2.YLim(1));
rectangle('Position', [x_network(1) + x_shift - width/5, ...
    z_network(1) + z_shift - height/2, ...
    width, height], 'EdgeColor', 'k', 'FaceColor', 'w');
text(ax2, x_network(1) + x_shift, z_network(1) + z_shift, '1');
width = 0.036 * (ax2.XLim(2) - ax2.XLim(1));
rectangle('Position', [x_network(end) + x_shift - width/8, ...
    z_network(end) + z_shift - height/2, ...
    width, height], 'EdgeColor', 'k', ...
    'FaceColor', 'w');
text(ax2, x_network(end) + x_shift, z_network(end) + z_shift, ...
    string(length(seisfiles)));

% plot decoration
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';

if sv
    figure(gcf)
    % print the figure
    % save the figure
    eps_filename = strcat(savedir, mfilename, '_', savename, '.epsc');
    pdf_filename = strcat(savedir, mfilename, '_', savename, '.pdf');
    print(eps_filename, '-depsc');
    system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
    system(sprintf('rm %s', eps_filename));
end
end
