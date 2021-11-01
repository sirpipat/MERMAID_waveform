function drawsetting(ddir, name, savedir, savename, sv)
% DRAWSETTING(ddir, name, savedir, savename, sv)
%
% Draws the background image of the simulation with sources and receivers.
%
% INPUT:
% ddir          main directory of the simulation
% name          name of the model               [Default: removepath(ddir(1:end-1))]
% savedir       directory for the saved file    [Default: $EPS]
% savename      name of the saved file
% sv            whether to save or not          [Default: true]
%
% Last modified by sirawich@princeton.edu, 10/25/2021

defval('name', removepath(ddir(1:end-1)))
interfacefile = [ddir 'DATA/interfaces_' name '.dat'];

figure
set(gcf, 'Units', 'inches', 'Position', [2 2 8 5]);
clf
ax = subplot('Position', [0.05 0.05 0.9 0.9]);
% draw background
ax = drawbackground(interfacefile, ax);

% plot sources
sources = loadsource([ddir 'DATA/SOURCE_' name]);
for ii = 1:length(sources)
    scatter(ax, sources{ii}.xs, sources{ii}.zs, 15, 's', ...
        'MarkerEdgeColor',  rgbcolor('k'), 'MarkerFaceColor', ...
        rgbcolor('yellow'));
end

% plot receivers
[~, ~, ~, x, z] = read_stations([ddir 'DATA/STATIONS']);
scatter(ax, x, z, 15, 's', 'MarkerEdgeColor', [0.1 0.6 0.1], ...
    'MarkerFaceColor', [0.1 0.6 0.1]);

ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

if sv
    % print the figure
    % save the figure
    eps_filename = strcat(savedir, mfilename, '_', savename, '.epsc');
    pdf_filename = strcat(savedir, mfilename, '_', savename, '.pdf');
    print(eps_filename, '-depsc');
    system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
    system(sprintf('rm %s', eps_filename));
end
end