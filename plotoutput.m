function plotoutput(ddir, name)
% PLOTOUTPUT(ddir)
%
% Makes various plots from the output. Then, saves to ddir/PLOTS/.
%
% INPUT:
% ddir          simulation directory
% name          name for the model
%
% Last modified by sirawich@princeton.edu, 10/26/2021

defval('name', removepath(ddir(1:end-1)));

% makes a directory for the plots if there is not one
savedir = [ddir 'PLOTS/'];
system(['mkdir ' savedir]);

% plot setting
drawsetting(ddir, name, savedir, name, true);

% reads receiversets
[~, ~, networks, ~, ~] = read_stations([ddir 'DATA/STATIONS']);
networks = unique(networks);

for ii = 1:length(networks)
    plotarray(ddir, name, networks{ii}, savedir, networks{ii}, true);
end

% make animation
animatepropagation(ddir, savedir);

% plot water sound speed profile if possible
try
    % load the supplementary file
    load(sprintf('%sDATA/supplementary_%s.mat', ddir, name), ...
         'water_model');
    [cz,z] = munk(water_model.zm, water_model.zc, water_model.dz, ...
        water_model.B);
    figure
    set(gcf, 'Units', 'inches', 'Position', [2 2 6 8]);
    clf
    ax = subplot('Position', [0.08 0.08 0.84 0.84], 'FontSize', 12);
    cla
    % plot the profile
    plot(ax, cz, z, 'LineWidth', 1, 'Color', 'k');
    axis ij
    grid on
    xlabel('sound speed (m/s)')
    ylabel('depth (m)')
    title('water sound speed profile')
    
    figure(gcf)
    % print the figure
    % save the figure
    eps_filename = strcat(savedir, mfilename, '_water_profile.epsc');
    pdf_filename = strcat(savedir, mfilename, '_water_profile.pdf');
    print(eps_filename, '-depsc');
    system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
    system(sprintf('rm %s', eps_filename));
catch
    fprintf('Cannot find the supplementary file\n');
end
end
