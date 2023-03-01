function plotpolaritymap(evla, evlo, evdp, model, M, savename, options)
% PLOTPOLARITYMAP(evla, evlo, evdp, model, M, savename, options)
%
% Plots the polarity of P, SV, and SH waves on the map of a given moment
% tensor.
%
% see Dahlen and Tromp, Theoretical Global Seismology, 1998 page 529
%
% INPUT:
% evla      event latitude
% evlo      event longitude
% evdp      event depth in km
% model     Earth model [Default: 'ak135']
% M         moment tensor either 6 elements or 3x3 matrix
%           M=[Mrr Mtt Mpp Mrt Mrp Mtp]
%              / Mrr Mrt Mrp \
%           M=|  Mrt Mtt Mtp  |
%              \ Mrp Mtp Mpp /
% savename  name of the output file of a saved figure
% options   optional parameters stored as a struct with following fields
%   - resolution : [Default: 'crude']
%       'crude'  - 5 degrees
%       'medium' - 2 degrees
%       'fine'   - 1 degree, might be slow
%   - shading : [Default: 'beachball']
%       'colormap'  - use a variety of colors from a colormap
%       'beachball' - use a single color for shading (default)
%   - shadecolor :
%       name of a colormap if shading=='colormap'
%       name of a color or a RGB triplet if shading=='beachball'
%   - numcolors: [Default: 3]
%       integer greater or equal to 2
%   - beachballcolor: [Default: [0 0 0]]
%       name of a color or a RGB triplet of the beachball
%
% SEE ALSO:
% FOCALMECH, CMTPOLARITY
%
% Last modified by sirawich-at-princeton.edu: 03/01/2023

if and(nargin == 1, strcmp(evla, 'demo'))
    evla = -7.4260;
    evlo = 119.8340;
    evdp = 529;
    model = 'ak135';
    M = [1.6800   -0.9090   -0.7660    4.1200    5.6500   -1.7800];
    savename = 'demo';
    plotpolaritymap(evla, evlo, evdp, model, M, savename);
    return
end

defval('model', 'ak135')

% check the optional parameters
if ~exist('options', 'var')
    options.resolution = 'crude';
    options.shading = 'beachball';
    options.shadecolor = [0.4824    0.4078    0.9333];
    options.numcolors = 3;
    options.beachballcolor = [0 0 0];
else
    if ~isfield(options, 'resolution')
        options.resolution = 'crude';
    end
    if ~isfield(options, 'shading')
        options.shading = 'beachball';
    end
    if ~isfield(options, 'shadecolor')
        if strcmpi(options.shading, 'colormap')
            options.shadecolor = 'parula';
        elseif strcmpi(options.shading, 'beachball')
            options.shadecolor = [0.4824    0.4078    0.9333];
        else
            error('Invalid value for options.shading. Use one of these values: ''beachball'' | ''colormap''')
        end
    end
    if ~isfield(options, 'numcolors')
        options.numcolors = 3;
    end
    if ~isfield(options, 'beachballcolor')
        options.beachballcolors = [0 0 0];
    end
end

% convert to 3x3 symmetric matrix
if length(M) == 6
    MM = [M(1) M(4) M(5); M(4) M(2) M(6); M(5) M(6) M(3)];
else
    MM = M;
end

% convert longitude to [0 360]
evlo = mod(evlo, 360);

sname = sprintf('%s_%s.mat', mfilename, ...
    hash([evla evlo evdp double(char(model)) ...
    reshape(MM, [1, numel(MM)]) double(char(savename)) ...
    double(char(options.resolution)) ...
    double(char(options.shading)) ...
    double(char(options.shadecolor)) options.numcolors ...
    double(char(options.beachballcolor))], 'SHA-1'));
pname = fullfile(getenv('IFILES'), 'HASHES', sname);

if ~exist(pname, 'file')
    if strcmpi(options.resolution, 'fine')
        resolution = 1;
    elseif strcmpi(options.resolution, 'medium')
        resolution = 2;
    else
        resolution = 5;
    end
    % station lon/lat grid
    stlo = (0:resolution:360)';
    stla = (-90:resolution:90)';
    [stlo, stla] = meshgrid(stlo, stla);

    % compute the azimuth and epicentral distance
    azim = nan(size(stlo));
    dist = nan(size(stlo));
    theta = nan(size(stlo));
    is_down = false(size(stlo));
    phase = cell(size(stlo));
    for ii = 1:size(stlo, 1)
        for jj = 1:size(stlo, 2)
            [azim(ii,jj), ~, dist(ii,jj)] = azimdist([evlo evla], ...
                [stlo(ii,jj) stla(ii,jj)]);

            tt = indeks(taupTime(model, evdp, ...
                'p,P,Pdiff,PKIKP', 'evt', [evla evlo], ...
                'sta', [stla(ii,jj) stlo(ii,jj)]), 1);
            theta(ii,jj) = tt.incidentDeg;
            is_down(ii,jj) = strcmp(tt.phaseName(1), 'P');
            phase{ii,jj} = tt.phaseName;
        end
    end

    % determine the take-off vector and the polarization vectors
    fp = nan(size(stlo));
    fsv = nan(size(stlo));
    fsh = nan(size(stlo));
    for ii = 1:size(stlo, 1)
        for jj = 1:size(stlo, 2)
            % determine "take-off" unit vector
            azim_rad = azim(ii,jj) * pi / 180;
            theta_rad = theta(ii,jj) * pi / 180;
            if is_down
                v = [-cos(theta_rad); -sin(theta_rad) * cos(azim_rad); ...
                    sin(theta_rad) * sin(azim_rad)];
            else
                v = [ cos(theta_rad); -sin(theta_rad) * cos(azim_rad); ...
                    sin(theta_rad) * sin(azim_rad)];
            end
            % polarization unit vectors
            np = v;
            nsh = cross(v, [1 0 0]') / norm(cross(v, [1 0 0]'), 2);
            nsv = cross(nsh, np);

            % compute radiation pattern
            % Dahlen and Tromp, Theoretical Global Seismology, 1998 page 529 Eq. 12.262
            M0 = sqrt(sum(sum(MM .* MM)) / 2);
            fp(ii,jj) = v' * (MM / M0) * np;
            fsv(ii,jj) = 1/2 * (v' * (MM / M0) * nsv + nsv' * (MM / M0) * v);
            fsh(ii,jj) = 1/2 * (v' * (MM / M0) * nsh + nsh' * (MM / M0) * v);
        end
    end
    % save
    fprintf('save the output to a file to %s ...\n', pname);
    save(pname, 'MM', 'azim', 'dist', 'evdp', 'evla', 'evlo', 'fp',  ...
        'fsh', 'fsv', 'is_down', 'model', 'options', 'phase', 'stla', ...
        'stlo', 'theta', 'savename');
else
    % load
    fprintf('found the save in a file in %s\n', pname);
    fprintf('load the variables ...\n');
    load(pname, 'MM', 'azim', 'dist', 'evdp', 'evla', 'evlo', 'fp',  ...
        'fsh', 'fsv', 'is_down', 'model', 'options', 'phase', 'stla', ...
        'stlo', 'theta', 'savename');
end

% plot
figure
set(gcf, 'Units', 'inches', 'Position', [0 1 6 9])
clf

% color map setting
num_colors = options.numcolors;
if strcmpi(options.shading, 'colormap')
    cmap_all = colormap(options.shadecolor);
    icmap = round(((1:num_colors) - 0.5) / num_colors * size(cmap_all, 1));
    cmap = cmap_all(icmap, :);
else
    c0 = [1 1 1];
    cn = options.shadecolor;
    icmap = (0:num_colors-1)' / (num_colors-1);
    cmap = icmap * cn + (1 - icmap) * c0;
end

% title
ax0 = subplot('Position', [0.1300 0.9500 0.7750 0.0220]);
title(savename) 
set(ax0, 'FontSize', 11, 'Color', 'none');
set(ax0.Title, 'FontSize', 11)
ax0.XAxis.Visible = 'off';
ax0.YAxis.Visible = 'off';
% plot fp
ax1 = subplot('Position', [0.10 0.68 0.84 0.27]);
imagesc(ax1,[0 360],[-90 90],fp);
axis xy
grid on
xticks(0:30:360)
yticks(-90:30:90)
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')
set(ax1, 'CLim', [-1 1], 'TickDir', 'out', 'Box', 'on', 'FontSize', 11)
hold on
addazimdistlines(ax1, evlo, evla);
[axlim,handl,XYZ] = plotcont([0 90], [360 -90], 1, 0);
[handlp, XYp] = plotplates([0 90], [360 -90], 1);
handlp.Color = 'r';
ax1c = colorbar(ax1);
colormap(ax1, cmap)
set(ax1c, 'TickDir', 'out', 'Box', 'on')
ax1b = doubleaxes(ax1);
axes(ax1b)
focalmech(ax1b, MM, evlo, evla, 12, options.beachballcolor);
ax1b.XAxis.Visible = 'off';
ax1b.YAxis.Visible = 'off';
set(ax1b, 'XLim', ax1.XLim, 'YLim', ax1.YLim, 'Color', 'None', ...
    'Position', ax1.Position)
title(ax1, 'P-radiation pattern')

% plot fsv
ax2 = subplot('Position', [0.10 0.36 0.84 0.27]);
imagesc(ax2,[0 360],[-90 90],fsv);
axis xy
grid on
xticks(0:30:360)
yticks(-90:30:90)
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')
set(ax2, 'CLim', [-1 1], 'TickDir', 'out', 'Box', 'on', 'FontSize', 11)
hold on
addazimdistlines(ax2, evlo, evla);
[axlim,handl,XYZ] = plotcont([0 90], [360 -90], 1, 0);
[handlp, XYp] = plotplates([0 90], [360 -90], 1);
handlp.Color = 'r';
ax2c = colorbar(ax2);
colormap(ax2, cmap)
set(ax2c, 'TickDir', 'out', 'Box', 'on')
ax2b = doubleaxes(ax2);
axes(ax2b)
focalmech(ax2b, MM, evlo, evla, 12, options.beachballcolor);
ax2b.XAxis.Visible = 'off';
ax2b.YAxis.Visible = 'off';
set(ax2b, 'XLim', ax2.XLim, 'YLim', ax2.YLim, 'Color', 'None', ...
    'Position', ax2.Position)
title(ax2, 'SV-radiation pattern')

% plot fsh
ax3 = subplot('Position', [0.10 0.04 0.84 0.27]);
imagesc(ax3,[0 360],[-90 90],fsh);
axis xy
grid on
xticks(0:30:360)
yticks(-90:30:90)
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')
set(ax3, 'CLim', [-1 1], 'TickDir', 'out', 'Box', 'on', 'FontSize', 11)
hold on
addazimdistlines(ax3, evlo, evla);
[axlim,handl,XYZ] = plotcont([0 90], [360 -90], 1, 0);
[handlp, XYp] = plotplates([0 90], [360 -90], 1);
handlp.Color = 'r';
ax3c = colorbar(ax3);
colormap(ax3, cmap)
set(ax3c, 'TickDir', 'out', 'Box', 'on')
ax3b = doubleaxes(ax3);
axes(ax3b)
focalmech(ax3b, MM, evlo, evla, 12, options.beachballcolor);
ax3b.XAxis.Visible = 'off';
ax3b.YAxis.Visible = 'off';
set(ax3b, 'XLim', ax3.XLim, 'YLim', ax3.YLim, 'Color', 'None', ...
    'Position', ax3.Position)
title(ax3, 'SH-radiation pattern')

% save figure
set(gcf, 'Renderer', 'painters')
savefile = sprintf('%s_%s.eps', mfilename, savename);
figdisp(savefile, [], [], 2, [], 'epstopdf')
end