function [x, z] = bathymetryprofile(width, npts, lonlat, az)
% [x, z] = bathymetryprofile(width, npts, [lon lat], az)
%
% Generate a bathymetry profile centered at [LON LAT] along a straight path
% in a specified azimuthal direction. The azimuth will be the direction
% from the left to the right of the profile.
%
% INPUT
% width         length of the path
% npts          number of points on the path
% [lon lat]     longitude and latitude of the center of the path
% az            azimuthal direction of the path
%
% OUTPUT
% x             x-coordinate along the path from 0 to LENGTH
% z             elevation
%
% EXAMPLE
% % demo
% bathymetryprofile('demo');
%
% % basic example call
% lon   = -171.9965;
% lat   = -12.0744;
% az    = 90;
% npts  = 401;
% width = 20000;
% [x, z] = bathymetryprofile(width, npts, [lon lat], az);
%
% % plot the bathymetry
% plot(x, z)
%
% SEE ALSO
% BATHYMETRY, BATHYMETRYPROFILE2D
% 
% Last modified by sirawich-at-princeton.edu, 07/01/2024

% Earth's radius in meter
R = 6371000;

% convert meters to degrees
m2deg = 180 / pi / R;

%% demo
if ischar(width) && strcmp(width, 'demo')
    lon = -171.9965;%-175.13;%-159.74;%-175.18;%-174.9;
    lat = -12.0744;%-13.75;%-8.88;%-25.925;%-14.425;
    azs = 0:45:135;
    npts = 401;
    width = 20000;
    % half length in degrees
    halflength = width / 2 * m2deg;
    hl = 1.02 * halflength;
    % half widith of the small map in degrees
    hw = 1;
    
    figure(11)
    clf
    set(gcf, 'Units', 'inches', 'Position', [12 8 8 8])
    % bathymetry
    ax1 = subplot(12,7,[1,21]);
    hold on

    % orientation
    ax2 = subplot(12,7,[29,52]);
    [~, ~, ~, ~, c, xoffset] = bathymetry([], [-hl hl] + lon, ...
        [-hl hl] + lat, true, ax2);
    hold on
    
    % plot small map
    ax3 = subplot(12,7,[33,56]);
    [~, ~, ~, ~, c3, xoffset3] = bathymetry([], [-hw hw] + lon, ...
        [-hw hw] + lat, true, ax3);
    hold on
    [xbox, ybox] = boxcorner(mod(1 * [-hl hl] + lon + xoffset3, 360), ...
        1 * [-hl hl] + lat);
    plot(xbox, ybox, 'LineWidth', 2, 'Color', 'y')
    xlabel('longitude (degrees)')
    ylabel('latitude (degrees)')
    
    % plot large map
    ax4 = subplot(12,7,[64,84]);
    [~, ~, ~, ~, c4, xoffset4] = bathymetry([], [120 270], ...
        [-30 10], true, ax4);
    hold on
    [xbox, ybox] = boxcorner(mod([-hw hw] + lon + xoffset4, 360), ...
        [-hw hw] + lat);
    plot(xbox, ybox, 'LineWidth', 2, 'Color', 'y')
    xlabel('longitude (degrees)')
    ylabel('latitude (degrees)')
    
    for ii = 1:length(azs)
        % plot the bathymetry (cross section)
        [x, z] = bathymetryprofile(width, npts, [lon lat], azs(ii));
        plot(ax1, x, z, 'LineWidth', 1);
        % plot the line in the map
        plot(ax2, lon + halflength * cos(pi/2-azs(ii)*pi/180) * ...
            [-1 1] + xoffset, lat + halflength * ...
            sin(pi/2-azs(ii)*pi/180) * [-1 1], 'LineWidth', 1, ...
            'Color', 'k');
    end
    
    % axes decoration
    axes(ax1)
    box on
    grid on
    xlabel('location (m)')
    ylabel('elevation below sea level (m)')
    title('bathymetry')
    set(ax1, 'TickDir', 'both', 'FontSize', 11)
    legend(ax1, '180 -- 0', '225 -- 45', '270 -- 90', '315 -- 135', ...
        'Location', 'best')
    
    axes(ax2)
    xlabel('longitude (degrees)')
    ylabel('latitude (degrees)')
    set(ax2, 'TickDir', 'both', 'FontSize', 11)
    c.Label.String = 'elevation (m)';
    c.Label.FontSize = 12;
    hold on
    colormap(ax2, flip(kelicol, 1));
    ax2.CLim = ax1.YLim;
    
    axes(ax3)
    c3.Label.String = 'elevation (m)';
    c3.Label.FontSize = 12;
    hold on
    set(ax3, 'TickDir', 'both', 'FontSize', 11)
    
    axes(ax4)
    c4.Label.String = 'elevation (m)';
    c4.Label.FontSize = 12;
    hold on
    set(ax4, 'TickDir', 'both', 'FontSize', 11)
    
    % save the figure
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_demo.eps', mfilename), [], [], 2, [], 'epstopdf');
    return
end

%% computation
% x-coordinate of the bathymetry
x = linspace(-width/2, width/2, npts)';

% latitude and longitude of the bathymetry
[lats, lons] = reckon(lonlat(2), lonlat(1), x * m2deg, az);

% This handles cases where lons span across 180 E/W longitude
% Determine if there is any jump in lons
difflons = lons(2:end) - lons(1:(end-1));
if ~or(all(difflons < 0), all(difflons > 0))
    lons = mod(lons, 360);
    is_across180 = true;
else
    is_across180 = false;
end

% request grided elevation
[longrid, latgrid, elev] = bathymetry([], ...
    [min(lons) max(lons)] + [-0.1 0.1], ...
    [min(lats) max(lats)] + [-0.1 0.1], ...
    false);

% convert longitude to [-180 180] when lons do not span across 180 E/W
if ~is_across180
    longrid = mod(longrid + 180, 360) - 180;
end

% handle edge cases where longrid span across 180 E/W longitude
if longrid(1) > longrid(end)
    longrid = mod(longrid, 360);
    lons = mod(lons, 360);
end

% get the elevation along the straight path
z = interp2(latgrid, longrid, double(elev), lats, lons);

% shift x-cooridnate to [0 width]
x = x - x(1);
end