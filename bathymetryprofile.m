function [x, z] = bathymetryprofile(width, npts, lonlat, az)
% [x, z] = bathymetryprofile(width, npts, lonlat, az)
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
% BATHYMETRY
% 
% Last modified by sirawich-at-princeton.edu, 03/18/2022

% Earth's radius in meter
R = 6371000;

% convert meters to degrees
m2deg = 180 / pi / R;

%% demo
if ischar(width) && strcmp(width, 'demo')
    lon = -171.9965;
    lat = -12.0744;
    azs = 0:45:135;
    npts = 401;
    width = 20000;
    % half length in degrees
    halflength = width / 2 * m2deg;
    
    figure(11)
    clf
    set(gcf, 'Units', 'inches', 'Position', [18 8 6 8])
    % bathymetry
    ax1 = subplot(2,1,1);
    hold on

    % orientation
    ax2 = subplot(2,1,2);
    [~, ~, ~, ~, c, xoffset] = bathymetry([], [-0.09 0.09] + lon, ...
        [-0.09 0.09] + lat, true, ax2);
    c.Label.String = 'elevation (m)';
    c.Label.FontSize = 12;
    hold on
    colormap(ax2, flip(kelicol, 1));
    ax2.CLim = [-4900 -4600];
    
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
    set(ax1, 'TickDir', 'both', 'FontSize', 12)
    legend(ax1, '180 -- 0', '225 -- 45', '270 -- 90', '315 -- 135', ...
        'Location', 'northwest')
    
    axes(ax2)
    xlabel('longitude (degrees)')
    ylabel('latitude (degrees)')
    set(ax2, 'TickDir', 'both', 'FontSize', 12)
    
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

% request grided elevation
[longrid, latgrid, elev] = bathymetry([], ...
    [min(lons) max(lons)] + [-0.1 0.1], ...
    [min(lats) max(lats)] + [-0.1 0.1], ...
    false);

% convert longitude to [-180 180]
longrid = mod(longrid + 180, 360) - 180;

% get the elevation along the straight path
z = interp2(latgrid, longrid, double(elev), lats, lons);

% shift x-cooridnate to [0 width]
x = x - x(1);
end