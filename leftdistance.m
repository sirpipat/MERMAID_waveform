function [x1, x2] = leftdistance(N, theta, xs, H, d)
% [x1, x2] = LEFTDISTANCE(N, theta, xs, H, d)
%
% Computes the minimum distance for P-wave to remain planar when it arrives
% the station after reflecting at the ocean bottom N times.
%
% INPUT:
% N         number of reflections at the ocean bottom
% theta     angle of incidence (in radians)
% xs        distance from left boundary to first point of contact
% H         thickness of the ocean
% d         station depth
%
% OUTPUT:
% x1        minimum distance from left boundary for upgoing wave
% x2        minimum distance form left boundary for downgoing wave
%
% EXAMPLE:
% % basic function call
% [x1, x2] = leftdistance(7, 3 * pi/180, 2000, 6000, 1500);
%
% % Demo 1: plotting the output above
% leftdistance('demo1');
%
% % Demo 2: plotting the minimum distance for N=0 to 10
% leftdistance('demo2');
% 
% Last modified by sirawich-at-princeton.edu, 01/07/2022

if ischar(N) && strcmp(N, 'demo1')
    N = 7;
    theta = (0:0.01:40)';    
    xs = 2000;
    H = 6000;
    d = 1500;
    [x1, x2] = leftdistance(N, theta * pi/180, xs, H, d);
    
    % ray parameter
    p = (6371000 - H) * sin(theta * pi/180) / 1500;
    
    figure(1)
    clf
    plot(theta, x1, 'LineWidth', 1)
    hold on
    plot(theta, x2, 'LineWidth', 1)
    grid on
    xlabel('Angle of incidence (^\circ)')
    ylabel('Minimum distance (m)')
    legend('upgoing', 'downgoing', 'Location', 'best')
    set(gca, 'FontSize', 14, 'TickDir', 'both')
    ax = gca;
    ax.YAxis.Exponent = 0;
    ylim([0 80000])
    
    figure(4)
    clf
    plot(p, x1, 'LineWidth', 1)
    hold on
    plot(p, x2, 'LineWidth', 1)
    grid on
    xlabel('Ray parameter (rad s)')
    ylabel('Minimum distance (m)')
    legend('upgoing', 'downgoing', 'Location', 'best')
    set(gca, 'FontSize', 14, 'TickDir', 'both')
    ax = gca;
    ax.YAxis.Exponent = 0;
    ylim([0 60000])
    xlim([0 2000])
    return
elseif ischar(N) && strcmp(N, 'demo2')
    N = 0:10;
    defval('theta', 6000);
    H = theta;
    theta = (0:0.01:40)';    
    xs = 2000;
    d = 1500;
    [~, x2] = leftdistance(N, theta * pi/180, xs, H, d);
    
    % ray parameter
    p = (6371000 - H) * sin(theta * pi/180) / 1500;
    
    figure(1)
    clf
    plot(theta, x2(:,1), 'LineWidth', 1)
    hold on
    plot(theta, x2(:,2), 'LineWidth', 1)
    plot(theta, x2(:,3), 'LineWidth', 1)
    plot(theta, x2(:,4), 'LineWidth', 1)
    plot(theta, x2(:,5), 'LineWidth', 1)
    plot(theta, x2(:,6), 'LineWidth', 1)
    plot(theta, x2(:,7), 'LineWidth', 1)
    plot(theta, x2(:,8), 'LineWidth', 1, 'LineStyle', '--')
    plot(theta, x2(:,9), 'LineWidth', 1, 'LineStyle', '--')
    plot(theta, x2(:,10), 'LineWidth', 1, 'LineStyle', '--')
    plot(theta, x2(:,11), 'LineWidth', 1, 'LineStyle', '--')
    grid on
    xlabel('Angle of incidence (^\circ)')
    ylabel('Minimum distance (m)')
    legend('N=0', 'N=1', 'N=2', 'N=3', 'N=4', 'N=5', 'N=6', 'N=7', ...
        'N=8', 'N=9', 'N=10', 'Location', 'best')
    set(gca, 'FontSize', 14, 'TickDir', 'both')
    ax = gca;
    ax.YAxis.Exponent = 3;
    %ylim([0 100000])
    
    figure(4)
    clf
    plot(p, x2(:,1), 'LineWidth', 1)
    hold on
    plot(p, x2(:,2), 'LineWidth', 1)
    plot(p, x2(:,3), 'LineWidth', 1)
    plot(p, x2(:,4), 'LineWidth', 1)
    plot(p, x2(:,5), 'LineWidth', 1)
    plot(p, x2(:,6), 'LineWidth', 1)
    plot(p, x2(:,7), 'LineWidth', 1)
    plot(p, x2(:,8), 'LineWidth', 1, 'LineStyle', '--')
    plot(p, x2(:,9), 'LineWidth', 1, 'LineStyle', '--')
    plot(p, x2(:,10), 'LineWidth', 1, 'LineStyle', '--')
    plot(p, x2(:,11), 'LineWidth', 1, 'LineStyle', '--')
    grid on
    xlabel('Ray parameter (rad s)')
    ylabel('Minimum distance (m)')
    legend('N=0', 'N=1', 'N=2', 'N=3', 'N=4', 'N=5', 'N=6', 'N=7', ...
        'N=8', 'N=9', 'N=10', 'Location', 'best')
    set(gca, 'FontSize', 14, 'TickDir', 'both')
    ax = gca;
    ax.YAxis.Exponent = 3;
    %ylim([0 100000])
    xlim([0 2000])
end

x1 = xs + ((2 * N + 1) .* H - d) .* tan(theta);
x2 = xs + ((2 * N + 1) .* H + d) .* tan(theta);
end