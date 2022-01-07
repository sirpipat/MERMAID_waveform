function d = rightdistance(vp, dt, h, theta)
% d = RIGHTDISTANCE(vp, dt, h, theta)
%
% Computes the minimum distance for the authentic first P-wave arrives
% before the spurious wave generated at the end of source array line by dt.
%
% INPUT:
% vp        P-wave velocity in the medium
% dt        separation time
% h         depth of the source from the station
% theta     incident angle of the P-wave (in radians)
%
% OUTPUT:
% d         minimum distance
%
% EXAMPLE:
% % basic function call
% d = rightdistance(3400, 3, 4080, 3 * pi/180);
%
% % Demo 1: making d vs theta plot for the above example
% rightdistance('demo1');
%
% % Demo 2: making d/h vs theta with various b = vp*dt/h plot for more 
% % general purpose
% rightdistance('demo2');
%
% Last modified by sirawich-at-princeton.edu, 01/07/2022

if ischar(vp) && strcmp(vp, 'demo1')
    vp = 3400;
    dt = [2 3 4];
    h = 4080;
    theta = (0:0.01:90)';
    
    d = rightdistance(vp, dt, h, theta * pi / 180);
    
    figure(1)
    clf
    plot(theta, d(:,1), 'LineWidth', 1)
    hold on
    plot(theta, d(:,2), 'LineWidth', 1)
    plot(theta, d(:,3), 'LineWidth', 1)
    grid on
    xlabel('Angle of incidence (^\circ)')
    ylabel('d (m)')
    legend('\Delta t = 2 s', '\Delta t = 3 s', '\Delta t = 4 s')
    set(gca, 'FontSize', 14, 'TickDir', 'both')
    ax = gca;
    ax.YAxis.Exponent = 0;
    ylim([0 20000])
    return
elseif ischar(vp) && strcmp(vp, 'demo2')
    b = [1 5 10];
    h = 1;
    vp = 3400;
    dt = b * h / vp;
    theta = (0:0.01:90)';
    d = rightdistance(vp, dt, h, theta * pi / 180);
    figure(1)
    clf
    plot(theta, d(:,1), 'LineWidth', 1)
    hold on
    plot(theta, d(:,2), 'LineWidth', 1)
    plot(theta, d(:,3), 'LineWidth', 1)
    grid on
    xlabel('Angle of incidence (^\circ)')
    ylabel('d / h')
    ylim([0 12])
    legend('b = 1', 'b = 5', 'b = 10')    
    set(gca, 'FontSize', 14, 'TickDir', 'both')
    return
end

b = vp .* dt ./ h;
A = sqrt(b.^2 + 2 * b .* cos(theta));
C = b + cos(theta);

d = (A .* C - sin(theta)) ./ (C + sin(theta) .* A) .* h;
end