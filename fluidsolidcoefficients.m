function [E, D, P] = fluidsolidcoefficients(rho_f, vp_f, rho_s, vp_s, vs_s, theta)
% [E, D, P] = fluidsolidcoefficients(rho_f, vp_f, rho_s, vp_s, vs_s, theta)
%
% Calculate reflection coefficient and transmission coefficients for P-wave
% and SV-wave given the properties of the media and the incident angle.
%
% INPUT:
% rho_f     density of the fluid
% vp_f      P-wave velocity in the fluid
% rho_s     density of the solid
% vp_s      P-wave velocity in the solid
% vs_s      S-wave velocity in the solid
% theta     incidence angle
%
% OUTPUT:
% E         energy flux coefficient             : E = [Er Etp Etsv] where
%               Er        reflection coefficient
%               Etp       transmission coefficient for P-wave
%               Etsv      transmission coefficient for SV-wave
% D         displacement amplitude coefficient  : D = [Dr Dtp Dtsv] where
%               Dr        reflection coefficient
%               Dtp       transmission coefficient for P-wave
%               Dtsv      transmission coefficient for SV-wave
% P         potential amplitude coefficient     : P = [Pr Ptp Ptsv] where
%               Pr        reflection coefficient
%               Ptp       transmission coefficient for P-wave
%               Ptsv      transmission coefficient for SV-wave
%
% EXAMPLE:
% % simple call
% rho_f = 1000;
% vp_f = 1500;
% rho_s = 3000;
% vp_s = 5000;
% vs_s = 3000;
% theta = 10 * pi/180;  % 10 degrees in radians
%
% [E, D, P] = fluidsolidcoefficients(rho_f, vp_f, rho_s, vp_s, vs_s, theta)
%
% % run a demo from Stein & Wysession (2003) chapter 2 page 83
% fluidsolidcoefficients('demo1');
%
% % run a demo for a P-wave with a ray parameter of 600 rad s and various
% % P-wave velocity in the crust
% fluidsolidcoefficients('demo2', 600);
%
% % run a demo for varying P-wave ray parameters with a realistic ocean 
% % moodel. You can specify the P-wave velocity in the solid in the second 
% % input
% fluidsolidcoefficients('demo3', 5000);
%
% Last modified by sirawich-at-princeton.edu, 04/03/2025

% demo
% run a demo from Stein & Wysession (2003) chapter 2 page 83
if ischar(rho_f) && strcmp(rho_f, 'demo1')
    rho_f = 1000;
    vp_f = 1500;
    rho_s = 3000;
    vp_s = 5000;
    vs_s = 3000;
    
    theta = (0:0.01:90)';
    P = zeros(size(theta, 1), 3);
    D = zeros(size(theta, 1), 3);
    E = zeros(size(theta, 1), 3);
    
    for ii = 1:size(theta,1)
        [E(ii,:),D(ii,:),P(ii,:)] = fluidsolidcoefficients(rho_f, ...
            vp_f, rho_s, vp_s, vs_s, theta(ii,1) * pi/180);
    end
    
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [18 2 6 8])
    ax1 = subplot(2,1,1);
    plot(theta, E(:,1), 'LineWidth', 1, 'Color', [0.1 0.8 0.1])
    hold on
    plot(theta, E(:,2), 'LineWidth', 1, 'Color', [0.1 0.1 0.8])
    plot(theta, E(:,3), 'LineWidth', 1, 'Color', [0.8 0.1 0.1])
    grid on
    xlabel('angle of incidence (\circ)')
    ylabel('Energy flux ratio')
    xlim([0 90])
    ylim([-0.1 1.1])
    legend('Reflected P-wave', 'Transmitted P-wave', 'Transmitted S-wave', 'Location', 'best')
    ax2 = subplot(2,1,2);
    plot(theta, E(:,1) + E(:,2) + E(:,3), 'LineWidth', 1, 'Color', [0 0 0])
    xlabel('angle of incidence (^{\circ})')
    ylabel('Sum of energy flux ratio')
    xlim([0 90])
    ylim([-0.1 1.1])
    grid on
    return
% run a demo for a P-wave with a ray parameter of 600 rad s and various
% P-wave velocity in the crust
elseif ischar(rho_f) && strcmp(rho_f, 'demo2')
    defval('vp_f', 600);
    p = vp_f;
    r = 6365000;
    rho_f = 1020;
    vp_f = 1500;
    rho_s = 2500;
    vp_s = (2000:10:min(10000, r/p))';
    vs_s = vp_s / sqrt(3);
    
    theta = asin(p/r * vp_f);
    P = zeros(size(vp_s, 1), 3);
    D = zeros(size(vp_s, 1), 3);
    E = zeros(size(vp_s, 1), 3);
    
    for ii = 1:size(vp_s,1)
        [E(ii,:),D(ii,:),P(ii,:)] = fluidsolidcoefficients(rho_f, ...
            vp_f, rho_s, vp_s(ii,1), vs_s(ii,1), theta);
    end
    
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [18 2 6 8])
    ax1 = subplot(2,1,1);
    plot(vp_s, E(:,1), 'LineWidth', 1, 'Color', [0.1 0.8 0.1])
    hold on
    plot(vp_s, E(:,2), 'LineWidth', 1, 'Color', [0.1 0.1 0.8])
    plot(vp_s, E(:,3), 'LineWidth', 1, 'Color', [0.8 0.1 0.1])
    grid on
    xlabel('P-wave velocity in solid (m/s)')
    ylabel('Energy flux ratio')
    title(sprintf('ray parameter = %.2f rad s', p))
    xlim([2000 10000])
    ylim([-0.1 1.1])
    legend('Reflected P-wave', 'Transmitted P-wave', 'Transmitted S-wave', 'Location', 'best')
    
    ax2 = subplot(2,1,2);
    plot(vp_s, D(:,1), 'LineWidth', 1, 'Color', [0.1 0.8 0.1])
    hold on
    plot(vp_s, D(:,2), 'LineWidth', 1, 'Color', [0.1 0.1 0.8])
    plot(vp_s, D(:,3), 'LineWidth', 1, 'Color', [0.8 0.1 0.1])
    grid on
    xlabel('P-wave velocity in solid (m/s)')
    ylabel('Displacement amplitude ratio')
    title(sprintf('ray parameter = %.2f rad s', p))
    xlim([2000 10000])
    legend('Reflected P-wave', 'Transmitted P-wave', 'Transmitted S-wave', 'Location', 'best')
    return
% run a demo for varying P-wave ray parameters 
elseif ischar(rho_f) && strcmp(rho_f, 'demo3')
    p = (0:10:4000)';
    r = 6365000;
    defval('vp_f', 5000)
    vp_s = vp_f;
    
    rho_f = 1020;
    vp_f = 1500;
    rho_s = 3340;
    vs_s = 4500;
    
    theta = asin(p/r * vp_f);
    P = zeros(size(p, 1), 3);
    D = zeros(size(p, 1), 3);
    E = zeros(size(p, 1), 3);
    
    for ii = 1:size(p,1)
        [E(ii,:),D(ii,:),P(ii,:)] = fluidsolidcoefficients(rho_f, ...
            vp_f, rho_s, vp_s, vs_s, theta(ii,1));
    end
    
    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [18 2 6 8])
    ax1 = subplot(2,1,1);
    plot(p, E(:,1), 'LineWidth', 1, 'Color', [0.1 0.8 0.1])
    hold on
    plot(p, E(:,2), 'LineWidth', 1, 'Color', [0.1 0.1 0.8])
    plot(p, E(:,3), 'LineWidth', 1, 'Color', [0.8 0.1 0.1])
    grid on
    xlabel('Ray parameter (rad s)')
    ylabel('Energy flux ratio')
    title(sprintf('\x03c1_s = %.2f kg/m^3, \x03b1_s = %.2f m/s, \x03b2_s = %.2f m/s', rho_s, vp_s, vs_s))
    xlim([0 4000])
    ylim([-0.1 1.1])
    legend('Reflected P-wave', 'Transmitted P-wave', 'Transmitted S-wave', 'Location', 'best')
    ax2 = subplot(2,1,2);
    plot(p, abs(D(:,1)), 'LineWidth', 1, 'Color', [0.1 0.8 0.1])
    hold on
    plot(p, abs(D(:,2)), 'LineWidth', 1, 'Color', [0.1 0.1 0.8])
    plot(p, abs(D(:,3)), 'LineWidth', 1, 'Color', [0.8 0.1 0.1])
    grid on
    xlabel('Ray parameter (rad s)')
    ylabel('Displacement amplitude ratio')
    title(sprintf('\x03c1_s = %.2f kg/m^3, \x03b1_s = %.2f m/s, \x03b2_s = %.2f m/s', rho_s, vp_s, vs_s))
    xlim([0 4000])
    legend('Reflected P-wave', 'Transmitted P-wave', 'Transmitted S-wave', 'Location', 'best')
    return
end

% incident angle
theta_i = theta;

% transmitted angle
theta_t = asin(vp_s / vp_f * sin(theta_i));
theta_s = asin(vs_s / vp_f * sin(theta_i));

% intermediate variables for simplification
a = cos(theta_i) / vp_f;
b = cos(theta_t) / vp_s;
c = sin(theta_s) / vs_s;
d = rho_f;
f = rho_s * cos(2 * theta_s);
g = rho_s * sin(2 * theta_s);
m = sin(2 * theta_t) / vp_s^2;
n = cos(2 * theta_s) / vs_s^2;
p = a * (f * n + g * m);
q = d * (b * n + c * m);

% potential amplitude ratio
r_potential_amp_ratio = (p - q) / (p + q);
tp_potential_amp_ratio = 2 * a * d * n / (p + q);
tsv_potential_amp_ratio = -2 * a * d * m / (p + q);

P = [r_potential_amp_ratio, tp_potential_amp_ratio, tsv_potential_amp_ratio];

% displacement amplitude ratio
r_amp_ratio = r_potential_amp_ratio;
tp_amp_ratio = tp_potential_amp_ratio * (vp_f / vp_s) * isreal(theta_t);
tsv_amp_ratio = tsv_potential_amp_ratio * (vp_f / vs_s) * isreal(theta_s);

D = [r_amp_ratio, tp_amp_ratio, tsv_amp_ratio];

% energy coefficient
r = abs(r_amp_ratio) ^ 2;
tp = abs(tp_amp_ratio) ^ 2 * real(rho_s * vp_s * cos(theta_t)) / ...
    (rho_f * vp_f * cos(theta_i));
tsv = abs(tsv_amp_ratio) ^ 2 * real(rho_s * vs_s * cos(theta_s)) / ...
    (rho_f * vp_f * cos(theta_i));

E = [r, tp, tsv];
end