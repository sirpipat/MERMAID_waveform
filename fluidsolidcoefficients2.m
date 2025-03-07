function [E, D, P] = fluidsolidcoefficients2(rho_f, vp_f, rho_s, vp_s, vs_s, theta)
% [E, D, P] = FLUIDSOLIDCOEFFCICIENTS2(rho_f, vp_f, rho_s, vp_s, vs_s, theta)
%
% Calculate reflection coefficient and transmission coefficients for P-wave
% and SV-wave given the properties of the media and the incident angle of
% the incoming P-wave in the elastic, solid medium.
%
% INPUT:
% rho_f     density of the fluid
% vp_f      P-wave velocity in the fluid
% rho_s     density of the solid
% vp_s      P-wave velocity in the solid
% vs_s      S-wave velocity in the solid
% theta     incident angle (in radians)
%
% OUTPUT:
% E         energy flux coefficient             : E = [Erp Ersv Et] where
%               Erp       reflection coefficient for P-wave
%               Ersv      reflection coefficient for SV-wave
%               Et        transmission coefficient
% D         displacement amplitude coefficient  : D = [Drp Drsv Dt] where
%               Drp       reflection coefficient for P-wave
%               Drsv      reflection coefficient for SV-wave
%               Dt        transmission coefficient
% P         potential amplitude coefficient     : P = [Prp Prsv Pt] where
%               Prp       reflection coefficient for P-wave
%               Prsv      reflection coefficient for SV-wave
%               Pt        transmission coefficient
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
% [E, D, P] = fluidsolidcoefficients2(rho_f, vp_f, rho_s, vp_s, vs_s, theta)
%
% % run a demo from Stein & Wysession (2003) chapter 2 page 83
% fluidsolidcoefficients2('demo1');
%
% % run a demo for a P-wave with a ray parameter of 600 rad s and various
% % P-wave velocity in the crust
% fluidsolidcoefficients('demo2', 600);
%
% SEE ALSO:
% FLUIDSOLIDCOEFFCICIENTS
%
% Last modified by sirawich-at-prineton.edu: 03/04/2025

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
        [E(ii,:),D(ii,:),P(ii,:)] = fluidsolidcoefficients2(rho_f, ...
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
    legend('Reflected P-wave', 'Reflected SV-wave', 'Transmitted P-wave', 'Location', 'best')
    ax2 = subplot(2,1,2);
    plot(theta, E(:,1) + E(:,2) + E(:,3), 'LineWidth', 1, 'Color', [0 0 0])
    xlabel('angle of incidence (^{\circ})')
    ylabel('Sum of energy flux ratio')
    xlim([0 90])
    ylim([-0.1 1.1])
    grid on
    return
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
        [E(ii,:),D(ii,:),P(ii,:)] = fluidsolidcoefficients2(rho_f, ...
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
    legend('Reflected P-wave', 'Reflected SV-wave', 'Transmitted P-wave', 'Location', 'best')
    
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
    legend('Reflected P-wave', 'Reflected SV-wave', 'Transmitted P-wave', 'Location', 'best')
    return
end

% ray parameter
p = sin(theta) / vp_s;

% intermediate variables for simplification
a = sqrt(1 / vp_s^2 - p^2);
b = p;
c = sqrt(1 / vp_f^2 - p^2);
d = rho_s * (1 - 2 * vs_s^2 * p^2);
e = 2 * rho_s * vs_s^2 * p * sqrt(1 / vs_s^2 - p^2);
f = rho_f;
g = 2 * p * a;
h = (1 / vs_s^2) - (2 * p^2);

U = (a * f * h) + (b * f * g) + (c * e * g);
V = c * d * h;

% potential amplitude ratio
rp_potential_amp_ratio = (U - V) / (U + V);
rsv_potential_amp_ratio = (2 * c * g * d) / (U + V);
t_potential_amp_ratio = 2 * d * (a * h + b * g) / (U + V);

% K = [a -b c; -d e f; g h 0] \ [a; d; g];
% rp_potential_amp_ratio = K(1);
% rsv_potential_amp_ratio = K(2);
% t_potential_amp_ratio = K(3);

P = [rp_potential_amp_ratio, rsv_potential_amp_ratio, t_potential_amp_ratio];

% displacement amplitude ratio
rp_amp_ratio = rp_potential_amp_ratio;
rsv_amp_ratio = rsv_potential_amp_ratio * (vp_s / vs_s) * isreal(a);
t_amp_ratio = t_potential_amp_ratio * (vp_s / vp_f) * isreal(c);

D = [rp_amp_ratio, rsv_amp_ratio, t_amp_ratio];

% energy coefficient
rp = abs(rp_potential_amp_ratio) ^ 2;
rsv = abs(rsv_potential_amp_ratio) ^ 2 * (vp_s / vs_s) * ...
    real(sqrt((1 - vs_s^2 * p^2) / (1 - vp_s^2 * p^2)));
t = abs(t_potential_amp_ratio) ^ 2 * (rho_f * vp_s) / (rho_s * vp_f) * ...
    real(sqrt((1 - vp_f^2 * p^2) / (1 - vp_s^2 * p^2)));

E = [rp rsv t];
end