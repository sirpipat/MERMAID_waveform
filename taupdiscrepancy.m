function taupdiscrepancy(depth, phase)
% taupdiscrepancy(depth, phase)
% 
% Illustrate the arrival time and ray parameter discrepancies between PREM,
% IASP91 and AK135 models.
%
% INPUT:
% depth         source depth in km      [default: 50]
% pahse         seismic phase name      [default: 'P']
%
% Last modified by sirawich-at-princeton.edu, 03/09/2022

defval('depth', 50)
defval('phase', 'P')

% epicentral distance in degrees
deg = 4:99;

% travel time
t_prem = zeros(size(deg));
t_iasp91 = zeros(size(deg));
t_ak135 = zeros(size(deg));

% ray parameter
p_prem = zeros(size(deg));
p_iasp91 = zeros(size(deg));
p_ak135 = zeros(size(deg));

% compute the travel times and ray parameters
for ii = 1:length(deg)
    tt = tauptime('mod', 'prem', 'dep', depth, 'ph', phase,'deg',deg(ii));  
    if length(tt) <= 0
        % set to NaN if there is no arrival
        t_prem(ii) = nan;
        p_prem(ii) = nan;            
    else
        % use only the first arrival
        t_prem(ii) = tt(1).time;
        % convert the output to rad s
        p_prem(ii) = tt(1).rayparameter * 180 / pi;
    end
    tt = tauptime('mod','iasp91', 'dep', depth,'ph', phase,'deg',deg(ii));
    if length(tt) <= 0
        % set to NaN if there is no arrival
        t_iasp91(ii) = nan;    
        p_iasp91(ii) = nan;            
    else
        % use only the first arrival
        t_iasp91(ii) = tt(1).time;
        % convert the output to rad s
        p_iasp91(ii) = tt(1).rayparameter * 180 / pi;
    end
    tt = tauptime('mod','ak135', 'dep', depth,'ph', phase,'deg',deg(ii)); 
    if length(tt) <= 0
        % set to NaN if there is no arrival
        t_ak135(ii) = nan;    
        p_ak135(ii) = nan;            
    else
        % use only the first arrival
        t_ak135(ii) = tt(1).time;
        % convert the output to rad s
        p_ak135(ii) = tt(1).rayparameter * 180 / pi;
    end
end

%% travel time discrepencies plot
figure(10)
set(gcf, 'Units', 'inches', 'Position', [18 8 6 8]);
clf

% plot the travel time according to each model
subplot('Position', [0.08 0.58 0.8 0.36])
cla
plot(deg, t_prem, 'LineWidth', 1)
hold on
plot(deg, t_iasp91, 'LineWidth', 1)
plot(deg, t_ak135, 'LineWidth', 1)
grid on 
legend('PREM', 'IASP91', 'AK135', 'Location', 'northwest')
title(sprintf('travel time, phase = ''%s'', source depth = %.2f km', phase, depth))
xlabel('epicentral distance (degrees)')
ylabel('travel time (s)')

% plot the travel time discrepencies between the models
subplot('Position', [0.08 0.08 0.8 0.36])
cla
plot(deg, t_iasp91 - t_prem, 'LineWidth', 1)
hold on
plot(deg, t_ak135 - t_iasp91, 'LineWidth', 1)
grid on
legend('IASP91 - PREM', 'AK135 - IASP91', 'Location', 'best')
title(sprintf('travel time difference, phase = ''%s'', source depth = %.2f km', phase, depth))
xlabel('epicentral distance (degrees)')
ylabel('travel time difference (s)')

% save the figure
set(gcf, 'Renderer', 'painters')
savename = sprintf('%s_depth-%d_phase-%s_time.eps', mfilename, depth, phase);
figdisp(savename, [], [], 2, [], 'epstopdf');

%% ray parameter discrepencies
figure(11)
set(gcf, 'Units', 'inches', 'Position', [12 8 6 8]);
clf

% plot the ray parameter for each model
subplot('Position', [0.08 0.58 0.8 0.36])
cla
plot(deg, p_prem, 'LineWidth', 1)
hold on
plot(deg, p_iasp91, 'LineWidth', 1)
plot(deg, p_ak135, 'LineWidth', 1)
grid on
legend('PREM', 'IASP91', 'AK135', 'Location', 'best')
title(sprintf('ray parameter, phase = ''%s'', source depth = %.2f km', phase, depth))
xlabel('epicentral distance (degrees)')
ylabel('ray parameter (rad s)')

% plot the ray parameter discrepencies between models
subplot('Position', [0.08 0.08 0.8 0.36])
cla
plot(deg, p_iasp91 - p_prem, 'LineWidth', 1)
hold on
plot(deg, p_ak135 - p_iasp91, 'LineWidth', 1)
grid on
legend('IASP91 - PREM', 'AK135 - IASP91', 'Location', 'best')
title(sprintf('ray parameter difference, phase = ''%s'', source depth = %.2f km', phase, depth))
xlabel('epicentral distance (degrees)')
ylabel('ray parameter difference (rad s)')

% save the figure
set(gcf, 'Renderer', 'painters')
savename = sprintf('%s_depth-%d_phase-%s_rayparam.eps', mfilename, depth, phase);
figdisp(savename, [], [], 2, [], 'epstopdf');
end