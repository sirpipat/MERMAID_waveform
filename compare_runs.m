function compare_runs(ddir_run1, ddir_run2, wh)


% read stations
station_fname1 = strcat(ddir_run1, 'DATA/STATIONS');
[n1, name1, network1, ~, z1] = read_stations(station_fname1);
station_fname2 = strcat(ddir_run1, 'DATA/STATIONS');
[n2, name2, network2, ~, z2] = read_stations(station_fname2);

% setting parameters
% TODO: figure out how to read the parameters from SOURCE, Par_file etc.
interface_fname = strcat(ddir_run1,'DATA/interfaces_fluid_flat.dat');
source_depth = 4080;
[water_depth, solid_depth] = get_fluid_solid_setting(interface_fname);
vp = 3400;
vw = 1500;

% add expected arrivals of each phase
depth = water_depth + solid_depth - z1;
ta1 = source_depth / vp + (water_depth - depth) / vw;
ta2 = source_depth / vp + (water_depth + depth) / vw;

% number of reflections
num_reflects = 5;
t_add = reshape(repmat(0:num_reflects, 2, 1), [1 2*(num_reflects+1)]) * ...
    (2 * water_depth / vw);

% expected arrival times
tas =  repmat([ta1, ta2], 1, num_reflects+1) + t_add;

for jj = 1:n1
    filename = strcat(ddir_run1, network1{jj}, '.', name1{jj}, '.PRE.semp');
    labels = removepath(filename);
    % read seismograms
    [t1, x1] = read_seismogram(filename);
    
    % plot the signal
    figure(1)
    set(gcf, 'Unit', 'inches', 'Position', [1 10 7 5])
    clf
    plot(t1, x1, 'Color', 'k', 'LineWidth', 1.5);
    
    filename = strcat(ddir_run2, network2{jj}, '.', name2{jj}, '.PRE.semp');
    labels = removepath(filename);
    % read seismograms
    [t2, x2] = read_seismogram(filename);
    
    % plot the signal
    hold on
    plot(t2, x2, 'Color', [0.9 0.15 0.15], 'LineWidth', 0.9);
    hold off
    grid on
    xlabel('time (s)')
    ylabel('pressure (counts)')
    legend('reference', 'Munk')
    
    for ii = 1:size(tas, 2)
        figure(1)
        xlim([-0.5 0.5] + tas(jj, ii))
        title(sprintf('depth = %.2f m, arrival # %02d', depth(jj), ii))
        set(gca, 'TickDir', 'both', 'FontSize', 11);
        figdisp(sprintf('%s_pressure_depth_%04d_arrival_%02d', ...
            mfilename, depth(jj), ii), [], [], 2, 'epsc', 'epstopdf');

        % find timeshift
        where = and(t2 >= tas(jj, ii) - 0.5, t2 < tas(jj, ii) + 0.5);
        x1_where = x1(where);
        x2_where = x2(where);
        [c, lags] = xcov(x1_where, x2_where, 'coeff');
        dt = t2(2) - t1(1);

        % max timeshift
        [~, i_max] = max(c);
        t_shift = lags(i_max) * dt;

        % plot the timeshift
        figure(2)
        set(gcf, 'Unit', 'inches', 'Position', [8 10 7 5])
        clf
        plot(lags * dt, c, 'Color', 'k', 'LineWidth', 1);
        grid on
        ylim([-1 1])
        xlabel('time shift (s)')
        ylabel('correlation coefficient')
        title(sprintf('depth = %.2f m, arrival # %02d, timeshift = %.4f s', ...
            depth(jj), ii, t_shift))
        set(gca, 'TickDir', 'both', 'FontSize', 11);
        figdisp(sprintf('%s_cc_depth_%04d_arrival_%02d', mfilename, ...
            depth(jj), ii), [], [], 2, 'epsc', 'epstopdf');
        pause(1.5)
    end
end
end