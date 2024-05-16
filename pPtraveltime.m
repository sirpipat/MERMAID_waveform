function pPtraveltime
% PPTRAVELTIME
%
% Demonstrates how pP arrivals are close to the first P arrival if the 
% source is shallow enough. Plots t_pP - tP curve at a few source depths.
%
% Last modified by sirawich-at-princeton.edu: 09/21/2023

depth = [10 30 50 70 100 150 200];
distance = (0:0.25:100)';

[ddepth, ddistance] = meshgrid(depth, distance);
tP = nan(size(ddepth));
tpP = nan(size(ddepth));

for ii = 1:size(ddepth,1)
    for jj = 1:size(ddepth,2)
        try
            tP(ii,jj) = getfield(indeks(tauptime('dep',ddepth(ii,jj),...
                'ph','P','deg',ddistance(ii,jj),'mod','ak135'),1),'time');
        catch
        end
        try
            tpP(ii,jj) = getfield(indeks(tauptime('dep',ddepth(ii,jj),...
                'ph','pP','deg',ddistance(ii,jj),'mod','ak135'),1),'time');
        catch
        end
    end
end

figure
set(gcf, 'Units', 'inches', 'Position', [0 2 8 5]);
hold on
for ii = 1:size(ddepth,2)
    plot(ddistance(:,ii), tpP(:,ii) - tP(:,ii), 'LineWidth', 1);
end
grid on
xlabel('epicentral distance (degrees)')
ylabel('time since first P-wave arrival (seconds)')
title('pP wave arrival time')
legend('10 km', '30 km', '50 km', '70 km', '100 km', '150 km', ...
    '200 km', 'Location', 'northwest')
set(gca, 'TickDir', 'out', 'FontSize', 12, 'Box', 'on')
end