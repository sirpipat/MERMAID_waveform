function plot_source_time_function(stffile)
% PLOT_SOURCE_TIME_FUNCTION(stffile)
%
% Plots source time function.
%
% INPUT
% stffile       filename of the source time function

sizeData = [2 Inf];
fid = fopen(stffile,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);
figure
plot(data(1,:),data(2,:),'LineWidth',1,'Color','k')
xlim([-0.05 0.05])
grid on
xlabel('time [s]')
title('source time function')
end