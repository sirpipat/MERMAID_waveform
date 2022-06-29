function figure_adaptivefreq_CCmax()
% figure_adaptivefreq_CCmax()
%
% Plot histograms of maximum correlation coefficients from correlation
% travel time measurement of the 2 cases:
% 1. using the same frequency band between 1--2 Hz
% 2. using adaptive frequency bands chosen by FREQSELECT
%
% It use the output files produced by BATHYMATTER saved in $MERMAID2/DATA/.
%
% SEE ALSO:
% FREQSELECT, BATHYMATTER, COMPARERESPONSEFUNCTIONS
%
% Last modified by sirawich-at-princeton.edu, 06/29/2022

fname_fixed = sprintf('%sDATA/bathymatter_fixed.mat', ...
    getenv('MERMAID2'));
fname_adaptive = sprintf('%sDATA/bathymatter_adaptive.mat', ...
    getenv('MERMAID2'));

fixed = load(fname_fixed);
adaptive = load(fname_adaptive);

figure(1)
set(gcf, 'Units', 'inches', 'Position', [0 4 12 6]);

subplot('Position',[0.06 0.1 0.4 0.8])
histogram(fixed.CCmaxs(:,2),'FaceAlpha',0.4,'FaceColor',rgbcolor('1'))
grid on
xlim([0 1])
ylim([0 160])
xticks(0:0.1:1)
xlabel('Correlation coefficients')
ylabel('Counts')
set(gca, 'TickDir', 'out', 'FontSize', 11)
title('Fixed frequency band 1--2 Hz', 'FontSize', 14)
vline(gca, median(fixed.CCmaxs(:,2)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-.');

subplot('Position',[0.56 0.1 0.4 0.8])
histogram(adaptive.CCmaxs(:,2),'FaceAlpha',0.4,'FaceColor',rgbcolor('1'))
grid on
xlim([0 1])
ylim([0 160])
xticks(0:0.1:1)
xlabel('Correlation coefficients')
ylabel('Counts')
set(gca, 'TickDir', 'out', 'FontSize', 11)
title('Adaptive frequency bands', 'FontSize', 14)
vline(gca, median(adaptive.CCmaxs(:,2)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-.');

set(gcf, 'Renderer', 'painters')
figdisp('CCmaxs_adaptive_vs_flat.pdf',[],[],2,[],'epstopdf');
end