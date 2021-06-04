function plot_station(ddir, network, station)
% PLOT_STATION(ddir, network, station)
% plot seismograms of a station
% 
% INPUT
% ddir      Directory of the seismogram files
% network   Network code
% station   Station number

coms = {'X','Z'};
labels = {'a','b'};
figure(2);
for ii = 1:length(coms)
    filename = strcat(ddir, network, sprintf('.S%s.BX%s.semd', ...
                      num2str(station,'%4.4i'), coms{ii}));
    ax = subplot(2, 1, ii);
    plot_seismogram(filename, ax);
    [x,y] = norm2trueposition(ax, 0.02, 0.9);
    text(x, y, labels{ii}, 'FontSize', 12)
end
% save the figure
eps_filename = strcat(network,sprintf('.S%s',num2str(station,'%4.4i')),...
                      '.epsc');
pdf_filename = strcat(network,sprintf('.S%s',num2str(station,'%4.4i')),...
                      '.pdf');
print(eps_filename, '-depsc');
system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
end