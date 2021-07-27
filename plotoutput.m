function plotoutput(ddir)
% PLOTOUTPUT(ddir)
%
% Makes various plots from the output. Then, saves to ddir/PLOTS/.
%
% INPUT:
% ddir          simulation directory
%
% Last modified by sirawich@princeton.edu, 07/27/2021

example = removepath(ddir(1:end-1));

% makes a directory for the plots if there is not one
savedir = [ddir 'PLOTS/'];
system(['mkdir ' savedir]);

% plot setting
drawsetting(ddir, savedir, example, true);

% reads receiversets
[~, ~, networks, ~, ~] = read_stations([ddir 'DATA/STATIONS']);
networks = unique(networks);

for ii = 1:length(networks)
    plotarray(ddir, networks{ii}, savedir, networks{ii}, true);
end
end