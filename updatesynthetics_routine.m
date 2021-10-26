function updatesynthetics_routine(synmasterdir, plt)
% UPDATESYNTHETICS_ROUTINE(synmasterdir, plt)
%
% A script to update headers of the synthetic files stored in synmasterdir
% and then plot the seismograms with expected arrival times.
%
% INPUT:
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% plt               whether to plot or not
%
% SEE ALSO:
% UPDATESYNTHETICS
%
% Last modified by sirawich-at-princeton.edu, 10/26/2021

defval('synmasterdir', '/Users/sirawich/research/processed_data/SYNTHETICS/')
defval('plt', false)

[allsyndirs, dndex] = allfile(synmasterdir);

for ii = 1:dndex
    [allsynfiles, sndex] = allfile([allsyndirs{ii} '/']);
    for jj = 1:sndex
        updatesynthetics(allsynfiles{jj}, 'ak135');
        if plt
            figure(1);
            clf
            [SeisData, HdrData] = readsac(allsynfiles{jj});
            ax = plotsac(SeisData, HdrData, gca, 'Color', 'k');
            ax.Box = 'on';
            savename = sprintf('plotsynthetic_%d_%s.eps', ...
                HdrData.USER7, replace(HdrData.KSTNM, ' ', ''));
            figdisp(savename,[],[],2,[],'epstopdf');
        end
    end
end
end