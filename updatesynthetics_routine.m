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
% Last modified by sirawich-at-princeton.edu, 02/16/2024

defval('synmasterdir', '/Users/sirawich/research/processed_data/SYNTHETICS/')
defval('plt', false)

[allsynfiles, sndex] = allfilen(synmasterdir, 2);

for ii = 1:sndex
    updatesynthetics(allsynfiles{ii}, 'ak135');
    if plt
        [SeisData, HdrData] = readsac(allsynfiles{ii});
        plotsac2(SeisData, HdrData, 'Color', 'k');
        savename = sprintf('%s_%d_%s.eps', mfilename, HdrData.USER7, ...
            HdrData.KSTNM(ismember(HdrData.KSTNM, 33:116)));
        figdisp(savename,[],[],2,[],'epstopdf');
    end
end
end