function [r, a] = presiduestat(sacfiles, fcs, plt)
% [r, a] = presiduestat(sacfiles, fcs, plt)
% 
% Compute the residual times between the pick arrvial and the assigned
% arrival time for P-phase. The arivals are picked by finding the first
% incident when the absolute signal exceeds 2 percent of the maximum of the
% absolute of the signal within 30 seconds centered at the assigned
% arrival.
%
% INPUT:
% sacfiles      cell array to SAC files
% fcs           corner frequency pairs [default: []]
%               - empty, no filtering
%               - 1x2 matrix, apply the same bandpass to all waveforms
%               - Nx2 matrix where N is the number of sacfiles, i-th
%                 waveform is filtered using the i-th corner frequency pair
% plt           whether to plot the picks and the stat or not 
%               [default: false]
%
% OUTPUT:
% r             residuals
% a             amplitude (signed absolute maximum)
%
% Last modified by sirawich-at-princeton.edu, 11/26/2023

defval('fcs', [])
defval('plt', false)

if all(size(fcs) == [1 2])
    fcs = repmat(fcs, length(sacfiles), 1);
end

n = length(sacfiles);
r = zeros(length(sacfiles), 1);
a = zeros(length(sacfiles), 1);

for ii = 1:n
    % read the sac file
    [x, HdrData] = readsac(sacfiles{ii});
    
    % gets the information from SAC header
    [dt_ref, dt_B, dt_E, fs, npts, dts, tims] = gethdrinfo(HdrData);
    t = seconds(tims - tims(1));
    
    % detect P-phase
    phaseNum = 0;
    phaseTime = HdrData.T0;
    phaseName = HdrData.KT0;
    
    while (~strcmpi(phaseName, 'P') && ~strcmpi(phaseName, 'PKP') && ...
            ~strcmp(phaseName, 'PKIKP') && ...
            ~strcmpi(phaseName, 'Pdiff')) && phaseTime == -12345 && ...
            phaseNum < 10
        phaseNum = phaseNum + 1;
        switch phaseNum
            case 1
                phaseTime = HdrData.T1;
                phaseName = HdrData.KT1;
            case 2
                phaseTime = HdrData.T2;
                phaseName = HdrData.KT2;
            case 3
                phaseTime = HdrData.T3;
                phaseName = HdrData.KT3;
            case 4
                phaseTime = HdrData.T4;
                phaseName = HdrData.KT4;
            case 5
                phaseTime = HdrData.T5;
                phaseName = HdrData.KT5;
            case 6
                phaseTime = HdrData.T6;
                phaseName = HdrData.KT6;
            case 7
                phaseTime = HdrData.T7;
                phaseName = HdrData.KT7;
            case 8
                phaseTime = HdrData.T8;
                phaseName = HdrData.KT8;
            case 9
                phaseTime = HdrData.T9;
                phaseName = HdrData.KT9;
            otherwise
                break
        end
    end
    
    % filter the waveform
    if ~isempty(fcs)
        x = bandpass(x, fs, fcs(ii, 1), fcs(ii, 2), 4, 2, ...
            'butter', 'linear');
    end
    
    % pick the arrival based on the rise of the signal
    wh = and(t - phaseTime > -15, t - phaseTime < 15);
    t_wh = t(wh);
    x_wh = x(wh);
    
    [absmax, whmax] = max(abs(x_wh));
    
    arrival = indeks(t_wh(abs(x_wh) > 2e-2 * absmax), 1);
    r(ii) = arrival - phaseTime;
    a(ii) = x_wh(whmax);
    
    if plt
        if HdrData.USER7 == -12345
            HdrData.USER7 = ii;
        end
        plotsac2(x, HdrData, 'Color', 'k');
        
        % add the pick arrival
        fig = gcf;
        axb = fig.Children(3);
        ax = fig.Children(4);
        hold on
        vline(ax, r(ii), 'LineWidth', 1, 'LineStyle', '-', 'Color', [0.1 0.4 0.9]);
        hold off
                
        % move the seismogram to the front
        ax.Children = ax.Children([2:end 1]);
        
        % add legend
        legend(ax.Children(3:-1:2), 'Instaseis', 'TauP', 'Location', 'southwest')

        % bring back plot label
        axes(axb)
        
        savename = sprintf('%s_seis_%d_%s.eps', mfilename, ...
            HdrData.USER7, replace(HdrData.KSTNM, ' ', ''));
        figdisp(savename,[],[],2,[],'epstopdf');
        delete(fig)
    end
end

if plt
    figure
    histogram(r, 'BinWidth', 0.1, 'FaceColor', [0.75 0.75 0.75])
    hold on
    xlim([-5, 1])
    [~,v1] = vline(gca, mean(r), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
    [~,v2] = vline(gca, median(r), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-');
    [~,v3] = vline(gca, median(r) + std(r) * [-1 1], 'Color', [0.1 0.4 0.9], 'LineWidth', 1.5, 'LineStyle', '-');
    grid on
    ax = gca;
    ax.Children = ax.Children([end 1:(end-1)]);
    set(ax, 'FontSize', 12, 'TickDir', 'out');
    xlabel('residual (s)')
    ylabel('counts')
    title(sprintf('n = %d, mean = %.2f, median = %.2f, std = %.2f', n, mean(r), median(r), std(r)));
    legend([v1 v2 v3(1)], 'mean', 'median', '1 std from median', ...
        'Location', 'northwest')
    
    % move the title up a little bit
    [ax.Title.Position(1), ax.Title.Position(2)] = ...
        norm2trueposition(ax, 0.5, 1.03);
    
    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_histogram.eps', mfilename);
    figdisp(savename,[],[],2,[],'epstopdf');
end
end