function [fc, s] = freqselect_routine(sacfiles, plt, option)
% [fc, s] = FREQSELECT_ROUTINE(sacfiles, plt, option)
%
% Figures out the frequency band where the signal stands out the most from
% the background noise in multiple seismograms.
% 
% INPUT:
% sacfiles      cell array to sacfiles
% plt           whether to plot and save figure or not
% option        how to select the best corner frequency [default: 4]
%               1 -- highest bandpass SNR
%               2 -- highest bandpass SNR to bandstop SNR ratio
%               3 -- widest bandwidth that keep bandpass SNR > 50% of the 
%               highest
%               4 -- widest bandwidth that keep SNR ratio > 50% of the
%               highest
%               5 -- argmax(bandpass SNR + 1/(1 - bandstop SNR))
%
% OUTPUT:
% fc            best corner frequency for each seismogram
% s             best signal-to-noise ratio for each seismogram
%
% Last modified by sirawich-at-princeton.edu, 02/12/2024

defval('option', 4)

badval = -12345;

fc = zeros(length(sacfiles), 2);
s = zeros(length(sacfiles), 3);
for ii = 1:length(sacfiles)
    [seisdata, hdrdata] = readsac(sacfiles{ii});
    [dt_ref, dt_B, ~, fs, ~, dts] = gethdrinfo(hdrdata);
    
    % determine arrival time if it does not exist
    if isnan(hdrdata.T0) || hdrdata.T0 == badval
        % first, try to get travel time if earthquake info is available
        if ~isnan(hdrdata.USER7) && hdrdata.USER7 ~= badval
            try
                ev = irisFetch.Events('eventID', string(hdrdata.USER7));
                tt = tauptime('mod', 'ak135', 'dep', hdrdata.EVDP, ...
                    'ph', 'p,P,Pdiff,PKP,PKIKP', 'deg', hdrdata.GCARC);
                hdrdata.T0 = seconds(datetime(ev(1).PreferredTime, 'Format', ...
                    'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', 'TimeZone', 'UTC') + ...
                    seconds(tt(1).time) - dt_ref);
            catch ME
                fprintf('%s\n', ME.getReport);
                keyboard
            end
        % else, pick the arrival time 
        else
            % TODO: Implement
            continue
        end
        
        if isnan(hdrdata.T0) || hdrdata.T0 == badval
            keyboard
        end
    end
    
    t = seconds(dts - dt_ref) - hdrdata.T0;
    pa = real(counts2pa(seisdata, fs, [0.01 0.02 10 20], [], 'sacpz', false));
    titlename = sprintf('Event ID: %d, Magnitude: %.2f, Distance: %.2f^{\\circ}, Station: %s', ...
        hdrdata.USER7, hdrdata.MAG, hdrdata.GCARC, strtrim(hdrdata.KSTNM));
    savename = sprintf('%d_%s', hdrdata.USER7, strtrim(hdrdata.KSTNM));
    try
        [fc(ii,:), s(ii,:)] = freqselect(t, pa, fs, plt, titlename, ...
            savename, option);
    catch ME
        fprintf('%s\n', ME.getReport);
        continue
    end
end
end