function [ft, tr, hdr] = extractfeatures(allsacfiles, arr, feat_select)
% Extracts features from seismograms, which includes
% 1. amplitude of the first peak
% 2. duration between arrial time and the first peak
% 3. amplitude of the noise

defval('plt', false)
defval('arr', 1)

numfiles = length(allsacfiles);
numfeatures = length(feat_select);

ft = zeros(numfiles, numfeatures+2);
tr = cell(numfiles, 1);
hdr = cell(numfiles, 1);

for ii = 1:numfiles
    [SeisData, HdrData, ~, ~, tims] = readsac(allsacfiles{ii});
    fs = 1/(tims(2) - tims(1));
    
    % gather dist and azimuth data
    ft(ii, 1) = HdrData.GCARC;
    ft(ii, 2) = HdrData.AZ;
    
    % convert digital counts to pressure
    x = real(counts2pa(SeisData, fs));
    % filter out the ambient noise below 1 Hz
    x = bandpass(x, fs, 1, 2, 2, 1, 'butter', 'linear');
    
    tr{ii, 1} = x;
    hdr{ii, 1} = HdrData;
    
    % split the section into before and after arrival
    if arr == 1
        wh = (tims < HdrData.T0);
    else
        wh = (tims < HdrData.T0 - HdrData.USER4);
    end
    tims_before = tims(wh);
    x_before = x(wh);
    tims_after = tims(~wh);
    x_after = x(~wh);
    [pks_after, locs_after] = findpeaks(abs(x_after));
    [pks_before, locs_before] = findpeaks(abs(x_before));
    pks_before = pks_before(4:end) .* sign(x_before(locs_before(4:end)));
    for jj = 1:numfeatures
        switch feat_select(jj)
            case 1
                ft(ii, 2+feat_select(jj)) = pks_after(1) * sign(x_after(locs_after(1)));
            case 2
                ft(ii, 2+feat_select(jj)) = tims_after(locs_after(1)) - HdrData.T0;
            case 3
                ft(ii, 2+feat_select(jj)) = indeks(maxk(abs(pks_before), 5), 1);
            otherwise
        end
    end
end
end