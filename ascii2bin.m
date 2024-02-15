function ascii2bin(ddir, precision)
% ASCII2BIN(ddir, precision)
%
% Converts SPECFEM2D output seismograms from ASCII format to binary at
% either single or double precision.
%
% INPUT:
% ddir          directory of a SPECFEM2D simulation output containing
%               OUTPUT_FILES/
% precision     floating-point precsion either
%               'single' [default] or 'double'
%
% OUTPUT:
% No output. Only binary files written at OUTPUT_FILES/
%
% Last modified by sirawich-at-princeton.edu, 02/15/2024

% seismotype 1-3: displacement, velocity, acceleration
seismotypes = {'d', 'v', 'a'};
channels = {'BXX', 'BXY', 'BXZ'};
channels_bin = {'x', 'y', 'z'};
writebinaryfile(ddir, precision, seismotypes, channels, channels_bin);

% seismotype 4: pressure
seismotypes = {'p'};
channels = {'PRE'};
channels_bin = {'p'};
writebinaryfile(ddir, precision, seismotypes, channels, channels_bin);

% seismotype 5: curl
seismotypes = {'c'};
channels = {'BXX', 'BXY', 'BXZ', 'cur'};
channels_bin = {'x', 'y', 'z', 'c'};
writebinaryfile(ddir, precision, seismotypes, channels, channels_bin);

% seismotype 6: fluid potential
seismotypes = {'x'};
channels = {'POT'};
channels_bin = {'p'};
writebinaryfile(ddir, precision, seismotypes, channels, channels_bin);
end

% Writes binary files from ASCII seismogram files.
%
% INPUT:
% ddir              directory of a SPECFEM2D simulation output containing
%                   OUTPUT_FILES/
% precision         floating-point precsion either
%                   'single' [default] or 'double'
% seismotypes       list of seimogram types e.g. {'d', 'v', 'a'} and {'p'}
% channels          list of channel names for ASCII seimogram file names
%                   e.g. {'BXX', 'BXY', 'BXZ'} and {'PRE'}
% channels_bin      list of channel names for binary file names
%                   e.g. {'x', 'y', 'z'} and {'p'}
%
% channels and channels_bin should have the same number of elements and
% match 1-by-1.
function writebinaryfile(ddir, precision, seismotypes, channels, ...
    channels_bin)
for ii = 1:length(seismotypes)
    for jj = 1:length(channels)
        try
            ascii_files = ls2cell(sprintf('%sOUTPUT_FILES/*%s.sem%s', ...
                ddir, channels{jj}, seismotypes{ii}), 1);
        catch
            continue
        end
        fid = fopen(sprintf('%sOUTPUT_FILES/U%s_file_%s_%s.bin', ddir, ...
            channels_bin{jj}, precision, seismotypes{ii}), 'wb');
        for kk = 1:length(ascii_files)
            [~, x] = read_seismogram(ascii_files{kk});
            fwrite(fid, x, precision);
        end
        fclose(fid);
    end
end
end