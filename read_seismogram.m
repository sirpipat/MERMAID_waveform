function [t, x] = read_seismogram(filename)
% [t, x] = READ_SEISMOGRAM(filename)
%
% Read a seismogram output (text file) from SPECFEM
%
% INPUT
% filename      name of the file
%
% OUTPUT
% t             time
% x             data
%
% Last modified by Sirawich Pipatprathanporn, 05/31/21

% read seismograms
sizeData = [2 Inf];
fid = fopen(filename,'r');
data = fscanf(fid, '%f %f', sizeData);
fclose(fid);

t = data(1,:);
x = data(2,:);
end