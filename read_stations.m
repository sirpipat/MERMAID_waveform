function [n, name, network, x, z] = read_stations(fname)
% [n, name, network, x, z] = read_stations(fname)
% 
% Reads a STATION file for SPECFEM2D simulation.
%
% INPUT
% fname         full filename of a STATION file
%
% OUTPUT
% n             the number of stations
% name          station names
% network       station network names
% x             x-coordinates
% z             z-coordinates
%
% Last modified by Sirawich Pipatprathanporn, 11/09/2021

% read the station file as a table
opts = detectImportOptions(fname, 'FileType', 'text');
T = readtable(fname, opts);

% if the table is empty, assume the table has only one entry
if isempty(T)
    fid = fopen(fname, 'r');
    line = fgetl(fid);
    words = split(line);
    T = struct('Var1', {words(1)}, 'Var2', {words(2)}, ...
        'Var3', str2double(words{3}), ...
        'Var4', str2double(words{4}));
end

name = T.Var1;
network = T.Var2;
x = T.Var3;
z = T.Var4;
n = size(T, 1);
end