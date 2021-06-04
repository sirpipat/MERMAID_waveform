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
% Last modified by Sirawich Pipatprathanporn, 02/24/2021

% read the station file as a table
opts = detectImportOptions(fname);
T = readtable(fname, opts);

name = T.Var1;
network = T.Var2;
x = T.Var3;
z = T.Var4;
n = size(T, 1);
end