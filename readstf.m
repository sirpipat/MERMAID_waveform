function stf = readstf(fname)
% stf = READSTF(fname)
%
% Reads a source-time function file with the format described in SCARDEC
% website. See http://scardec.projects.sismo.ipgp.fr
%
% INPUT:
% fname         name of the source-time function file
% 
% OUTPUT:
% stf           a struct represnting a source-time function with metadata
%
% Last modified by sirawich-at-princeton.edu, 07/14/2022

fid = fopen(fname, 'r');

stf.year     = fscanf(fid, '%d', 1);
stf.month    = fscanf(fid, '%d', 1);
stf.day      = fscanf(fid, '%d', 1);
stf.hour     = fscanf(fid, '%d', 1);
stf.minute   = fscanf(fid, '%d', 1);
stf.second   = fscanf(fid, '%f', 1);
stf.lat      = fscanf(fid, '%f', 1);
stf.lon      = fscanf(fid, '%f', 1);

stf.depth_km = fscanf(fid, '%f', 1);
stf.M0       = fscanf(fid, '%g', 1);
stf.Mw       = fscanf(fid, '%f', 1);
stf.strike1  = fscanf(fid, '%f', 1);
stf.dip1     = fscanf(fid, '%f', 1);
stf.rake1    = fscanf(fid, '%f', 1);
stf.strike2  = fscanf(fid, '%f', 1);
stf.dip2     = fscanf(fid, '%f', 1);
stf.rake2    = fscanf(fid, '%f', 1);

data         = fscanf(fid, '%f', [2 Inf])';
fclose(fid);

stf.t           = data(:,1);
stf.moment_rate = data(:,2);
stf.dt          = data(2,1) - data(1,1);
end