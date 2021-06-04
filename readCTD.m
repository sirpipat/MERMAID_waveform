function ctd_profile = readCTD(ctdfile)
z = ncread(ctdfile, 'z');
T = ncread(ctdfile, 'Temperature');
S = ncread(ctdfile, 'Salinity');

ctd_profile = struct('z', z, 'T', T, 'S', S);
end