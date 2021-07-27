function params = writestations(params, fname)
% WRITESTATIONS(params, fname)
%
% Writes a list of stations from a STATION file. The receiver sets are 
% named as 'AA', 'AB', ... . In each set, the receivers are named as S0001, 
% S0002, ... . It sets use_existing_STATIONS to true.
%
% INPUT:
% params        parameters
% fname         name of the STATION file
%
% SEE ALSO:
% READ_STATIONS
%
% Last modified by sirawich@princeton.edu, 07/22/2021

fid = fopen(fname, 'w');
for ii = 1:length(params.RECEIVERS)
    r = params.RECEIVERS{ii};
    x = linspace(r.xdeb, r.xfin, r.nrec);
    z = linspace(r.zdeb, r.zfin, r.nrec);
    % array name
    n2 = mod(ii-1, 26);
    n1 = mod(floor((ii-1) / 26), 26);
    network = char([n1+65, n2+65]);
    for jj = 1:r.nrec
        fprintf(fid, ['S%4.4i    %s %20.7f %20.7f       0.0         ' ...
            '0.0\n'], jj, network, x(jj), z(jj));
    end
end
fclose(fid);

% sets use_existing_STATIONS to true
params.use_existing_STATIONS = true;
end