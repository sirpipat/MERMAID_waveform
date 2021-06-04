function [iex, iez, igllx, igllz] = unpack_iGLOB(iglob, NEX, NEZ, NGLLX, NGLLZ)
% [iex, iez, igllx, igllz] = UNPACK_IGLOB(iglob, NEX, NEZ, NGLLX, NGLLZ)
%
% calculates element indices and GLL indices from the global index.
%
% INPUT:
% iglob         global index
% NEX           number of elements in X direction
% NEZ           number of elements in Z direction (unused)
% NGLLX         number of GLL points in X direction
% NGLLZ         number of GLL points in Z direction
%
% OUTPUT:
% iex           X element index     [1, NEX]
% iez           Z element index     [1, NEZ]
% igllx         GLLX index          [1, NGLLX]
% igllz         GLLZ index          [1, NGLLZ]
%
% Last modified by Sirawich Pipatprathanporn, 03/15/2021

igllx = zeros(size(iglob));
igllz = zeros(size(iglob));
iex = zeros(size(iglob));
iez = zeros(size(iglob));
for ii = 1:length(iglob)
    jj = iglob(ii);
    if or(jj < 2, jj > NEX * NEZ * NGLLX * NGLLZ + 1)
        igllx(ii) = nan;
        igllz(ii) = nan;
        iex(ii) = nan;
        iez(ii) = nan;
    else
        igllx(ii) = 1 + mod(jj - 2, NGLLX);
        igllz(ii) = 1 + mod(floor((jj - 2) / NGLLX), NGLLZ);
        iex(ii) = 1 + mod(floor((jj - 2) / (NGLLX * NGLLZ)), NEX);
        iez(ii) = 1 + floor((jj - 2) / (NGLLX * NGLLZ * NEX));
    end
end
end