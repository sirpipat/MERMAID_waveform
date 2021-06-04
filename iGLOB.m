function iglob = iGLOB(iex, iez, igllx, igllz, NEX, NEZ, NGLLX, NGLLZ)
% iglob = IGLOB(iex, iez, igllx, igllz, NEX, NEX, NGLLX, NGLLZ)
%
% returns a global index for any given element indices and the numbers of
% GLL points and elements in X and Z directions.
%
% INPUT:
% iex           X element index     [1, NEX]
% iez           Z element index     [1, NEZ]
% igllx         GLLX index          [1, NGLLX]
% igllz         GLLZ index          [1, NGLLZ]
% NEX           number of elements in X direction
% NEZ           number of elements in Z direction (unused)
% NGLLX         number of GLL points in X direction
% NGLLZ         number of GLL points in Z direction
%
% OUTPUT:
% iglob         global index
%
% Last modified by Sirawich Pipatprathanporn, 03/15/2021

iglob = 2 + (igllx - 1) + NGLLX * (igllz - 1 + NGLLZ  * (iex - 1 + NEX * (iez - 1)));
end