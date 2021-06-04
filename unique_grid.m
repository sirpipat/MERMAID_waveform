function grid = unique_grid(NEX, NEZ, NGLLX, NGLLZ)
% grid = UNIQUE_GRID(NEX, NEZ, NGLLX, NGLLZ)
%
% returns an array of indices for generating a uniquie grid from model
% files for SPECFEM2D.
%
% INPUT:
% NEX           number of elements in X direction
% NEZ           number of elements in Z direction (unused)
% NGLLX         number of GLL points in X direction
% NGLLZ         number of GLL points in Z direction
%
% OUTPUT:
% grid          array of indices for generating a uniquie grid
%
% Last modified by Sirawich Pipatprathanporn, 03/16/2021

grid = zeros(NEZ * (NGLLZ - 1) + 1, NEX * (NGLLX - 1) + 1);

for iex = 1:NEX
    for iez = 1:NEZ
        for igllx = 1:NGLLX
            for igllz = 1:NGLLZ
                if and(or(igllx > 1, iex == 1), or(igllz > 1, iez == 1))
                    ii = (iez - 1) * (NGLLZ - 1) + igllz;
                    jj = (iex - 1) * (NGLLX - 1) + igllx;
                    grid(ii, jj) = iGLOB(iex, iez, igllx, igllz, NEX, NEZ, NGLLX, NGLLZ);
                end
            end
        end
    end
end

end