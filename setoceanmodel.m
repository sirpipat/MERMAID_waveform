function grid = setoceanmodel(data, grid, N, Ninterface)
% grid = SETOCEANMODEL(data, grid, N, Ninterface)
%
% SET 
% data.d
% data.vp
% data.rho
% grid.x
% grid.z
% grid.vp
% grid.vs
% grid.rho
% N.NEX
% N.NEZ
% N.NGLLX
% N.NGLLZ

% create an interpolator
% vp(z), rho(z)
z_water = max(grid.z, [], 'all') - data.d;
% make z_water monotonically increasing
[z_water, iz] = sort(z_water, 'ascend');
data.vp = data.vp(iz);
data.rho = data.rho(iz);
vp_interp = griddedInterpolant(z_water, data.vp, 'makima');
rho_interp = griddedInterpolant(z_water, data.rho, 'makima');

% set the grid
for iglob = 1:length(grid.x)
    [~, iez, ~, ~] = unpack_iGLOB(iglob, N.NEX, N.NEZ, N.NGLLX, N.NGLLZ);
    if iez > Ninterface
        grid.vp(iglob) = vp_interp(grid.z(iglob));
        grid.rho(iglob) = rho_interp(grid.z(iglob));
    end
end
end