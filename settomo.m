function grid = settomo(data, grid)
% grid = SETTOMO(data, grid)
%
% SET 
% data.z
% data.vp
% data.rho
% grid.x
% grid.z
% grid.vp
% grid.vs
% grid.rho

% create an interpolator
% vp(z), rho(z)
% TODO: make them dependent on both x and z
vp_interp = griddedInterpolant(data.z, data.vp, 'makima');
rho_interp = griddedInterpolant(data.z, data.rho, 'makima');

% set the grid
for iGLOB = 1:length(data.x)
    grid.vp(iGLOB) = vp_interp(grid.z(iGLOB));
    grid.rho(iGLOB) = rho_interp(grid.z(iGLOB));
end
end