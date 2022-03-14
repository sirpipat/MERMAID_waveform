function specfem2d_input_setup(name, topo, tparams, water, solid, source, freq, angle, Par_file_base, outputdir, branch)
% SPECFEM2D_INPUT_SETUP(name, topo, tparams, water, solid, source, freq, angle, Par_file_base, outputdir, branch)
%
% Generates Par_file, source file, and interface file for a fluid-solid
% simulation.
%
% INPUT:
% name              name for the model
% topo              topography end member. Options are the followings:
%                       'flat'          z = B
%                       'incline'       z = A * (x - x0) / x1 + B
%                       'sloped'        z = A/2 / (1 + exp((x - x0) / x1)) + B  |||| A tanh(-(x - x0) / x1) + B
%                       'sinusoidal'    z = \sum (A sin(2*pi*k*x + c)) + B
%                       'parabolic'     z = A (x - x0)^2 + B
%                       'hill'          z = \sum (A exp(-(x - x0)^2 / x1^2)) + B
%                       'custom'        Z = f(X)
% tparams           struct containing parameters to define the topography
%                       N  -- the number of such feature
%                       A  -- amplitude(s) of the feature
%                       B  -- flat elevation or offset
%                       x0 -- central location(s) of the feature
%                       x1 -- length scale(s)
%                       k  -- wave number(s)
%                       c  -- phase(s)
%                       X  -- x-coordinate for 'custom' topography
%                       Z  -- z-coordinate for 'custom' topography
%                   Refer to the equations above for the specific usages.
%                   If the parameters are not defined, they will be assined
%                   to random values.
% water             sound speed profile of the water. Options are the followings:
%                       'homogenous'
%                       'Munk'          Munk sound speed profile
% solid             seismic structure of the solid. Options are the followings:
%                       'homogeneous'
%                       'layered'       2 layers with flat surface
% source            proximity of the source
%                       'shallow'
%                       'deep' -or- 'distant'
% freq              source frequency [Default: 10 Hz]
% angle             incident angle in degrees [Default: 0]
% Par_file_base     base Par_file to setting up Par_file
% outputdir         directory for the input files
% branch            SPECFEM2D branch [Default: 'master']
%                   'master' (commit: e937ac2f74f23622f6ebbc8901d30fb33c1a2c38)
%                   'devel'  (commit: cf89366717d9435985ba852ef1d41a10cee97884)
%
% Last modified by sirawich-at-princeton.edu, 03/14/2022

defval('topo', 'flat')
defval('water', 'homogeneous')
defval('solid', 'homogeneous')
defval('source', 'shallow')
defval('freq', 10)
defval('angle', 0)
defval('branch', 'master')

% load baseline Par_file
if isempty(Par_file_base)
    params = makeparams();
else
    params = loadparfile(Par_file_base);
end

% create the DATA/ folder for the files generated by this function
system(sprintf('mkdir %s', outputdir));
system(sprintf('mkdir %sDATA/', outputdir));

%% define models
% define material value
material1 = struct(...
    'model_number'      , 1         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'elastic' , ...
    'rho'               , 2500      , ...
    'vp'                , 3400      , ...
    'vs'                , 1963      , ...
    'QKappa'            , 9999      , ...
    'Qmu'               , 9999        ...
);

material2 = struct(...
    'model_number'      , 2         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'acoustic', ...
    'rho'               , 1020      , ...
    'vp'                , 1500      , ...
    'vs'                , 0         , ...
    'QKappa'            , 9999      , ...
    'Qmu'               , 9999        ...
);

material3 = struct(...
    'model_number'      , 3         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'elastic' , ...
    'rho'               , 3300      , ...
    'vp'                , 4400      , ...
    'vs'                , 2540      , ...
    'QKappa'            , 9999      , ...
    'Qmu'               , 9999        ...
);

% set models
switch lower(solid)
    % TODO: add the second solid layer and
    case 'layered'
        params.nbmodels = 3;
        params.MODELS = {material1, material2, material3};
    otherwise
        params.nbmodels = 2;
        params.MODELS = {material1, material2};
end

%% define geography of the problem

% set xlimit of the simulation
xmin = 0;
xmax = 20000;
width = xmax - xmin;
elemsize = 400/9; % mean element size is 44.44 m (for T_phase, it is 100 m)
nx = width / elemsize;

% set interfaces
% bottom
itf1.npts = 2;
itf1.pts = [xmin 0; xmax 0];
% ocean bottom
itf2.npts = 401;
x = linspace(xmin, xmax, itf2.npts)';
switch lower(topo)
    case 'incline'
        if ~isfield(tparams, 'A') || isempty(tparams.A)
            tparams.A = 1200 * random('unif', 0.9, 1.1);
        end
        if ~isfield(tparams, 'B') || isempty(tparams.B)
            tparams.B = 4800 * random('unif', 0.9, 1.1);
        end
        if ~isfield(tparams, 'x0') || isempty(tparams.x0)
            tparams.x0 = (xmin + xmax) / 2;
        end
        if ~isfield(tparams, 'x1') || isempty(tparams.x1)
            tparams.x1 = (xmax - xmin) / 2;
        end
        z = tparams.B + tparams.A * (x - tparams.x0) / tparams.x1;
    case 'sloped'
        if ~isfield(tparams, 'A') || isempty(tparams.A)
            tparams.A = 4200 * random('unif', 0.99, 1.01);
        end
        if ~isfield(tparams, 'B') || isempty(tparams.B)
            tparams.B = 4800 * random('unif', 0.99, 1.01);
        end
        if ~isfield(tparams, 'x0') || isempty(tparams.x0)
            tparams.x0 = 5400 * random('unif', 0.9, 1.1);
        end
        if ~isfield(tparams, 'x1') || isempty(tparams.x1)
            tparams.x1 = 1000 * random('unif', 0.9, 1.1);
        end
        z = tparams.B + tparams.A/2 / ...
            (1 + exp((x - tparams.x0) / tparams.x1));
        %z = B + A/2 * (1 + tanh(- (x - x0) / x1));
    case 'sinusoidal'
        if ~isfield(tparams, 'N') || isempty(tparams.N)
            tparams.N = unidrnd(3);
        end
        if ~isfield(tparams, 'A') || isempty(tparams.A)
            tparams.A = unifrnd(0, 400, [1 tparams.N]);
        end
        if ~isfield(tparams, 'B') || isempty(tparams.B)
            tparams.B = 4800 * random('unif', 0.9, 1.1);
        end
        if ~isfield(tparams, 'k') || isempty(tparams.k)
            tparams.k = unifrnd(1/10000, 1/1000, [1 tparams.N]);
        end
        if ~isfield(tparams, 'c') || isempty(tparams.c)
            tparams.c = unifrnd(0, 2*pi, [1 tparams.N]);
        end
        z = tparams.B + sum(tparams.A .* sin(2 * pi * ...
            tparams.k .* x + tparams.c), 2);
    case 'parabolic'
        if ~isfield(tparams, 'x0') || isempty(tparams.x0)
            tparams.x0 = (xmin + xmax) / 2;
        end
        if ~isfield(tparams, 'A') || isempty(tparams.A)
            tparams.A = 4000 * (1/(xmax - tparams.x0)^2) * random('unif', 0.2, 1);
        end
        if ~isfield(tparams, 'B') || isempty(tparams.B)
            tparams.B = 4800 * random('unif', 0.9, 1.1);
        end
        z = tparams.A * (x - tparams.x0) .^ 2 + tparams.B;
    case 'hill'
        if ~isfield(tparams, 'N') || isempty(tparams.N)
            tparams.N = unidrnd(3);
        end
        if ~isfield(tparams, 'A') || isempty(tparams.A)
            tparams.A = unifrnd(200, 1000, [1 tparams.N]);
        end
        if ~isfield(tparams, 'B') || isempty(tparams.B)
            tparams.B = 4800 * random('unif', 0.9, 1.1);
        end
        if ~isfield(tparams, 'x0') || isempty(tparams.x0)
            tparams.x0 = unifrnd(xmin, xmax, [1 tparams.N]);
        end
        if ~isfield(tparams, 'x1') || isempty(tparams.x1)
            tparams.x1 = unifrnd(400, width/10, [1 tparams.N]);
        end
        z = tparams.B + sum(tparams.A .* exp(- (x - tparams.x0) .^ 2 ./ ...
            (tparams.x1 .^ 2)), 2);
    case 'custom'
        if ~isfield(tparams, 'X') || ~isfield(tparams, 'Z') || ...
                ~isempty(tparams.X) || ~isempty(tparams.Z) || ...
                size(tparams.X) ~= size(tparams.Z)
            z = 4800 * random('unif', 0.9, 1.1) * ones(size(x));
        else
            z = interp1(tparams.X, tparams.Z, x, 'linear', 'extrap');
        end
    otherwise
        z = 4800 * random('unif', 0.9, 1.1) * ones(size(x));
end
itf2.pts = [x , z];
% sea surface
zmax = 9600;
itf3.npts = 2;
itf3.pts = [xmin zmax; xmax zmax];



switch lower(solid)
    case 'layered'
        % solid-solid interface for layered crust model
        itf4.npts = 401;
        x = linspace(xmin, xmax, itf4.npts)';
        switch lower(topo)
            case 'sinusoidal'
                N = unidrnd(10);
                A = unifrnd(0, 100, [1 N]);
                B = 2400 * random('unif', 0.9, 1.1);
                k = unifrnd(1/10000, 1/1000, [1 N]);
                c = unifrnd(0, 2*pi, [1 N]);
                z = B + sum(A .* sin(2 * pi * k .* x + c), 2);
            otherwise
                z = 2400 * random('unif', 0.9, 1.1) * ones(size(x));
        end
        itf4.pts = [x,z];
        itfs = {itf1, itf4, itf2, itf3};
        % set layers [crust1 crust2 ocean]
        nz = zmax / elemsize;
        layers = [nz/4 nz/4 nz/2];
    otherwise
        itfs = {itf1, itf2, itf3};
        % set layers [crust ocean]
        nz = zmax / elemsize;
        layers = [nz/2 nz/2];
end

% write interfaces file
% set parameters for internal meshing
params.interfacesfile = sprintf('interfaces_%s.dat', name);
writeinterfacefile(itfs, layers, sprintf('%sDATA/interfaces_%s.dat', outputdir, name));

params.read_external_mesh = false;
params.xmin = xmin;
params.xmax = xmax;
params.nx = nx;

% define regions
region1 = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , nx        , ...
    'nzmin'             , 1         , ...
    'nzmax'             , nz/2      , ...
    'material_number'   , 1           ...
);

region2 = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , nx        , ...
    'nzmin'             , nz/2 + 1  , ...
    'nzmax'             , nz        , ...
    'material_number'   , 2           ...
);

region1a = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , nx        , ...
    'nzmin'             , 1         , ...
    'nzmax'             , nz/4      , ...
    'material_number'   , 3           ...
);

region1b = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , nx        , ...
    'nzmin'             , nz/4 + 1  , ...
    'nzmax'             , nz/2      , ...
    'material_number'   , 1           ...
);


switch lower(solid)
    case 'layered'
        regions = {region1a, region1b, region2};
        params.nbregions = 3;
        params.REGIONS = regions;
    otherwise
        regions = {region1, region2};
        params.nbregions = 2;
        params.REGIONS = regions;
end

%% define water model
switch lower(water)
    case 'munk'
        water_model.name = 'Munk';
        water_model.epsilon = 0.00737 * random('unif', 0.9, 1.1);
        water_model.zc = 1300 * random('unif', 0.5, 2);
        water_model.zm = zmax - (min(z) - 500);
        water_model.dz = 10;
        water_model.B = 1300 * random('unif', 0.5, 2);
    otherwise
        water_model.name = 'homogeneous';
end

%% define SOURCE(s)
% local, point-like source
switch lower(source)
    case 'shallow'
        source = struct(...
            'source_surf'           , false     , ...   % inside the medium
            'xs'                    , 1000      , ...
            'zs'                    , 720       , ...
            'source_type'           , 2         , ...   % moment tensor
            'time_function_type'    , 1         , ...   % Ricker
            'name_of_source_file'   , '""'      , ...   % blank for now
            'burst_band_width'      , 0         , ...
            'f0'                    , freq      , ...   % dominant frequency
            'tshift'                , 0         , ...
            'anglesource'           , 0         , ...
            'Mxx'                   , 1.0       , ...   % explosion
            'Mzz'                   , 1.0       , ...   % explosion
            'Mxz'                   , 0.0       , ...   % explosion
            'factor'                , 1e-9      , ...
            'vx'                    , 0.0       , ...
            'vz'                    , 0.0         ...
            );
        sources{1} = source;
    otherwise
        % use multiple sources to imitate plane wave
        sources = cell(1, 197);
        for ii = 1:197
            vp = params.MODELS{params.REGIONS{1}.material_number}.vp;
            tshift = (ii - 1) * 100 * sin(angle * pi / 180) / vp;
            source = struct(...
                'source_surf'           , false     , ...   % inside the medium
                'xs'                    , (ii+1) * 100  , ...
                'zs'                    , 720       , ...
                'source_type'           , 2         , ...   % moment tensor
                'time_function_type'    , 1         , ...   % Ricker
                'name_of_source_file'   , '""'      , ...   % blank for now
                'burst_band_width'      , 0         , ...
                'f0'                    , freq      , ...   % dominant frequency
                'tshift'                , tshift    , ...
                'anglesource'           , 0         , ...
                'Mxx'                   , 1.0       , ...   % explosion
                'Mzz'                   , 1.0       , ...   % explosion
                'Mxz'                   , 0.0       , ...   % explosion
                'factor'                , 1e-9 * cos(angle * pi / 180), ...
                'vx'                    , 0.0       , ...
                'vz'                    , 0.0         ...
            );
            sources{ii} = source;
        end
        % add more sources if angle is 5 degrees or above
        if angle >= 5
            tshift = 0;
            pf = depthprofile(200, itfs, 'linear');
            for ii = 1:50
                zs = 720 + ii * 100;
                layer = sum(pf < zs);
                % stop if the source is above the bottom layer
                % TODO: figure out how to adjust tshift for sources in
                % different materials
                if layer > 1
                    break
                end
                vp = params.MODELS{params.REGIONS{layer}.material_number}.vp;
                tshift = tshift + 100 * cos(angle * pi / 180) / vp;
                source = struct(...
                    'source_surf'           , false     , ...   % inside the medium
                    'xs'                    , 200       , ...
                    'zs'                    , zs        , ...
                    'source_type'           , 2         , ...   % moment tensor
                    'time_function_type'    , 1         , ...   % Ricker
                    'name_of_source_file'   , '""'      , ...   % blank for now
                    'burst_band_width'      , 0         , ...
                    'f0'                    , freq      , ...   % dominant frequency
                    'tshift'                , tshift    , ...
                    'anglesource'           , 0         , ...
                    'Mxx'                   , 1.0       , ...   % explosion
                    'Mzz'                   , 1.0       , ...   % explosion
                    'Mxz'                   , 0.0       , ...   % explosion
                    'factor'                , 1e-9 * sin(angle * pi / 180), ...
                    'vx'                    , 0.0       , ...
                    'vz'                    , 0.0         ...
                );
                sources{ii+197} = source;
            end
        end
end


%sources = {source1, source2, source3, source4, source5, source6, source7};
    
params.NSOURCES = length(sources);

% write source file
writesource(sources, sprintf('%sDATA/SOURCE_%s', outputdir, name), branch);

%% defeine set(s) of RECEIVER(s)
receiverset1 = struct(...
    'nrec'                              , 29    , ...
    'xdeb'                              , 10000 , ...
    'zdeb'                              , 4950  , ...
    'xfin'                              , 10000 , ...
    'zfin'                              , 9450  , ...
    'record_at_surface_same_vertical'   , false   ...
);

receiverset2 = struct(...
    'nrec'                              , 39    , ...
    'xdeb'                              , 500  , ...
    'zdeb'                              , 8100  , ...
    'xfin'                              , 19500 , ...
    'zfin'                              , 8100  , ...
    'record_at_surface_same_vertical'   , false   ...
);

% receiverset2 = struct(...
%     'nrec'                              , 56    , ...
%     'xdeb'                              , 16500  , ...
%     'zdeb'                              , 8100  , ...
%     'xfin'                              , 99000 , ...
%     'zfin'                              , 8100  , ...
%     'record_at_surface_same_vertical'   , false   ...
% );
% 
% receiverset3 = struct(...
%     'nrec'                              , 29    , ...
%     'xdeb'                              , 50000 , ...
%     'zdeb'                              , 4950  , ...
%     'xfin'                              , 50000 , ...
%     'zfin'                              , 9450  , ...
%     'record_at_surface_same_vertical'   , false   ...
% );
% 
% receiverset4 = struct(...
%     'nrec'                              , 29    , ...
%     'xdeb'                              , 95000 , ...
%     'zdeb'                              , 4950  , ...
%     'xfin'                              , 95000 , ...
%     'zfin'                              , 9450  , ...
%     'record_at_surface_same_vertical'   , false   ...
% );

receiversets = {receiverset2, receiverset1};

params.seismotype = 2;
params.nreceiversets = 2;
params.RECEIVERS = receiversets;

% writes STATIONS file to include receiverset information when there are
% multiple receiversets.
params = writestations(params, sprintf('%sDATA/STATIONS', outputdir));
%% define other parameters
params.title = sprintf('fluid/solid interface : %s', name);
params.time_stepping_scheme = 1;
params.NSTEP = 25000;   % T_phase 100000
params.DT = 5e-4;       % T_phase 1e-3
%% write Par_file
writeparfile(params, sprintf('%sDATA/Par_file_%s', outputdir, name), branch);
%% write a supplementary file for runthisexample.m
%  It is not used by specfem2d.
save(sprintf('%sDATA/supplementary_%s.mat', outputdir, name), 'water_model');
end

% ocean bottom is modeled as a modified logistic function
% add more parameters for 4800, 4000 etc. or even random parameters
function z = func(x)
z = 4800 + 4000 * (1 ./ (1 + exp((x-5400)/1000)));
end