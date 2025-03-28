function outputdirs = specfem2d_input_setup_response(name, topo, tparams, depth, water, solid, freq, theta, stf_type, Par_file_base, outputdir, saveimage, branch, gpu_mode)
% outputdirs = SPECFEM2D_INPUT_SETUP_RESPONSE(name, topo, tparams, depth, water, solid, freq, theta, stf_type, Par_file_base, outputdir, saveimage, branch, gpu_mode)
% Generates Par_file, source file, and interface file for a fluid-solid
% simulation.
%
% INPUT:
% name              name for the model
% topo              topography end member. Options are the followings:
%                       'flat'          z = B
%                       'incline'       z = A * (x - x0) / x1 + B
%                       'sloped'        z = A/2 / (1 + exp((x - x0) / x1)) + B
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
% depth             depth of the hydrophone
% water             sound speed profile of the water. Options are the followings:
%                       'homogeneous'
%                       'Munk'          Munk sound speed profile
% solid             seismic structure of the solid. Options are the followings:
%                       'homogeneous'   [Default]
%                       'layered'       2 layers with flat surface
% freq              source frequency [Default: 10 Hz]
% theta             incident angle in degrees [Default: 0]
% stf_type          source-time-function type [Default: 1]
%                        1 = second derivative of a Gaussian (a.k.a. Ricker),
%                        2 = first derivative of a Gaussian,
%                        3 = Gaussian,
%                        4 = Dirac,
%                        5 = Heaviside,
%                        6 = ocean acoustics type I,
%                        7 = ocean acoustics type II,
%          'stf_file_name' = external source time function (source read from file),
%                        9 = burst,
%                       10 = Sinus source time function,
%                       11 = Marmousi Ormsby wavelet
% Par_file_base     base Par_file to setting up Par_file
% outputdir         directory for the input files
% saveimage         whether to save the snapshots of not [Default: true]
% branch            SPECFEM2D branch [Default: 'master']
%                   'master' (commit: e937ac2f74f23622f6ebbc8901d30fb33c1a2c38)
%                   'devel'  (commit: cf89366717d9435985ba852ef1d41a10cee97884)
% gpu_mode          whether to enable GPU MODE [Default: false]
%
% OUTPUT:
% outputdirs        output directories for the fluid-solid simulation
%       --- 'master' branch ---
%       outputdirs{1}       pressure hydrophone in the fluid at depth
%       outputdirs{2}       displacement OBS at the ocean floor
%       --- 'devel' branch ---
%       outputdirs          both pressures and displacements at OBS and
%                           hydrophone
%
% SEE ALSO:
% SPECFEM2D_INPUT_SETUP, SPECFEM2D_INPUT_SETUP_FLAT, 
% SPECFEM2D_INPUT_SETUP_TPHASE, RUNFLATSIM
%
% Last modified by sirawich-at-princeton.edu, 10/16/2024

defval('bottom', 4800)
defval('depth', 1500)
defval('water', 'homogeneous')
defval('solid', 'homogeneous')
defval('freq', 10)
defval('theta', 0)
defval('stf_type', 1)
defval('saveimage', true)
defval('branch', 'master')
defval('gpu_mode', false)

% 'master' branch --- indices 1--2
if strcmpi(branch, 'master')
    outputdirs = cell(2,1);
    % run 2 loops: one for 3-component displacement OBS
    %              another for pressure MERMAID
    index_list = 1:2;
% 'devel'  branch --- indices -1
else
    % run 1 loop only for both displacement and pressure
    outputdirs = [];
    index_list = -1;
end

% Handle stf_file input
% If it is a string (file name), it is the type 8 source-time function.
if isstring(stf_type) || ischar(stf_type)
    stf_file = stf_type;
    stf_type = 8;
else
    stf_file = '""';
end

for kk = index_list
    % load baseline Par_file
    if isempty(Par_file_base)
        params = makeparams();
    else
        params = loadparfile(Par_file_base);
    end

    % for 'master' branch
    if kk > 0
        outputdir_kk = sprintf('%s_%d/', outputdir(1:end-1), kk);
        outputdirs{kk, 1} = outputdir_kk;
    % for 'devel' branch
    else
        outputdir_kk = outputdir;
        outputdirs = outputdir;
    end
    % create the DATA/ folder for the files generated by this function
    system(sprintf('mkdir %s', outputdir_kk));
    system(sprintf('mkdir %sDATA/', outputdir_kk));
    
    % save input argument for future reproducibility
    save(sprintf('%sDATA/input_arguments.mat', outputdir_kk), ...
        'name', 'topo', 'tparams', 'depth', 'water', 'solid', 'freq', ...
        'theta', 'Par_file_base', 'outputdir', 'saveimage', 'branch', ...
        'gpu_mode');

    %% define models
    % define material value
    switch lower(solid)
        case 'layered'
            material1 = struct(...
                'model_number'      , 1         , ...
                'type_number'       , 1         , ...
                'type_name'         , 'elastic' , ...
                'rho'               , 3000      , ...
                'vp'                , 6000      , ...
                'vs'                , 3400      , ...
                'QKappa'            , 9999      , ...
                'Qmu'               , 9999        ...
            );
        
            material2 = struct(...
                'model_number'      , 2         , ...
                'type_number'       , 1         , ...
                'type_name'         , 'elastic' , ...
                'rho'               , 2500      , ...
                'vp'                , 3400      , ...
                'vs'                , 1963      , ...
                'QKappa'            , 9999      , ...
                'Qmu'               , 9999        ...
            );


            material3 = struct(...
                'model_number'      , 3         , ...
                'type_number'       , 1         , ...
                'type_name'         , 'acoustic', ...
                'rho'               , 1020      , ...
                'vp'                , 1500      , ...
                'vs'                , 0         , ...
                'QKappa'            , 9999      , ...
                'Qmu'               , 9999        ...
            );

            % set models
            params.nbmodels = 3;
            params.MODELS = {material1, material2, material3};
        otherwise
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

            % set models
            params.nbmodels = 2;
            params.MODELS = {material1, material2};
    end

    %% define geography of the problem

    % set xlimit of the simulation
    xmin = 0;
    xmax = 20000;
    width = xmax - xmin;
    elemsize = 40;
    nx = width / elemsize;

    % set interfaces
    % bottom
    itf1.npts = 2;
    itf1.pts = [xmin 0; xmax 0];
    % ocean bottom
    itf2.npts = 501;
    x = linspace(xmin, xmax, itf2.npts)';
    
    switch lower(topo)
        case 'flat'
            if ~isfield(tparams, 'B') || isempty(tparams.B)
                tparams.B = 4800 * random('unif', 0.9, 1.1);
            end
            z = tparams.B * ones(size(x));
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
                    isempty(tparams.X) || isempty(tparams.Z) || ...
                    all(size(tparams.X) ~= size(tparams.Z))
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
            itf4.npts = 501;
            x = linspace(xmin, xmax, itf4.npts)';
            switch lower(topo)
                case 'sinusoidal'
                    N = unidrnd(10);
                    A = unifrnd(0, 100, [1 N]);
                    B = 2400 * random('unif', 0.9, 1.1);
                    k = unifrnd(1/10000, 1/1000, [1 N]);
                    c = unifrnd(0, 2*pi, [1 N]);
                    zs = B + sum(A .* sin(2 * pi * k .* x + c), 2);
                otherwise
                    zs = 2400 * random('unif', 0.9, 1.1) * ones(size(x));
            end
            itf4.pts = [x,zs];
            itfs = {itf1, itf4, itf2, itf3};
            % set layers [crust1 crust2 ocean]
            nz = zmax / elemsize;
            layers = [round(nz*mean(zs)/zmax) ...
                round(nz*mean(z)/zmax)-round(nz*mean(zs)/zmax) ...
                nz-round(nz*mean(z)/zmax)];
            %layers = [nz/4 nz/4 nz/2];
        otherwise
            itfs = {itf1, itf2, itf3};
            % set layers [crust ocean]
            nz = zmax / elemsize;
            layers = [round(nz*mean(z)/zmax) nz-round(nz*mean(z)/zmax)];
    end

    % write interfaces file
    % set parameters for internal meshing
    params.interfacesfile = sprintf('interfaces_%s.dat', name);
    writeinterfacefile(itfs, layers, sprintf('%sDATA/interfaces_%s.dat', outputdir_kk, name));

    params.read_external_mesh = false;
    params.xmin = xmin;
    params.xmax = xmax;
    params.nx = nx;

    % define regions
    switch lower(solid)
        case 'layered'
            region1 = struct(...
                'nxmin'             , 1         , ...
                'nxmax'             , nx        , ...
                'nzmin'             , 1         , ...
                'nzmax'             , layers(1) , ...
                'material_number'   , 1           ...
            );

            region2 = struct(...
                'nxmin'             , 1         , ...
                'nxmax'             , nx        , ...
                'nzmin'             , layers(1) + 1  , ...
                'nzmax'             , layers(1) + layers(2)      , ...
                'material_number'   , 2           ...
            );
        
            region3 = struct(...
                'nxmin'             , 1         , ...
                'nxmax'             , nx        , ...
                'nzmin'             , layers(1) + layers(2) + 1  , ...
                'nzmax'             , nz        , ...
                'material_number'   , 3           ...
            );
        
            regions = {region1, region2, region3};
            params.nbregions = 3;
        otherwise
            region1 = struct(...
                'nxmin'             , 1         , ...
                'nxmax'             , nx        , ...
                'nzmin'             , 1         , ...
                'nzmax'             , layers(1) , ...
                'material_number'   , 1           ...
            );

            region2 = struct(...
                'nxmin'             , 1         , ...
                'nxmax'             , nx        , ...
                'nzmin'             , layers(1) + 1  , ...
                'nzmax'             , nz        , ...
                'material_number'   , 2           ...
            );
        
            regions = {region1, region2};
            params.nbregions = 2;
    end

    
    params.REGIONS = regions;

    %% define water model
    switch lower(water)
        case 'munk'
            water_model.name = 'Munk';
            water_model.epsilon = 0.00737;
            water_model.zc = 1300;
            water_model.zm = zmax - (min(z) - 500);
            water_model.dz = 10;
            water_model.B = 1300;
        otherwise
            water_model.name = 'homogeneous';
    end

    %% define SOURCE(s)

    % use multiple sources to imitate plane wave
    sources = cell(1, 197);
    for ii = 1:197
        vp = params.MODELS{params.REGIONS{1}.material_number}.vp;
        tshift = (ii - 1) * 100 * sin(theta * pi / 180) / vp;
        source = struct(...
            'source_surf'           , false     , ...   % inside the medium
            'xs'                    , (ii+1) * 100  , ...
            'zs'                    , 720       , ...
            'source_type'           , 2         , ...   % moment tensor
            'time_function_type'    , stf_type  , ...
            'name_of_source_file'   , stf_file  , ...
            'burst_band_width'      , 0         , ...
            'f0'                    , freq      , ...   % dominant frequency
            'tshift'                , tshift    , ...
            'anglesource'           , 0         , ...
            'Mxx'                   , 1.0       , ...   % explosion
            'Mzz'                   , 1.0       , ...   % explosion
            'Mxz'                   , 0.0       , ...   % explosion
            'factor'                , 1e-9 * cos(theta * pi / 180), ...
            'vx'                    , 0.0       , ...
            'vz'                    , 0.0         ...
        );
        sources{ii} = source;
    end
    % add more sources if angle is 5 degrees or above
    if theta >= 5
        tshift = 0;
        pf = depthprofile(200, itfs, 'linear');
        for ii = 1:100
            zs = 720 + ii * 100;
            % for safety purpose as the meshed fluid-solid boundary is not
            % exactly at the interface
            layer = sum(pf - 0 * elemsize < zs);
            % stop if the source is above the bottom layer
            % TODO: figure out how to adjust tshift for sources in
            % different materials
            if layer > 1
                break
            end
            vp = params.MODELS{params.REGIONS{layer}.material_number}.vp;
            tshift = tshift + 100 * cos(theta * pi / 180) / vp;
            source = struct(...
                'source_surf'           , false     , ...   % inside the medium
                'xs'                    , 200       , ...
                'zs'                    , zs        , ...
                'source_type'           , 2         , ...   % moment tensor
                'time_function_type'    , stf_type  , ...
                'name_of_source_file'   , stf_file  , ...
                'burst_band_width'      , 0         , ...
                'f0'                    , freq      , ...   % dominant frequency
                'tshift'                , tshift    , ...
                'anglesource'           , 0         , ...
                'Mxx'                   , 1.0       , ...   % explosion
                'Mzz'                   , 1.0       , ...   % explosion
                'Mxz'                   , 0.0       , ...   % explosion
                'factor'                , 1e-9 * sin(theta * pi / 180), ...
                'vx'                    , 0.0       , ...
                'vz'                    , 0.0         ...
            );
            sources{ii+197} = source;
        end
    end

    params.NSOURCES = length(sources);

    % write source file
    % writesource(sources, sprintf('%sDATA/SOURCE_%s', outputdir_kk, name), branch);
    writesource(sources, sprintf('%sDATA/SOURCE', outputdir_kk), branch);

    %% defeine set(s) of RECEIVER(s)
    bottom = interp1(x, z, 10000, 'linear');
    receiverset1 = struct(...
        'nrec'                              , 1     , ...
        'xdeb'                              , 10000 , ...
        'zdeb'                              , zmax - depth  , ...
        'xfin'                              , 10000 , ...
        'zfin'                              , zmax - depth  , ...
        'record_at_surface_same_vertical'   , false   ...
    );
    receiverset2 = struct(...
        'nrec'                              , 1    , ...
        'xdeb'                              , 10000 , ...
        'zdeb'                              , bottom  , ...
        'xfin'                              , 10000 , ...
        'zfin'                              , bottom  , ...
        'record_at_surface_same_vertical'   , false   ...
    );

    % 'master' branch : hydrophone
    if kk == 1
        params.seismotype = 4;
    % 'master' branch : OBS
    elseif kk == 2
        params.seismotype = 1;
    % 'devel' branch
    else
        params.seismotype = [1 4];
    end
    receiversets = {receiverset1, receiverset2};
    params.nreceiversets = 2;
    params.RECEIVERS = receiversets;

    % writes STATIONS file to include receiverset information when there are
    % multiple receiversets.
    params = writestations(params, sprintf('%sDATA/STATIONS', outputdir_kk));
    %% define other parameters
    params.title = sprintf('fluid/solid interface : %s -- on %s branch', ...
        name, upper(branch));
    params.time_stepping_scheme = 1;
    params.NSTEP = 65000;
    params.DT = 5e-4;
    params.subsamp_seismos = round(1 /params.DT / 100); % reduce to 100 Hz
    params.save_ASCII_seismograms = 0;
    params.save_binary_seismograms_single = 1;
    params.save_binary_seismograms_double = 0;
    params.SU_FORMAT = 0;
    params.NTSTEP_BETWEEN_OUTPUT_ENERGY = 1000;
    params.PML_BOUNDARY_CONDITIONS = false;
    params.STACEY_ABSORBING_CONDITIONS = true;
    params.GPU_MODE = gpu_mode;
    params.factor_subsample_image = 5.0;
    params.MODEL = 'default';
    % do not save the images when SAVEIMAGE is set to false
    params.output_color_image = saveimage;
    params.output_postscript_snapshot = false; %saveimage;
    %% write Par_file
    % writeparfile(params, sprintf('%sDATA/Par_file_%s', outputdir_kk, name), branch);
    writeparfile(params, sprintf('%sDATA/Par_file', outputdir_kk), branch);
    %% write a supplementary file for runthisexample.m
    %  It is not used by specfem2d.
    save(sprintf('%sDATA/supplementary_%s.mat', outputdir_kk, name), 'water_model'); 
end
end