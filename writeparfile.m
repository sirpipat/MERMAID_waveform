function writeparfile(params, fname, branch)
% WRITEPARFILE(params, fname, branch)
%
% Writes parameters from a Par_file.
%
% DISCLAIMER: This is not the official way to read/write Par_file. I just
% go through comments and parameters in an instant of Par_file and
% read/write accordingly.
%
% INPUT:
% params        parameters
% fname         name of the Par_file
% branch        SPECFEM2D branch [Default: 'master']
%               'master' (commit: e937ac2f74f23622f6ebbc8901d30fb33c1a2c38)
%               'devel'  (commit: cf89366717d9435985ba852ef1d41a10cee97884)
%               'pp-devel' (commit: TBD) -- calls WRITEPARFILEPP
%
% SEE ALSO:
% LOADPARFILE, MAKEPARAMS, WRITEPARFILEPP
%
% Last modified by Sirawich Pipatprathanporn, 11/10/2024

defval('branch', 'master')

if strcmp(branch(1:2), 'pp')
    % TODO: Implement writeparfile pp
    % writeparfilepp(params, fname, branch);
    return
end

if isempty(fname)
    fid = 1;
else
    fid = fopen(fname, 'w');
end

%% simulation input parameters (+Mesh for 'devel' branch)
writetitle(fid, 'simulation input parameters');
writeblank(fid);

writecomment(fid, '# title of job');
writestring(fid, 'title', params.title, []);
writeblank(fid);

writecomment(fid, '# forward or adjoint simulation');
writecomment(fid, '# 1 = forward, 2 = adjoint, 3 = both simultaneously');
writecomment(fid, ['# note: 2 is purposely UNUSED ' ...
    '(for compatibility with the numbering of our 3D codes)']);
writeint(fid, 'SIMULATION_TYPE', params.SIMULATION_TYPE, []);
writecomment(fid, ['# 0 = regular wave propagation simulation, ' ...
    '1/2/3 = noise simulation']);
writeint(fid, 'NOISE_TOMOGRAPHY', params.NOISE_TOMOGRAPHY, []);
writecomment(fid, '# save the last frame, needed for adjoint simulation');
writebool(fid, 'SAVE_FORWARD', params.SAVE_FORWARD, []);
writeblank(fid);

writecomment(fid, '# parameters concerning partitioning');
writeint(fid, 'NPROC', params.NPROC, 'number of processes');
if strcmpi(branch, 'master')
    % PARAMETER NAME CHANGE
    writeint(fid, 'partitioning_method', params.partitioning_method, ...
        'SCOTCH = 3, ascending order (very bad idea) = 1');
    writeblank(fid);

    % PARAMETER NAME CHANGE: ngnod is equivalent to NGNOD in 'devel' branch
    writecomment(fid, '# number of control nodes per element (4 or 9)');
    writeint(fid, 'ngnod', params.ngnod, []);
    writeblank(fid);

    writecomment(fid, '# time step parameters');
    writecomment(fid, '# total number of time steps');
    writeint(fid, 'NSTEP', params.NSTEP, []);
    writecomment(fid, ['# duration of a time step (see section "How to ' ...
        'choose the time step" of the manual for how to do this)']);
    writefloat(fid, 'DT', params.DT, []);
    writeblank(fid);

    writecomment(fid, '# time stepping');
    writecomment(fid, ['# 1 = Newmark (2nd order), ' ...
        '2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta), ' ... 
        '3 = classical RK4 4th-order 4-stage Runge-Kutta']);
    writeint(fid, 'time_stepping_scheme', params.time_stepping_scheme, []);
    writeblank(fid);

    writecomment(fid, ['# axisymmetric (2.5D) or ' ...
                 'Cartesian planar (2D) simulation']);
    writebool(fid, 'AXISYM', params.AXISYM, []);
    fprintf(fid, '\n');

    writecomment(fid, ['# set the type of calculation ' ...
                 '(P-SV or SH/membrane waves)']);
    writebool(fid, 'P_SV', params.P_SV, []);
    writeblank(fid);

    writecomment(fid, '# set to true to use GPUs');
    writebool(fid, 'GPU_MODE', params.GPU_MODE, []);
    writeblank(fid);
else
    writeblank(fid);
    
    writecomment(fid, '# time step parameters');
    writecomment(fid, '# total number of time steps');
    writeint(fid, 'NSTEP', params.NSTEP, []);
    writeblank(fid);
    
    writecomment(fid, ['# duration of a time step (see section "How to ' ...
        'choose the time step" of the manual for how to do this)']);
    writefloat(fid, 'DT', params.DT, []);
    writeblank(fid);

    writecomment(fid, '# time stepping');
    writecomment(fid, ['# 1 = Newmark (2nd order), ' ...
        '2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta), ' ... 
        '3 = classical RK4 4th-order 4-stage Runge-Kutta']);
    writeint(fid, 'time_stepping_scheme', params.time_stepping_scheme, []);
    writeblank(fid);
    
    writecomment(fid, ['# set the type of calculation ' ...
             '(P-SV or SH/membrane waves)']);
    writebool(fid, 'P_SV', params.P_SV, []);
    writeblank(fid);
    
    writecomment(fid, ['# axisymmetric (2.5D) or ' ...
                 'Cartesian planar (2D) simulation']);
    writebool(fid, 'AXISYM', params.AXISYM, []);
    writeblank(fid);
    
    writetitle(fid, 'Mesh');
    writeblank(fid);
    
    % PARAMETER NAME CHANGE
    writeint(fid, 'PARTITIONING_TYPE', params.partitioning_method, ...
        'SCOTCH = 3, ascending order (very bad idea) = 1');
    writeblank(fid);
    
    % PARAMETER NAME CHANGE: NGNOD is equivalent to ngnod in 'master' branch
    writecomment(fid, '# number of control nodes per element (4 or 9)');
    writeint(fid, 'NGNOD', params.ngnod, []);
    writeblank(fid);
end

writecomment(fid, ['# creates/reads a binary database that allows to ' ...
    'skip all time consuming setup steps in initialization']);
writecomment(fid, '# 0 = does not read/create database');
writecomment(fid, '# 1 = creates database');
writecomment(fid, '# 2 = reads database');
writeint(fid, 'setup_with_binary_database', ...
    params.setup_with_binary_database, []);
writeblank(fid);

writecomment(fid, '# available models');
writecomment(fid, '#   default       - define model using nbmodels below');
writecomment(fid, '#   ascii         - read model from ascii database file');
writecomment(fid, '#   binary        - read model from binary databse file');
writecomment(fid, '#   binary_voigt  - read Voigt model from binary database file');
writecomment(fid, '#   external      - define model using define_external_model subroutine');
writecomment(fid, '#   gll           - read GLL model from binary database file');
writecomment(fid, '#   legacy        - read model from model_velocity.dat_input');
writestring(fid, 'MODEL', params.MODEL, []);
writeblank(fid);

writecomment(fid, '# Output the model with the requested type, does not save if turn to default or .false.');
writecomment(fid, '# (available output formats: ascii,binary,gll,legacy)');
writestring(fid, 'SAVE_MODEL', params.SAVE_MODEL, []);
writeblank(fid);
writeblank(fid);

%% attenuation
writetitle(fid, 'attenuation');
writeblank(fid);

writecomment(fid, '# attenuation parameters');
writebool(fid, 'ATTENUATION_VISCOELASTIC', ...
    params.ATTENUATION_VISCOELASTIC, ...
    ['turn attenuation (viscoelasticity) on or off for ' ... 
    'non-poroelastic solid parts of the model']);
writebool(fid, 'ATTENUATION_VISCOACOUSTIC', ...
    params.ATTENUATION_VISCOACOUSTIC, ...
    ['turn attenuation (viscoacousticity) on or off for ' ... 
    'non-poroelastic fluid parts of the model']);
writeblank(fid);

writecomment(fid, '# for viscoelastic or viscoacoustic attenuation');
writeint(fid, 'N_SLS', params.N_SLS, ...
    ['number of standard linear solids for attenuation ' ... 
    '(3 is usually the minimum)']);
writefloat(fid, 'ATTENUATION_f0_REFERENCE', ...
    params.ATTENUATION_f0_REFERENCE, ...
    ['(Hz) relevant only if source is a Dirac or a Heaviside, ' ...
    'otherwise it is f0 the dominant frequency of the source in ' ...
    'the DATA/SOURCE file']);
writebool(fid, 'READ_VELOCITIES_AT_f0', params.READ_VELOCITIES_AT_f0, ...
    ['shift velocities to account for physical dispersion ' ...
    '(see user manual for more information)']);
writebool(fid, 'USE_SOLVOPT', params.USE_SOLVOPT, ...
    ['use more precise but much more expensive way of determining ' ...
    'the Q factor relaxation times, as in ' ...
    'https://doi.org/10.1093/gji/ggw024']);
writeblank(fid);

writecomment(fid, '# for poroelastic attenuation');
writebool(fid, 'ATTENUATION_PORO_FLUID_PART', params.ATTENUATION_PORO_FLUID_PART, ...
    ['turn viscous attenuation on or off for the fluid part of ' ...
    'poroelastic parts of the model']);
writefloat(fid, 'Q0_poroelastic', params.Q0_poroelastic, ...
    ['quality factor for viscous attenuation (ignore it if you are ' ...
    'not using a poroelastic material)']);
writefloat(fid, 'freq0_poroelastic', params.freq0_poroelastic, ...
    ['frequency for viscous attenuation (ignore it if you are not ' ...
    'using a poroelastic material)']);
writeblank(fid);

writecomment(fid, ['# to undo attenuation and/or PMLs for sensitivity ' ...
    'kernel calculations or forward runs with SAVE_FORWARD']);
writecomment(fid, ['# use the flag below. It performs undoing of ' ...
    'attenuation and/or of PMLs in an exact way for sensitivity ' ... 
    'kernel calculations']);
writecomment(fid, ['# but requires disk space for temporary storage, ' ...
    'and uses a significant amount of memory used as buffers for ' ... 
    'temporary storage.']);
writecomment(fid, ['# When that option is on the second parameter ' ...
    'indicates how often the code dumps restart files to disk ' ...
    '(if in doubt, use something between 100 and 1000).']);
writebool(fid, 'UNDO_ATTENUATION_AND_OR_PML', ...
    params.UNDO_ATTENUATION_AND_OR_PML, []);
writeint(fid, 'NT_DUMP_ATTENUATION', params.NT_DUMP_ATTENUATION , []);
writecomment(fid, ['# Instead of reconstructing the forward ' ...
    'wavefield, this option reads it from the disk using ' ...
    'asynchronous I/O.']);
writecomment(fid, ['Outperforms conventional mode using a value ' ...
    'of NSTEP_BETWEEN_COMPUTE_KERNELS high enough.']);
writebool(fid, 'NO_BACKWARD_RECONSTRUCTION', ...
    params.NO_BACKWARD_RECONSTRUCTION, []);
writeblank(fid);
writeblank(fid);

%% sources
writetitle(fid, 'sources');
writeblank(fid);

writecomment(fid, '# source parameters');
writeint(fid, 'NSOURCES', params.NSOURCES, ...
    ['number of sources (source information is then read from the ' ...
    'DATA/SOURCE file)']);
writebool(fid, 'force_normal_to_surface', ...
    params.force_normal_to_surface, ...
    'angleforce normal to surface (external mesh and curve file needed)');
writeblank(fid);

writecomment(fid, ['# use an existing initial wave field as source ' ...
    'or start from zero (medium initially at rest)']);
writebool(fid, 'initialfield', params.initialfield, []);
writebool(fid, 'add_Bielak_conditions_bottom', ...
    params.add_Bielak_conditions_bottom, ...
    'add Bielak conditions or not if initial plane wave');
writebool(fid, 'add_Bielak_conditions_right', ...
    params.add_Bielak_conditions_right, []);
writebool(fid, 'add_Bielak_conditions_top', ...
    params.add_Bielak_conditions_top, []);
writebool(fid, 'add_Bielak_conditions_left', ...
    params.add_Bielak_conditions_left, []);
writeblank(fid);

writecomment(fid, '# acoustic forcing');
writebool(fid, 'ACOUSTIC_FORCING', params.ACOUSTIC_FORCING , ...
    'acoustic forcing of an acoustic medium with a rigid interface');
writeblank(fid);

if strcmpi(branch, 'master')
    writeblank(fid);
else
    % NEW PARAMETER
    if ~isfield(params, 'noise_source_time_function_type')
        params.noise_source_time_function_type = 4;
    end
    writecomment(fid, ['# noise simulations - type of noise source ' ...
        'time function:']);
    writecomment(fid, ['# 0=external (S_squared), ' ...
        '1=Ricker(second derivative), 2=Ricker(first derivative), ' ...
        '3=Gaussian, 4=Figure 2a of Tromp et al. 2010']);
    writecomment(fid, ['# (default value 4 is chosen to reproduce the ' ...
        'time function from Fig 2a of "Tromp et al., 2010, ' ...
        'Noise Cross-Correlation Sensitivity Kernels")']);
    writeint(fid, 'noise_source_time_function_type', ...
        params.noise_source_time_function_type, []);
    writeblank(fid);
    
    % NEW PARAMETER
    if ~isfield(params, 'write_moving_sources_database')
        params.write_moving_sources_database = false;
    end
    writecomment(fid, '# moving sources');
    writecomment(fid, ['# Set write_moving_sources_database to .true. ' ...
        'if the generation of moving source databases takes']);
    writecomment(fid, ['# a long time. Then the simulation is done in ' ...
        'two steps: first you run the code and it writes the ' ...
        'databases to file']);
    writecomment(fid, ['# (in DATA folder by default). Then you rerun ' ...
        'the code and it will read the databases in there directly ' ...
        'possibly']);
    writecomment(fid, '# saving a lot of time.');
    writecomment(fid, '# This is only useful for GPU version (for now)');
    writebool(fid, 'write_moving_sources_database', ...
        params.write_moving_sources_database, []);
    writeblank(fid);
end

%% receivers
writetitle(fid, 'receivers');
writeblank(fid);

writecomment(fid, ['# receiver set parameters for recording ' ...
    'stations (i.e. recording points)']);
if strcmpi(branch, 'master')
    writeint(fid, 'seismotype', params.seismotype, ...
        ['record 1=displ 2=veloc 3=accel 4=pressure 5=curl of displ ' ...
        '6=the fluid potential']);
    writeblank(fid);

    % PARAMETER NAME CHANGE
    writecomment(fid, ['# subsampling of the seismograms to create ' ...
        'smaller files (but less accurately sampled in time)']);
    writeint(fid, 'subsamp_seismos', params.subsamp_seismos, []);
    writeblank(fid);

    writecomment(fid, ['# so far, this option can only be used if all the ' ...
        'receivers are in acoustic elements']);
    writebool(fid, 'USE_TRICK_FOR_BETTER_PRESSURE', ...
        params.USE_TRICK_FOR_BETTER_PRESSURE, []);
    writeblank(fid);

    writecomment(fid, ['# every how many time steps we save the ' ...
        'seismograms']);
    writecomment(fid, ['# (costly, do not use a very small value; if ' ...
        'you use a very large value that is larger than the total ' ...
        'number']);
    writecomment(fid, ['#  of time steps of the run, the seismograms ' ...
        'will automatically be saved once at the end of the run anyway)']);
    writeint(fid, 'NSTEP_BETWEEN_OUTPUT_SEISMOS ', ...
        params.NSTEP_BETWEEN_OUTPUT_SEISMOS, []);
    writeblank(fid);
else
    % TODO: make params.seismotype an array e.g. [1 2 4]
    writecomment(fid, ['# seismotype : record 1=displ 2=veloc 3=accel ' ...
        '4=pressure 5=curl of displ 6=the fluid potential']);
    % convert params.seismotype (int array) to string
    str = indeks(sprintf('%s,', string(params.seismotype)), '1:end-1');
    writestring(fid, 'seismotype', str, ...
        'several values can be chosen. For example : 1,2,4');
    writeblank(fid);
    
    writecomment(fid, ['# interval in time steps for writing of ' ...
        'seismograms']);
    writecomment(fid, ['# every how many time steps we save the ' ...
        'seismograms']);
    writecomment(fid, ['# (costly, do not use a very small value; if ' ...
        'you use a very large value that is larger than the total ' ...
        'number']);
    writecomment(fid, ['#  of time steps of the run, the seismograms ' ...
        'will automatically be saved once at the end of the run anyway)']);
    writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS ', ...
        params.NSTEP_BETWEEN_OUTPUT_SEISMOS, []);
    writeblank(fid);
    
    % PARAMETER NAME CHANGE
    writecomment(fid, ['# set to n to reduce the sampling rate of ' ...
        'output seismograms by a factor of n']);
    writecomment(fid, '# defaults to 1, which means no down-sampling');
    writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_SAMPLE', ...
        params.subsamp_seismos, []);
    writeblank(fid);
    
    writecomment(fid, ['# so far, this option can only be used if all the ' ...
        'receivers are in acoustic elements']);
    writebool(fid, 'USE_TRICK_FOR_BETTER_PRESSURE', ...
        params.USE_TRICK_FOR_BETTER_PRESSURE, []);
    writeblank(fid);
end

writecomment(fid, ['# use this t0 as earliest starting time rather ' ...
    'than the automatically calculated one']);
writefloat(fid, 'USER_T0', params.USER_T0, []);
writeblank(fid);

writecomment(fid, '# seismogram formats');
writebool(fid, 'save_ASCII_seismograms', params.save_ASCII_seismograms, ...
    'save seismograms in ASCII format or not');
writebool(fid, 'save_binary_seismograms_single', ...
    params.save_binary_seismograms_single, ...
    ['save seismograms in single precision binary format or not ' ...
    '(can be used jointly with ASCII above to save both)']);
writebool(fid, 'save_binary_seismograms_double', ...
    params.save_binary_seismograms_double, ...
    ['save seismograms in double precision binary format or not ' ...
    '(can be used jointly with both flags above to save all)']);
writebool(fid, 'SU_FORMAT', params.SU_FORMAT, ...
    ['output single precision binary seismograms in Seismic Unix ' ...
    'format (adjoint traces will be read in the same format)']);
writeblank(fid);

writecomment(fid, ['# use an existing STATION file found in ./DATA or ' ...
    'create a new one from the receiver positions below in this ' ...
    'Par_file']);
writebool(fid, 'use_existing_STATIONS', params.use_existing_STATIONS, []);
writeblank(fid);

writecomment(fid, ['# number of receiver sets (i.e. number of ' ...
    'receiver lines to create below)']);
writeint(fid, 'nreceiversets', params.nreceiversets, []);
writeblank(fid);

writecomment(fid, '# orientation');
writefloat(fid, 'anglerec', params.anglerec, ...
    'angle to rotate components at receivers');
writebool(fid, 'rec_normal_to_surface', params.rec_normal_to_surface, ...
    ['base anglerec normal to surface (external mesh and curve file ' ...
    'needed)']);
writeblank(fid);

for ii = 1:params.nreceiversets
    % comment labeling n-th receiver set
    writecomment(fid, sprintf(['# %s receiver set (repeat these 6 ' ...
        'lines and adjust nreceiversets accordingly)'], orderize(ii)));
    writeint(fid, 'nrec', params.RECEIVERS{ii}.nrec, ...
        'number of receivers');
    writefloat(fid, 'xdeb', params.RECEIVERS{ii}.xdeb, ...
        'first receiver x in meters');
    writefloat(fid, 'zdeb', params.RECEIVERS{ii}.zdeb, ...
        'first receiver z in meters');
    writefloat(fid, 'xfin', params.RECEIVERS{ii}.xfin, ...
        'last receiver x in meters (ignored if only one receiver)');
    writefloat(fid, 'zfin', params.RECEIVERS{ii}.zfin, ...
        'last receiver z in meters (ignored if only one reciever)');
    writebool(fid, 'record_at_surface_same_vertical', ...
        params.RECEIVERS{ii}.record_at_surface_same_vertical, ...
        ['receivers inside the medium or at the surface (z values are ' ...
        'ignored if this is set to true, they are replaced with the ' ...
        'topography height)']);
    writeblank(fid);
end
writeblank(fid);

%% adjoint kenel outputs
writetitle(fid, 'adjoint kernel outputs');
writeblank(fid);

writecomment(fid, ['# save sensitivity kernels in ASCII format ' ...
    '(much bigger files, but compatible with current GMT scripts) or ' ...
    'in binary format']);
writebool(fid, 'save_ASCII_kernels', params.save_ASCII_kernels, []);
writeblank(fid);

if strcmpi(branch, 'master')
    % PARAMETER NAME CHANGE:
    writecomment(fid, ['# since the accuracy of kernel integration may ' ...
        'not need to respect the CFL, this option permits to save ' ...
        'computing time, and memory with UNDO_ATTENUATION_AND_OR_PML mode']);
    writeint(fid, 'NSTEP_BETWEEN_COMPUTE_KERNELS', ...
        params.NSTEP_BETWEEN_COMPUTE_KERNELS, []);
    writeblank(fid);
    writeblank(fid);
else
    % PARAMETER NAME CHANGE:
    writecomment(fid, ['# since the accuracy of kernel integration may ' ...
        'not need to respect the CFL, this option permits to save ' ...
        'computing time, and memory with UNDO_ATTENUATION_AND_OR_PML mode']);
    writeint(fid, 'NTSTEP_BETWEEN_COMPUTE_KERNELS', ...
        params.NSTEP_BETWEEN_COMPUTE_KERNELS, []);
    writeblank(fid);
    
    % NEW PARAMETER
    if ~isfield(params, 'APPROXIMATE_HESS_KL')
        params.APPROXIMATE_HESS_KL = false;
    end
    writecomment(fid, '# outputs approximate Hessian for preconditioning');
    writebool(fid, 'APPROXIMATE_HESS_KL', params.APPROXIMATE_HESS_KL, []);
    writeblank(fid);
end

%% boundary conditions
writetitle(fid, 'boundary conditions');
writeblank(fid);

writecomment(fid, '# Perfectly Matched Layer (PML) boundaries');
writecomment(fid, '# absorbing boundary active or not');
writebool(fid, 'PML_BOUNDARY_CONDITIONS', ...
    params.PML_BOUNDARY_CONDITIONS, []);
writeint(fid, 'NELEM_PML_THICKNESS', params.NELEM_PML_THICKNESS, []);
writebool(fid, 'ROTATE_PML_ACTIVATE', params.ROTATE_PML_ACTIVATE, []);
writefloat(fid, 'ROTATE_PML_ANGLE', params.ROTATE_PML_ANGLE, []);
writecomment(fid, ['# change the four parameters below only if you ' ...
    'know what you are doing; they change the damping profiles inside ' ...
    'the PMLs']);
writefloat(fid, 'K_MIN_PML', params.K_MIN_PML, 'from Gedney page 8.11');
writefloat(fid, 'K_MAX_PML', params.K_MAX_PML, []);
writefloat(fid, 'damping_change_factor_acoustic', ...
    params.damping_change_factor_acoustic, []);
writefloat(fid, 'damping_change_factor_elastic', ...
    params.damping_change_factor_elastic, []);
writecomment(fid, ['# set the parameter below to .false. unless you ' ...
    'know what you are doing; this implements automatic adjustment of ' ...
    'the PML parameters for elongated models.']);
writecomment(fid, ['# The goal is to improve the absorbing efficiency ' ...
    'of PML for waves with large incidence angles, but this can lead ' ...
    'to artefacts.']);
writecomment(fid, ['# In particular, this option is efficient only ' ...
    'when the number of sources NSOURCES is equal to one.']);
writebool(fid, 'PML_PARAMETER_ADJUSTMENT', ...
    params.PML_PARAMETER_ADJUSTMENT, []);
writeblank(fid);

writecomment(fid, '# Stacey ABC');
writebool(fid, 'STACEY_ABSORBING_CONDITIONS', ...
    params.STACEY_ABSORBING_CONDITIONS, []);
writeblank(fid);

writecomment(fid, '# periodic boundaries');
writebool(fid, 'ADD_PERIODIC_CONDITIONS', ...
    params.ADD_PERIODIC_CONDITIONS, []);
writefloat(fid, 'PERIODIC_HORIZ_DIST', params.PERIODIC_HORIZ_DIST, []);
writeblank(fid);

%% velocity and density moels
writetitle(fid, 'velocity and density models');
writeblank(fid);

writecomment(fid, '# number of model materials');
writeint(fid, 'nbmodels', params.nbmodels, []);
if strcmpi(branch, 'master')
    writecomment(fid, '# available material types (see user manual for more information)');
    writecomment(fid, '#   acoustic:    model_number 1 rho Vp 0  0 0 QKappa Qmu 0 0 0 0 0 0 (for QKappa and Qmu use 9999 to ignore them)');
    writecomment(fid, '# when viscoelasticity is turned on, the Vp and Vs values that are read here are the UNRELAXED ones i.e. the values at infinite frequency');
    writecomment(fid, '# unless the READ_VELOCITIES_AT_f0 parameter above is set to true, in which case they are the values at frequency f0.');
    writecomment(fid, '# Please also note that Qmu is always equal to Qs, but Qkappa is in general not equal to Qp.');
    writecomment(fid, '# To convert one to the other see doc/Qkappa_Qmu_versus_Qp_Qs_relationship_in_2D_plane_strain.pdf and');
    writecomment(fid, '# utils/attenuation/conversion_from_Qkappa_Qmu_to_Qp_Qs_from_Dahlen_Tromp_959_960.f90.');
    writecomment(fid, '#   elastic:     model_number 1 rho Vp Vs 0 0 QKappa Qmu 0 0 0 0 0 0 (for QKappa and Qmu use 9999 to ignore them)');
    writecomment(fid, '#   anisotropic: model_number 2 rho c11 c13 c15 c33 c35 c55 c12 c23 c25 0 0 0');
    writecomment(fid, '#   poroelastic: model_number 3 rhos rhof phi c kxx kxz kzz Ks Kf Kfr etaf mufr Qmu');
    writecomment(fid, '#   tomo:        model_number -1 0 0 A 0 0 0 0 0 0 0 0 0 0');
else
    writecomment(fid, '# available material types (see user manual for more information)');
    writecomment(fid, '#   acoustic:              model_number 1 rho Vp 0  0 0 QKappa 9999 0 0 0 0 0 0 (for QKappa use 9999 to ignore them)');
    writecomment(fid, '#   elastic:               model_number 1 rho Vp Vs 0 0 QKappa Qmu  0 0 0 0 0 0 (for QKappa and Qmu use 9999 to ignore them)');
    writecomment(fid, '#   anisotropic:           model_number 2 rho c11 c13 c15 c33 c35 c55 c12 c23 c25   0 QKappa Qmu');
    writecomment(fid, '#   anisotropic in AXISYM: model_number 2 rho c11 c13 c15 c33 c35 c55 c12 c23 c25 c22 QKappa Qmu');
    writecomment(fid, '#   poroelastic:           model_number 3 rhos rhof phi c kxx kxz kzz Ks Kf Kfr etaf mufr Qmu');
    writecomment(fid, '#   tomo:                  model_number -1 0 0 A 0 0 0 0 0 0 0 0 0 0');
    writecomment(fid, '#');
    writecomment(fid, '# note: When viscoelasticity or viscoacousticity is turned on,');
    writecomment(fid, '#       the Vp and Vs values that are read here are the UNRELAXED ones i.e. the values at infinite frequency');
    writecomment(fid, '#       unless the READ_VELOCITIES_AT_f0 parameter above is set to true, in which case they are the values at frequency f0.');
    writecomment(fid, '#');
    writecomment(fid, '#       Please also note that Qmu is always equal to Qs, but QKappa is in general not equal to Qp.');
    writecomment(fid, '#       To convert one to the other see doc/Qkappa_Qmu_versus_Qp_Qs_relationship_in_2D_plane_strain.pdf and');
    writecomment(fid, '#       utils/attenuation/conversion_from_Qkappa_Qmu_to_Qp_Qs_from_Dahlen_Tromp_959_960.f90.');
end

for ii = 1:params.nbmodels
    model = params.MODELS{ii};
    switch model.type_number
        case 1 % acoustic / elastic
            if or(model.QKappa == 9999, model.Qmu == 9999)
                s = sprintf('%d 1 %.3e %.3e %.3e 0 0 9999 9999 0 0 0 0 0 0 \n', ...
                    ii, model.rho, model.vp, model.vs);
            else
                s = sprintf('%d 1 %.3e %.3e %.3e 0 0 %.3e %.3e 0 0 0 0 0 0 \n', ...
                    ii, model.rho, model.vp, model.vs, model.QKappa, model.Qmu);
            end
        case 2 % anisotropic
            s = sprintf('%d 2 %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e 0 0\n', ...
                ii, model.rho, model.c11, model.c13, model.c15, ...
                model.c33, model.c35, model.c55, model.c12, model.c23, ...
                model.c25, model.c22);
        case 3 % poroelastic
            s = sprintf('%d 3 %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %3.e\n', ...
                ii, model.rhos, model.rhof, model.phi, model.c, ...
                model.kxx, model.kxz, model.kzz, model.Ks, model.Kf, ....
                model.etaf, model.mufr, model.Qmu);
        case -1 % tomo
            s = sprintf('%d -1 0 0 %.3e 0 0 0 0 0 0 0 0 0 0\n', ii, model.A);
        otherwise
            fprintf('Invalid model: cannot write into Par_file. Skip!\n');
            s = [];
    end
    if ~isempty(s)
        % don't forget to replace e to d for exponent notation
        s = replace(s, 'e', 'd');
        fprintf(fid, s);
    end
end
writeblank(fid);

writecomment(fid, '# external tomography file');
writestring(fid, 'TOMOGRAPHY_FILE', params.TOMOGRAPHY_FILE, []);
writeblank(fid);

writecomment(fid, ['# use an external mesh created by an external ' ...
    'meshing tool or use the internal mesher']);
writebool(fid, 'read_external_mesh', params.read_external_mesh, []);
writeblank(fid);

%% PARAMETERS FOR EXTERNAL MESHING
writetitle(fid, 'PARAMETERS FOR EXTERNAL MESHING');
writeblank(fid);

writecomment(fid, ['# data concerning mesh, when generated using ' ...
    'third-party app (more info in README)']);
writecomment(fid, '# (see also absorbing_conditions above)');
writestring(fid, 'mesh_file', params.mesh_file, 'file containing the mesh');
writestring(fid, 'nodes_coords_file', params.nodes_coords_file, ...
    'file containing the nodes coordinates');
writestring(fid, 'materials_file', params.materials_file, ...
    'file containing the material number for each element');
writestring(fid, 'free_surface_file', params.free_surface_file, ...
    'file containing the free surface');
writestring(fid, 'axial_elements_file', params.axial_elements_file, ...
    'file containing the axial elements if AXISYM is true');
writestring(fid, 'absorbing_surface_file', ...
    params.absorbing_surface_file, ...
    'file containing the absorbing surface');
% This parameter is not used in many examples but a few ones
if isfield(params, 'params.CPML_element_file')
    writestring(fid, 'CPML_element_file', params.CPML_element_file, ...
        'file containing the CPML element numbers');
end
writestring(fid, 'acoustic_forcing_surface_file', ...
    params.acoustic_forcing_surface_file, ...
    'file containing the acoustic forcing surface');
writestring(fid, 'absorbing_cpml_file', params.absorbing_cpml_file, ...
    'file containing the CPML element numbers');
writestring(fid, 'tangential_detection_curve_file', ...
    params.tangential_detection_curve_file, ...
    'file containing the curve delimiting the velocity model');
writeblank(fid);

%% PARAMETERS FOR INTERNAL MESHING
writetitle(fid, 'PARAMETERS FOR INTERNAL MESHING');
writeblank(fid);

writecomment(fid, '# file containing interfaces for internal mesh');
writestring(fid, 'interfacesfile', params.interfacesfile, []);
writeblank(fid);

writecomment(fid, ['# geometry of the model ' ...
    '(origin lower-left corner = 0,0) and mesh description']);
writefloat(fid, 'xmin', params.xmin, 'abscissa of left side of the model');
writefloat(fid, 'xmax', params.xmax, 'abscissa of right side of the model');
writeint(fid, 'nx', params.nx, 'number of elements along X');
writeblank(fid);

writecomment(fid, ['# absorbing boundary parameters ' ...
    '(see absorbing_conditions above)']);
writebool(fid, 'absorbbottom', params.absorbbottom, []);
writebool(fid, 'absorbright', params.absorbright, []);
writebool(fid, 'absorbtop', params.absorbtop, []);
writebool(fid, 'absorbleft', params.absorbleft, []);
writeblank(fid);

writecomment(fid, ['# define the different regions of the model in ' ...
    'the (nx,nz) spectral-element mesh']);
writeint(fid, 'nbregions', params.nbregions, ...
    'then set below the different regions and model number for each region');
writecomment(fid, ['# format of each line: nxmin nxmax nzmin nzmax ' ...
    'material_number']);

for ii = 1:params.nbregions
    region = params.REGIONS{ii};
    fprintf(fid, '%d %d %d %d %d\n', region.nxmin, region.nxmax, ...
        region.nzmin, region.nzmax, region.material_number);
end
writeblank(fid);

%% display parameters
writetitle(fid, 'display parameters');
writeblank(fid);

writecomment(fid, ['# interval at which we output time step info and ' ...
    'max of norm of displacement']);
writecomment(fid, ['# every how many time steps we display ' ...
    'information about the simulation (costly, do not use a very ' ...
    'small value)']);
% PARAMETER NAME CHANGE
if strcmpi(branch, 'master')
    writeint(fid, 'NSTEP_BETWEEN_OUTPUT_INFO', ...
        params.NSTEP_BETWEEN_OUTPUT_INFO, []);
else
    writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_INFO', ...
        params.NSTEP_BETWEEN_OUTPUT_INFO, []);
end
writeblank(fid);

writecomment(fid, '# meshing output');
writebool(fid, 'output_grid_Gnuplot', params.output_grid_Gnuplot, ...
    'generate a GNUPLOT file containing the grid, and a script to plot it');
writebool(fid, 'output_grid_ASCII', params.output_grid_ASCII, ...
    ['dump the grid in an ASCII text file consisting of a set of '...
    'X,Y,Z points or not']);
writeblank(fid);

writecomment(fid, ['# to plot total energy curves, for instance to ' ...
    'monitor how CPML absorbing layers behave;']);
writecomment(fid, ['# should be turned OFF in most cases because ' ...
    'a bit expensive']);
writebool(fid, 'OUTPUT_ENERGY', params.OUTPUT_ENERGY, []);
writeblank(fid);

writecomment(fid, ['# every how many time steps we compute energy ' ...
    '(which is a bit expensive to compute)']);
writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_ENERGY', ...
    params.NTSTEP_BETWEEN_OUTPUT_ENERGY, []);
writeblank(fid);

writecomment(fid, ['# Compute the field int_0^t v^2 dt for a set of ' ...
    'GLL points and write it to file. Use']);
writecomment(fid, ...
    ['# the script utils/visualisation/plotIntegratedEnergyFile.py ' ...
    'to watch. It is refreshed at the same time than the seismograms']);
writebool(fid, 'COMPUTE_INTEGRATED_ENERGY_FIELD', ...
    params.COMPUTE_INTEGRATED_ENERGY_FIELD, []);
writeblank(fid);

%% movies/images/snapshots
writetitle(fid, 'movies/images/snapshots');
writeblank(fid);

writecomment(fid, ['# every how many time steps we draw JPEG or ' ...
    'PostScript pictures of the simulation']);
writecomment(fid, ['# and/or we dump results of the simulation as ' ...
    'ASCII or binary files (costly, do not use a very small value)']);
if strcmpi(branch, 'master')
    % PARAMETER NAME CHANGE
    writeint(fid, 'NSTEP_BETWEEN_OUTPUT_IMAGES', ...
        params.NSTEP_BETWEEN_OUTPUT_IMAGES, []);
else
    % PARAMETER NAME CHANGE
    writeint(fid, 'NTSTEP_BETWEEN_OUTPUT_IMAGES', ...
        params.NSTEP_BETWEEN_OUTPUT_IMAGES, []);
end
writeblank(fid);

writecomment(fid, ['# minimum amplitude kept in % for the JPEG and ' ...
    'PostScript snapshots; amplitudes below that are muted']);
writefloat(fid, 'cutsnaps', params.cutsnaps, []);
writeblank(fid);

writecomment(fid, '#### for JPEG color images ####');
if strcmpi(branch, 'master')
    writebool(fid, 'output_color_image', params.output_color_image, ...
        ['output JPEG color image of the results every ' ...
        'NSTEP_BETWEEN_OUTPUT_IMAGES time steps or not']);
else
    writebool(fid, 'output_color_image', params.output_color_image, ...
        ['output JPEG color image of the results every ' ...
        'NTSTEP_BETWEEN_OUTPUT_IMAGES time steps or not']);
end
writeint(fid, 'imagetype_JPEG', params.imagetype_JPEG, ...
    ['display 1=displ_Ux 2=displ_Uz 3=displ_norm 4=veloc_Vx ' ...
    '5=veloc_Vz 6=veloc_norm 7=accel_Ax 8=accel_Az 9=accel_norm ' ...
    '10=pressure']);
writefloat(fid, 'factor_subsample_image', params.factor_subsample_image, ...
    ['(double precision) factor to subsample or oversample ' ...
    '(if set to e.g. 0.5) color images output by the code ' ...
    '(useful for very large models, or to get nicer looking denser ' ...
    'pictures)']);
writebool(fid, 'USE_CONSTANT_MAX_AMPLITUDE', ...
    params.USE_CONSTANT_MAX_AMPLITUDE, ...
    ['by default the code normalizes each image independently to its ' ...
    'maximum; use this option to use the global maximum below instead']);
writefloat(fid, 'CONSTANT_MAX_AMPLITUDE_TO_USE', ...
    params.CONSTANT_MAX_AMPLITUDE_TO_USE, ...
    ['constant maximum amplitude to use for all color images if the ' ...
    'above USE_CONSTANT_MAX_AMPLITUDE option is true']);
writefloat(fid, 'POWER_DISPLAY_COLOR', params.POWER_DISPLAY_COLOR, ...
    'non linear display to enhance small amplitudes in JPEG color images');
writebool(fid, 'DRAW_SOURCES_AND_RECEIVERS', ...
    params.DRAW_SOURCES_AND_RECEIVERS, ...
    ['display sources as orange crosses and receivers as green ' ...
    'squares in JPEG images or not']);
writebool(fid, 'DRAW_WATER_IN_BLUE', params.DRAW_WATER_IN_BLUE, ...
    ['display acoustic layers as constant blue in JPEG images, ' ...
    'because they likely correspond to water in the case of ocean ' ...
    'acoustics or in the case of offshore oil industry experiments ' ...
    '(if off, display them as greyscale, as for elastic or ' ...
    'poroelastic elements, for instance for acoustic-only oil ' ...
    'industry models of solid media)']);
writebool(fid, 'USE_SNAPSHOT_NUMBER_IN_FILENAME', ...
    params.USE_SNAPSHOT_NUMBER_IN_FILENAME, ...
    ['use snapshot number in the file name of JPEG color snapshots ' ...
    'instead of the time step (for instance to create movies in an ' ...
    'easier way later)']);
writeblank(fid);

writecomment(fid, '#### for PostScript snapshots ####');
if strcmpi(branch, 'master')
    writebool(fid, 'output_postscript_snapshot', ...
        params.output_postscript_snapshot, ...
        ['output Postscript snapshot of the results every ' ...
        'NSTEP_BETWEEN_OUTPUT_IMAGES time steps or not']);
else
    writebool(fid, 'output_postscript_snapshot', ...
        params.output_postscript_snapshot, ...
        ['output Postscript snapshot of the results every ' ...
        'NTSTEP_BETWEEN_OUTPUT_IMAGES time steps or not']);
end
writeint(fid, 'imagetype_postscript', params.imagetype_postscript, ...
    ['display 1=displ vector 2=veloc vector 3=accel vector; small ' ...
    'arrows are displayed for the vectors']);
writebool(fid, 'meshvect', params.meshvect, ...
    'display mesh on PostScript plots or not');
writebool(fid, 'modelvect', params.modelvect, ...
    'display velocity model on PostScript plots or not');
writebool(fid, 'boundvect', params.boundvect, ...
    'display boundary conditions on PostScript plots or not');
writebool(fid, 'interpol', params.interpol, ...
    ['interpolation of the PostScript display on a regular grid ' ...
    'inside each spectral element, or use the non-evenly spaced GLL ' ...
    'points']);
writeint(fid, 'pointsdisp', params.pointsdisp, ...
    ['number of points in each direction for interpolation of ' ...
    'PostScript snapshots (set to 1 for lower-left corner only)']);
writeint(fid, 'subsamp_postscript', params.subsamp_postscript, ...
    'subsampling of background velocity model in PostScript snapshots');
writefloat(fid, 'sizemax_arrows', params.sizemax_arrows, ...
    'maximum size of arrows on PostScript plots in centimeters');
writebool(fid, 'US_LETTER', params.US_LETTER, ...
    'use US letter or European A4 paper for PostScript plots');
writeblank(fid);

writecomment(fid, '#### for wavefield dumps ####');
writebool(fid, 'output_wavefield_dumps', params.output_wavefield_dumps, ...
    'output wave field to a text file (creates very big files)');
writeint(fid, 'imagetype_wavefield_dumps', ...
    params.imagetype_wavefield_dumps, ...
    'display 1=displ vector 2=veloc vector 3=accel vector 4=pressure');
writebool(fid, 'use_binary_for_wavefield_dumps', ...
    params.use_binary_for_wavefield_dumps, ...
    ['use ASCII or single-precision binary format for the wave ' ...
    'field dumps']);
writeblank(fid);

writecomment(fid, '#-----------------------------------------------------------');
writeblank(fid);

very_long_comment1 = {'# Ability to run several calculations (several earthquakes)', ...
    '# in an embarrassingly-parallel fashion from within the same run;', ...
    '# this can be useful when using a very large supercomputer to compute', ...
    '# many earthquakes in a catalog, in which case it can be better from', ...
    '# a batch job submission point of view to start fewer and much larger jobs,', ...
    '# each of them computing several earthquakes in parallel.', ...
    '# To turn that option on, set parameter NUMBER_OF_SIMULTANEOUS_RUNS to a value greater than 1.', ...
    '# To implement that, we create NUMBER_OF_SIMULTANEOUS_RUNS MPI sub-communicators,', ...
    '# each of them being labeled "my_local_mpi_comm_world", and we use them', ...
    '# in all the routines in "src/shared/parallel.f90", except in MPI_ABORT() because in that case', ...
    '# we need to kill the entire run.', ...
    '# When that option is on, of course the number of processor cores used to start', ...
    '# the code in the batch system must be a multiple of NUMBER_OF_SIMULTANEOUS_RUNS,', ...
    '# all the individual runs must use the same number of processor cores,', ...
    '# which as usual is NPROC in the Par_file,', ...
    '# and thus the total number of processor cores to request from the batch system', ...
    '# should be NUMBER_OF_SIMULTANEOUS_RUNS * NPROC.', ...
    '# All the runs to perform must be placed in directories called run0001, run0002, run0003 and so on', ...
    '# (with exactly four digits).', ...
    '#', ...
    '# Imagine you have 10 independent calculations to do, each of them on 100 cores; you have three options:', ...
    '#', ...
    '# 1/ submit 10 jobs to the batch system', ...
    '#', ...
    '# 2/ submit a single job on 1000 cores to the batch, and in that script create a sub-array of jobs to start 10 jobs,', ...
    '# each running on 100 cores (see e.g. http://www.schedmd.com/slurmdocs/job_array.html )', ...
    '#', ...
    '# 3/ submit a single job on 1000 cores to the batch, start SPECFEM2D on 1000 cores, create 10 sub-communicators,', ...
    '# cd into one of 10 subdirectories (called e.g. run0001, run0002,... run0010) depending on the sub-communicator', ...
    '# your MPI rank belongs to, and run normally on 100 cores using that sub-communicator.', ...
    '#', ...
    '# The option below implements 3/.', ...
    '#'};

for ii = 1:length(very_long_comment1)
    writecomment(fid, very_long_comment1{ii});
end
writeint(fid, 'NUMBER_OF_SIMULTANEOUS_RUNS', ...
    params.NUMBER_OF_SIMULTANEOUS_RUNS, []);
writeblank(fid);

very_long_comment2 = {'# if we perform simultaneous runs in parallel, if only the source and receivers vary between these runs', ...
    '# but not the mesh nor the model (velocity and density) then we can also read the mesh and model files', ...
    '# from a single run in the beginning and broadcast them to all the others; for a large number of simultaneous', ...
    '# runs for instance when solving inverse problems iteratively this can DRASTICALLY reduce I/Os to disk in the solver', ...
    '# (by a factor equal to NUMBER_OF_SIMULTANEOUS_RUNS), and reducing I/Os is crucial in the case of huge runs.', ...
    '# Thus, always set this option to .true. if the mesh and the model are the same for all simultaneous runs.', ...
    '# In that case there is no need to duplicate the mesh and model file database (the content of the DATABASES_MPI', ...
    '# directories) in each of the run0001, run0002,... directories, it is sufficient to have one in run0001', ...
    '# and the code will broadcast it to the others)'};

for ii = 1:length(very_long_comment2)
    writecomment(fid, very_long_comment2{ii});
end
writebool(fid, 'BROADCAST_SAME_MESH_AND_MODEL', ...
    params.BROADCAST_SAME_MESH_AND_MODEL, []);
writeblank(fid);

if strcmpi(branch, 'devel')
    writecomment(fid, '#-----------------------------------------------------------------------------');
    writeblank(fid);
    
    writecomment(fid, '# set to true to use GPUs');
    writebool(fid, 'GPU_MODE', params.GPU_MODE, []);
end

%% close the file
if fid >= 3
    fclose(fid);
end
end

%% Helper functions
function writetitle(fid, title)
fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#\n');
fprintf(fid, '# %s\n', title);
fprintf(fid, '#\n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');
end

function writeblank(fid)
fprintf(fid, '\n');
end

function writecomment(fid, comment)
if isempty(comment)
    writeblank(fid);
else
    if ~strcmp(comment(1), '#')
        fprintf(fid, '# %s\n', comment);
    else
        fprintf(fid, '%s\n', comment);
    end
end
end

function writebool(fid, name, value, comment)
if value == 0
    var_string = '.false.';
else
    var_string = '.true.';
end
writestring(fid, name, var_string, comment);
end


function writefloat(fid, name, value, comment)
if value == 0
    var_string = '0';
elseif and(abs(value) >= 1, abs(value) < 100)
    var_string = sprintf('%.3f', value);
else
    var_string = sprintf('%.3e', value);
    var_string = replace(var_string, 'e', 'd');
end
writestring(fid, name, var_string, comment);
end

function writeint(fid, name, value, comment)
if isempty(comment)
    fprintf(fid, '%-31s = %d\n', name, value);
else
    fprintf(fid, '%-31s = %-14d # %s\n', name, value, comment);
end
end

function writestring(fid, name, value, comment)
if isempty(comment)
    fprintf(fid, '%-31s = %s\n', name, value);
else
    fprintf(fid, '%-31s = %-14s # %s\n', name, value, comment);
end
end

% turn numeric to order e.g. 1 --> first, 2 --> second etc.
function order = orderize(number)
switch number
    case 1
        order = 'first';
    case 2
        order = 'second';
    case 3
        order = 'third';
    case 4
        order = 'fourth';
    case 5
        order = 'fifth';
    case 6
        order = 'sixth';
    case 7
        order = 'seventh';
    case 8
        order = 'eighth';
    case 9
        order = 'ninth';
    case 10
        order = 'tenth';
    otherwise
        order = char(string(number));
        if and(mod(number, 100) >= 11, mod(number, 100) <= 19)
            order = strcat(order, 'th');
        elseif mod(number, 10) == 1
            order = strcat(order, 'st');
        elseif mod(number, 10) == 2
            order = strcat(order, 'nd');
        elseif mod(number, 10) == 3
            order = strcat(order, 'rd');
        else
            order = strcat(order, 'th');
        end
end
end