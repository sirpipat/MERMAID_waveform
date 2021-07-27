function params = makeparams()
% params = MAKEPARAMS
%
% Makes a sensible parameters for Par_file.
%
% SEE ALSO:
% WRITEPARFILE
%
% Last modified by Sirawich Pipatprathanporn, 07/22/2021

% a default set of receivers
receiverset1 = struct(...
    'nrec'                              , 11    , ...
    'xdeb'                              , 300   , ...
    'zdeb'                              , 2200  , ...
    'xfin'                              , 3700  , ...
    'zfin'                              , 2200  , ...
    'record_at_surface_same_vertical'   , true    ...
);

receiverset2 = struct(...
    'nrec'                              , 11    , ...
    'xdeb'                              , 2500  , ...
    'zdeb'                              , 2500  , ...
    'xfin'                              , 2500  , ...
    'zfin'                              , 0     , ...
    'record_at_surface_same_vertical'   , false   ...
);

receiversets = {{receiverset1, receiverset2}};

% default models
model1 = struct(...
    'model_number'      , 1         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'elastic' , ...
    'rho'               , 2700      , ...
    'vp'                , 3000      , ...
    'vs'                , 1732.051  , ...
    'QKappa'           , 9999      , ...
    'Qmu'               , 9999        ...
);

model2 = struct(...
    'model_number'      , 2         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'acoustic', ...
    'rho'               , 2500      , ...
    'vp'                , 2700      , ...
    'vs'                , 0         , ...
    'QKappa'           , 9999      , ...
    'Qmu'               , 9999        ...
);

model3 = struct(...
    'model_number'      , 3         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'elastic' , ...
    'rho'               , 2200      , ...
    'vp'                , 2500      , ...
    'vs'                , 1443.375  , ...
    'QKappa'           , 9999      , ...
    'Qmu'               , 9999        ...
);

model4 = struct(...
    'model_number'      , 4         , ...
    'type_number'       , 1         , ...
    'type_name'         , 'elastic' , ...
    'rho'               , 2200      , ...
    'vp'                , 2200      , ...
    'vs'                , 1343.375  , ...
    'QKappa'           , 9999      , ...
    'Qmu'               , 9999        ...
);

models = {{model1, model2, model3, model4}};

% default regions
region1 = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , 80        , ...
    'nzmin'             , 1         , ...
    'nzmax'             , 20        , ...
    'material_number'   , 1           ...
);

region2 = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , 80        , ...
    'nzmin'             , 21        , ...
    'nzmax'             , 40        , ...
    'material_number'   , 2           ...
);

region3 = struct(...
    'nxmin'             , 1         , ...
    'nxmax'             , 80        , ...
    'nzmin'             , 41        , ...
    'nzmax'             , 60        , ...
    'material_number'   , 3           ...
);

region4 = struct(...
    'nxmin'             , 60        , ...
    'nxmax'             , 70        , ...
    'nzmin'             , 21        , ...
    'nzmax'             , 40        , ...
    'material_number'   , 4           ...
);

regions = {{region1, region2, region3, region4}};
    
params = struct(...
    'title'                          , 'DEFAULT TITLE'              , ...
    'SIMULATION_TYPE'                , 1                            , ...
    'NOISE_TOMOGRAPHY'               , 0                            , ...
    'SAVE_FORWARD'                   , false                        , ...
    'NPROC'                          , 1                            , ...
    'partitioning_method'            , 3                            , ...
    'ngnod'                          , 4                            , ...
    'NSTEP'                          , 2000                         , ...
    'DT'                             , 1.1e-3                       , ...
    'time_stepping_scheme'           , 1                            , ...
    'AXISYM'                         , false                        , ...
    'P_SV'                           , true                         , ...
    'GPU_MODE'                       , false                        , ...
    'setup_with_binary_database'     , 0                            , ...
    'MODEL'                          , 'binary'                     , ...
    'SAVE_MODEL'                     , 'binary'                     , ...
    'ATTENUATION_VISCOELASTIC'       , false                        , ...
    'ATTENUATION_VISCOACOUSTIC'      , false                        , ...
    'N_SLS'                          , 3                            , ...
    'ATTENUATION_f0_REFERENCE'       , 5.1960                       , ...
    'READ_VELOCITIES_AT_f0'          , false                        , ...
    'USE_SOLVOPT'                    , false                        , ...
    'ATTENUATION_PORO_FLUID_PART'    , false                        , ...
    'Q0_poroelastic'                 , 1                            , ...
    'freq0_poroelastic'              , 10                           , ...
    'UNDO_ATTENUATION_AND_OR_PML'    , false                        , ...
    'NT_DUMP_ATTENUATION'            , 500                          , ...
    'NO_BACKWARD_RECONSTRUCTION'     , false                        , ...
    'NSOURCES'                       , 1                            , ...
    'force_normal_to_surface'        , false                        , ...
    'initialfield'                   , false                        , ...
    'add_Bielak_conditions_bottom'   , false                        , ...
    'add_Bielak_conditions_right'    , false                        , ...
    'add_Bielak_conditions_top'      , false                        , ...
    'add_Bielak_conditions_left'     , false                        , ...
    'ACOUSTIC_FORCING'               , false                        , ...
    'seismotype'                     , 1                            , ...
    'subsamp_seismos'                , 1                            , ...
    'USE_TRICK_FOR_BETTER_PRESSURE'  , false                        , ...
    'NSTEP_BETWEEN_OUTPUT_SEISMOS'   , 1000                         , ...
    'USER_T0'                        , 0                            , ...
    'save_ASCII_seismograms'         , true                         , ...
    'save_binary_seismograms_single' , true                         , ...
    'save_binary_seismograms_double' , false                        , ...
    'SU_FORMAT'                      , false                        , ...
    'use_existing_STATIONS'          , false                        , ...
    'nreceiversets'                  , 2                            , ...
    'anglerec'                       , 0                            , ...
    'rec_normal_to_surface'          , false                        , ...
    'RECEIVERS'                      , receiversets                 , ...
    'save_ASCII_kernels'             , true                         , ...
    'NSTEP_BETWEEN_COMPUTE_KERNELS'  , 1                            , ...
    'PML_BOUNDARY_CONDITIONS'        , true                         , ...
    'NELEM_PML_THICKNESS'            , 3                            , ...
    'ROTATE_PML_ACTIVATE'            , false                        , ...
    'ROTATE_PML_ANGLE'               , 30                           , ...
    'K_MIN_PML'                      , 1.0                          , ...
    'K_MAX_PML'                      , 1.0                          , ...
    'damping_change_factor_acoustic' , 0.5                          , ...
    'damping_change_factor_elastic'  , 1.0                          , ...
    'PML_PARAMETER_ADJUSTMENT'       , false                        , ...
    'STACEY_ABSORBING_CONDITIONS'    , false                        , ...
    'ADD_PERIODIC_CONDITIONS'        , false                        , ...
    'PERIODIC_HORIZ_DIST'            , 4000                         , ...
    'nbmodels'                       , 4                            , ...
    'MODELS'                         , models                       , ...
    'TOMOGRAPHY_FILE'                , './tomo_file.xyz'            , ...
    'read_external_mesh'             , false                        , ...
    'mesh_file'                      , './mesh_file'                , ...
    'nodes_coords_file'              , './nodes_coords_file'        , ...
    'materials_file'                 , './materials_file'           , ...
    'free_surface_file'              , './free_surface_file'        , ...
    'axial_elements_file'            , './axial_elements_file'      , ...
    'absorbing_surface_file'         , './absorbing_surface_file'   , ...
    'CPML_element_file'              , './CPML_element_file'        , ...
    'acoustic_forcing_surface_file'  , './acoustic_forcing_surface_file', ...
    'absorbing_cpml_file'            , './absorbing_cmpl_file'      , ...
    'tangential_detection_curve_file', './tangential_detection_curve_file', ...
    'interfacesfile'                 , './interfaces_file.dat'      , ...
    'xmin'                           , 0.0                          , ...
    'xmax'                           , 4000.0                       , ...
    'nx'                             , 80                           , ...
    'absorbbottom'                   , true                         , ...
    'absorbright'                    , true                         , ...
    'absorbtop'                      , false                        , ...
    'absorbleft'                     , true                         , ...
    'nbregions'                      , 4                            , ...
    'REGIONS'                        , regions                      , ...
    'NSTEP_BETWEEN_OUTPUT_INFO'      , 1000                         , ...
    'output_grid_Gnuplot'            , false                        , ...
    'output_grid_ASCII'              , false                        , ...
    'OUTPUT_ENERGY'                  , false                        , ...
    'NTSTEP_BETWEEN_OUTPUT_ENERGY'   , 200                          , ...
    'COMPUTE_INTEGRATED_ENERGY_FIELD', false                        , ...
    'NSTEP_BETWEEN_OUTPUT_IMAGES'    , 200                          , ...
    'cutsnaps'                       , 1.0                          , ...
    'output_color_image'             , true                         , ...
    'imagetype_JPEG'                 , 6                            , ...
    'factor_subsample_image'         , 1.0                          , ...
    'USE_CONSTANT_MAX_AMPLITUDE'     , false                        , ...
    'CONSTANT_MAX_AMPLITUDE_TO_USE'  , 1.17e4                       , ...
    'POWER_DISPLAY_COLOR'            , 3.00e-1                      , ...
    'DRAW_SOURCES_AND_RECEIVERS'     , true                         , ...
    'DRAW_WATER_IN_BLUE'             , true                         , ...
    'USE_SNAPSHOT_NUMBER_IN_FILENAME', false                        , ...
    'output_postscript_snapshot'     , true                         , ...
    'imagetype_postscript'           , 1                            , ...
    'meshvect'                       , true                         , ...
    'modelvect'                      , false                        , ...
    'boundvect'                      , true                         , ...
    'interpol'                       , true                         , ...
    'pointsdisp'                     , 6                            , ...
    'subsamp_postscript'             , 1                            , ...
    'sizemax_arrows'                 , 1.0                          , ...
    'US_LETTER'                      , false                        , ...
    'output_wavefield_dumps'         , false                        , ...
    'imagetype_wavefield_dumps'      , 1                            , ...
    'use_binary_for_wavefield_dumps' , false                        , ...
    'NUMBER_OF_SIMULTANEOUS_RUNS'    , 1                            , ...
    'BROADCAST_SAME_MESH_AND_MODEL'  , true                          ...
    );
end