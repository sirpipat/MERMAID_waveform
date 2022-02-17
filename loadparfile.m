function params = loadparfile(fname)
% params = LOADPARFILE(fname)
%
% Reads parameters from a Par_file.
%
% DISCLAIMER: This is not the official way to read/write Par_file. I just
% go through comments and parameters in an instant of Par_file and
% read/write accordingly.
%
% DISCLAMIER 2: It is designed to read Par_file properly only in the two
% versions (commits) of SPECFEM2D:
%      'master' (commit: e937ac2f74f23622f6ebbc8901d30fb33c1a2c38)
%       'devel' (commit: cf89366717d9435985ba852ef1d41a10cee97884)
% The fields in params refers to the parameter names in 'master' version
% except those exist only in 'devel' version.
%
% INPUT:
% fname         name of the Par_file
%
% OUTPUT:
% params        parameters
%
% SEE ALSO:
% WRITEPARFILE, MAKEPARAMS
%
% Last modified by Sirawich Pipatprathanporn, 02/17/2022

names = {};
values = {};
numvar = 0;

fid = fopen(fname, 'r');

line = fgetl(fid);
while ischar(line)
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    numvar = numvar + 1;
    [name, value] = readgeneric(line);
    names{numvar} = name;
    values{numvar} = value;
    
    % handle sets of reveivers
    if strcmp(name, 'nreceiversets')
        nreceivers = value;
    elseif strcmp(name, 'rec_normal_to_surface')
        receivers = cell(1, nreceivers);
        % read sets of receivers
        for ii = 1:nreceivers
            line = fgetl(fid);
            % skip the comment
            while isempty(line) || strcmp(line(1), '#')
                line = fgetl(fid);
            end
            if ~ischar(line)
                break;
            end
            [~, receiver.nrec] = readgeneric(line);
            line = fgetl(fid);
            [~, receiver.xdeb] = readgeneric(line);
            line = fgetl(fid);
            [~, receiver.zdeb] = readgeneric(line);
            line = fgetl(fid);
            [~, receiver.xfin] = readgeneric(line);
            line = fgetl(fid);
            [~, receiver.zfin] = readgeneric(line);
            line = fgetl(fid);
            [~, receiver.record_at_surface_same_vertical] = readgeneric(line);
            receivers{ii} = receiver;
        end
        numvar = numvar + 1;
        names{numvar} = 'RECEIVERS';
        values{numvar} = receivers;
    % handle reading the velocity and density models
    elseif strcmp(name, 'nbmodels')
        nbmodels = value;
        % skip the comment (instruction)
        line = fgetl(fid);
        while isempty(line) || strcmp(line(1), '#')
            line = fgetl(fid);
        end
        if ~ischar(line)
            break;
        end
        % stores each parameter for a model (stored as a struct, see below)
        models = cell(1, nbmodels);
        % read models
        for ii = 1:nbmodels
            nums = sscanf(replace(line, 'd', 'e'), '%f');
            model.model_number = nums(1);
            model.type_number = nums(2);
            if nums(2) == 1
                if nums(5) > 0
                    model.type_name = 'elastic';
                else
                    model.type_name = 'acoustic';
                end
                model.rho = nums(3);
                model.vp = nums(4);
                model.vs = nums(5);
                model.QKappa = nums(8);
                model.Qmu = nums(9);
            elseif nums(2) == 2
                model.type_name = 'anisotropic';
                model.rho = nums(3);
                model.c11 = nums(4);
                model.c13 = nums(5);
                model.c15 = nums(6);
                model.c33 = nums(7);
                model.c35 = nums(8);
                model.c55 = nums(9);
                model.c12 = nums(10);
                model.c23 = nums(11);
                model.c25 = nums(12);
                model.c22 = nums(13);
            elseif nums(2) == 3
                model.type_name = 'poroelastic';
                model.rhos = nums(3);
                model.rhof = nums(4);
                model.phi = nums(5);
                model.c = nums(6);
                model.kxx = nums(7);
                model.kxz = nums(8);
                model.kzz = nums(9);
                model.Ks = nums(10);
                model.Kf = nums(11);
                model.Kfr = nums(12);
                model.etaf = nums(13);
                model.mufr = nums(14);
                model.Qu = nums(15);
            elseif nums(2) == -1
                model.A = nums(5);
            else
                model.type_name = 'invalid';
                fprintf('Invaid velocity and density model.\n')
            end
            models{ii} = model;
            
            % read the next model if the current is not the last one.
            if ii < nbmodels
                line = fgetl(fid);
            end
        end
        numvar = numvar + 1;
        names{numvar} = 'MODELS';
        values{numvar} = models;
    % handle reading the velocity and density models
    elseif strcmp(name, 'nbregions')
        nbregions = value;
        % stores each parameter for a model (stored as a struct, see below)
        regions = cell(1, nbregions);
        
        % skip the comment (instruction)
        line = fgetl(fid);
        while isempty(line) || strcmp(line(1), '#')
            line = fgetl(fid);
        end
        if ~ischar(line)
            break;
        end
        
        % read regions
        for ii = 1:nbregions
            nums = sscanf(replace(line, 'd', 'e'), '%d');
            region.nxmin = nums(1);
            region.nxmax = nums(2);
            region.nzmin = nums(3);
            region.nzmax = nums(4);
            region.material_number = nums(5);
            regions{ii} = region;
            % read the next model if the current is not the last one.
            if ii < nbregions
                line = fgetl(fid);
            end
        end
        numvar = numvar + 1;
        names{numvar} = 'REGIONS';
        values{numvar} = regions;
    end
    
    % read the next line
    line = fgetl(fid);
end

fclose(fid);

params = cell2struct(values, names, 2);
end

function [name, value] = readgeneric(line)
% list of variables with known data type
bool_var   = {'SAVE_FORWARD', ...
              'AXISYM', ...
              'P_SV', ...
              'GPU_MODE', ...
              'ATTENUATION_VISCOELASTIC', ...
              'ATTENUATION_VISCOACOUSTIC', ...
              'READ_VELOCITIES_AT_f0', ...
              'USE_SOLVOPT', ...
              'ATTENUATION_PORO_FLUID_PART', ...
              'UNDO_ATTENUATION_AND_OR_PML', ...
              'NO_BACKWARD_RECONSTRUCTION', ...
              'force_normal_to_surface', ...
              'initialfield', ...
              'add_Bielak_conditions_bottom', ...
              'add_Bielak_conditions_right', ...
              'add_Bielak_conditions_top', ...
              'add_Bielak_conditions_left', ...
              'ACOUSTIC_FORCING', ...
              'write_moving_sources_database', ...
              'USE_TRICK_FOR_BETTER_PRESSURE', ...
              'save_ASCII_seismograms', ...
              'save_binary_seismograms_single', ...
              'save_binary_seismograms_double', ...
              'SU_FORMAT', ...
              'use_existing_STATIONS', ...
              'rec_normal_to_surface', ...
              'record_at_surface_same_vertical', ...
              'save_ASCII_kernels', ...
              'APPROXIMATE_HESS_KL', ...
              'PML_BOUNDARY_CONDITIONS', ...
              'ROTATE_PML_ACTIVATE', ...
              'PML_PARAMETER_ADJUSTMENT', ...
              'STACEY_ABSORBING_CONDITIONS', ...
              'ADD_PERIODIC_CONDITIONS', ...
              'read_external_mesh', ...
              'absorbbottom', ...
              'absorbright', ...
              'absorbtop', ...
              'absorbleft', ...
              'output_grid_Gnuplot', ...
              'output_grid_ASCII', ...
              'OUTPUT_ENERGY', ...
              'COMPUTE_INTEGRATED_ENERGY_FIELD', ...
              'output_color_image', ...              
              'USE_CONSTANT_MAX_AMPLITUDE', ...
              'DRAW_SOURCES_AND_RECEIVERS', ...
              'DRAW_WATER_IN_BLUE', ...
              'USE_SNAPSHOT_NUMBER_IN_FILENAME', ...
              'output_postscript_snapshot', ...
              'meshvect', ...
              'modelvect', ...
              'boundvect', ...
              'interpol', ...
              'US_LETTER', ...
              'output_wavefield_dumps', ...
              'use_binary_for_wavefield_dumps', ...
              'BROADCAST_SAME_MESH_AND_MODEL'
              };
float_var  = {'DT', ...
              'ATTENUATION_f0_REFERENCE', ...
              'Q0_poroelastic', ...
              'freq0_poroelastic', ...
              'USER_T0', ...
              'anglerec', ...
              'xdeb', ...
              'zdeb', ...
              'xfin', ...
              'zfin', ...
              'ROTATE_PML_ANGLE', ...
              'K_MIN_PML', ...
              'K_MAX_PML', ...
              'damping_change_factor_acoustic', ...
              'damping_change_factor_elastic', ...
              'PERIODIC_HORIZ_DIST', ...
              'xmin', ...
              'xmax', ...
              'cutsnaps', ...
              'factor_subsample_image', ...
              'CONSTANT_MAX_AMPLITUDE_TO_USE', ...
              'POWER_DISPLAY_COLOR', ...
              'sizemax_arrows'
              };
int_var    = {'SIMULATION_TYPE', ...
              'NOISE_TOMOGRAPHY', ...
              'NPROC', ...
              'partitioning_method', ...
              'ngnod', ...
              'NSTEP', ...
              'time_stepping_scheme', ...
              'setup_with_binary_database', ...
              'N_SLS', ...
              'NT_DUMP_ATTENUATION', ...
              'NSOURCES', ...
              'noise_source_time_function_type', ...
              'seismotype', ...
              'subsamp_seismos', ...
              'NSTEP_BETWEEN_OUTPUT_SEISMOS', ...
              'nreceiversets', ...
              'nrec', ...
              'NSTEP_BETWEEN_COMPUTE_KERNELS', ...
              'NELEM_PML_THICKNESS', ...
              'nbmodels', ...
              'nx', ...
              'nbregions', ...
              'NSTEP_BETWEEN_OUTPUT_INFO', ...
              'NTSTEP_BETWEEN_OUTPUT_ENERGY', ...
              'NSTEP_BETWEEN_OUTPUT_IMAGES', ...
              'imagetype_JPEG', ...
              'imagetype_postscript', ...
              'pointsdisp', ...
              'subsamp_postscript', ...
              'imagetype_wavefield_dumps', ...
              'NUMBER_OF_SIMULTANEOUS_RUNS'
              };
string_var = {'title', ...
              'MODEL', ...
              'SAVE_MODEL', ...
              'TOMOGRAPHY_FILE', ...
              'mesh_file', ...
              'nodes_coords_file', ...
              'materials_file', ...
              'free_surface_file', ...
              'axial_elements_file', ...
              'absorbing_surface_file', ...
              'CPML_element_file', ...
              'acoustic_forcing_surface_file', ...
              'absorbing_cpml_file', ...
              'tangential_detection_curve_file', ...
              'interfacesfile'
              }; 

name = sscanf(line, '%s', 1);

% Handle the different variable naming in 'master' and 'devel' branches.
% If both branches have two different variable names referring to the same
% thing, the name in 'master' branch is used here.
if strcmp(name, 'PARTITIONING_TYPE')
    name = 'partitioning_method';
end
if strcmp(name, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS')
    name = 'NSTEP_BETWEEN_OUTPUT_SEISMOS';
end
if strcmp(name, 'NTSTEP_BETWEEN_COMPUTE_KERNELS')
    name = 'NSTEP_BETWEEN_COMPUTE_KERNELS';
end
if strcmp(name, 'NTSTEP_BETWEEN_OUTPUT_SAMPLE')
    name = 'subsamp_seismos';
end
if strcmp(name, 'NTSTEP_BETWEEN_OUTPUT_IMAGES')
    name = 'NSTEP_BETWEEN_OUTPUT_IMAGES';
end
if strcmp(name, 'NTSTEP_BETWEEN_OUTPUT_INFO')
    name = 'NSTEP_BETWEEN_OUTPUT_INFO';
end
if strcmp(name, 'NGNOD')
    name = 'ngnod';
end

if any(strcmp(bool_var, name))
    value = readbool(line);
elseif any(strcmp(int_var, name))
    value = readint(line);
elseif any(strcmp(float_var, name))
    value = readfloat(line);
elseif any(strcmp(string_var, name))
    value = readstring(line);
else
    fprintf('unable to determine type of "%s" variable. Read as string\n', ...
        name);
    keyboard;
    value = readstring(line);
end
end

function value = readstring(line)
% find equal sign
where_start = strfind(line, '=');

% find # where the comment starts
where_end = strfind(line, '#');
if isempty(where_end)
    where_end = length(line) + 1;
end

% read the value
value = strip(sscanf(line((where_start+1):(where_end-1)), '%c'));
end

function value = readbool(line)
value = readstring(line);
if strcmp(value, '.true.')
    value = true;
elseif strcmp(value, '.false.')
    value = false;
else
    % do not know what to do
    error(strcat('ValueError: cannot read a boolean\n', ...
                 sprintf('line >> %s\n', line)));
end
end

function value = readint(line)
% find equal sign
where = strfind(line, '=');
% read the value
value = sscanf(line((where+1):end), '%d', 1);
end

function value =  readfloat(line)
% find equal sign
where = strfind(line, '=');
% change the exponent notation syntax from 'd' to 'e'
line = replace(line, 'd', 'e');
% read the value
value = sscanf(line((where+1):end), '%f', 1);
end