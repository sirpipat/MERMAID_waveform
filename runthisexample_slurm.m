function runthisexample_slurm(masterdir)
% RUNTHISEXAMPLE_SLURM(masterdir)
%
% Executes RUNTHISEXAMPLE on all examples within the master directory. It 
% must be used with SLURM scheduler to function properly.
%
% INPUT:
% masterdir         full path to the directory containing all examples
%                   (subdirectories) for the simulation
%
% OUTPUT:
% no output
%
% SEE ALSO:
% RUNTHISEXAMPLE
%
% Last modified by sirawich-at-princeton.edu, 02/14/2022

% list of all examples
[listdirs, N] = allfile(masterdir);

% exit if the task ID is not assigned
if isempty(getenv("SLURM_ARRAY_TASK_ID"))
    fprintf('No valid SLURM_ARRAY_TASK_ID. Exit.\n');
    return
end

% get array task ID
idx = uint16(str2num(getenv("SLURM_ARRAY_TASK_ID")));

% figure out which example to run
ddir = listdirs(idx);
example = [ddir(1:end-2) '/'];

% execution
runthisexample(example, ddir, getenv("SPECFEM2D"));
end