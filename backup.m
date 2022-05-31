function backup
% BACKUP
%
% Save all variables in the working space to an autosaved file. Call it
% regularly to prevent loss due to a MATLAB crash or a machine crash.
%
% The backup file is saved to $MFILES/backup/
%
% It keeps up to 3 back up files with 'backup_01.mat' as the most recent 
% backup file and 'backup_03.mat' as the most ancient backup file. Note
% that if 'backup_03.mat' already exists, it will be overwritten.
%
% SEE ALSO:
% LOADBACKUP
%
% Last modified by sirawich-at-princeton.edu, 05/31/2022

% backup file name
fname = sprintf('%sbackup/backup_01.mat', getenv('MFILES'));

% Pushed the backup file number if it is already exist.
% The newly saved file is saved as backup_01.mat nonetheless, but the
% existing backup_01.mat is changed to backup_02.mat, and the
% existing backup_02.mat is changed to backup_03.mat. The
% existing backup_03.mat is overwritten.
if exist(fname, 'file')
    fname2 = sprintf('%sbackup/backup_02.mat', getenv('MFILES'));
    if exist(fname2, 'file')
        fname3 = sprintf('%sbackup/backup_03.mat', getenv('MFILES'));
        system(sprintf('mv %s %s', fname2, fname3));
    end
    system(sprintf('mv %s %s', fname, fname2));
end

% save the variable into the backup file
evalin('base', sprintf('save(''%s'')', fname));

end

