function animatepropagation(ddir, savedir)
% ANIMATEPROPAGATION(ddir, savedir)
%
% Takes simulation's snapshots and turns to a video.
%
% INPUT:
% ddir          simulation directory
% savedir       directory of the saved video
%               The video is saved as 'images2video_savename.mp4'
%
% Last modified by sirawich@princeton.edu, 07/25/2021

% identifies all snapshots
snapshots = ls2cell([ddir 'OUTPUT_FILES/forward_image*'], 1);

p = loadparfile([ddir 'DATA/Par_file']);

fr = 1 / p.NSTEP_BETWEEN_OUTPUT_IMAGES / p.DT;

savename = strcat(mfilename, '_', removepath(ddir(1:end-1)));
images2video(snapshots, 'jpg', fr, savename, savedir);
end