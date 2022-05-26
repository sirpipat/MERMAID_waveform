function setimagenan(ax, im, c)
% SETIMAGENAN(ax, im, c)
%
% Sets NaN value in the imageplot to a target color. It modifies the image.
%
% INPUT:
% ax            target axes
% im            image object
% c             color for NaN [default: [1 1 1] (white)]
%
% Last modified by sirawich-at-princeton.edu, 05/26/2022

defval('c', [1 1 1])

% get the current colormap
cmap = colormap(ax);

% add the target color
cmap = [c; cmap];
colormap(ax, cmap);

% set NaN to the target color
minval = min(im.CData, [], 'all');
maxval = max(im.CData, [], 'all');
im.CData(isnan(im.CData)) = minval - 2 * (maxval - minval) / size(cmap, 1);
end