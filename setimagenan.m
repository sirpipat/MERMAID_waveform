function setimagenan(ax, im, c, minval, maxval)
% SETIMAGENAN(ax, im, c, minval, maxval)
%
% Sets NaN value in the imageplot to a target color. It modifies the image.
%
% INPUT:
% ax            target axes
% im            image object
% c             color for NaN [default: [1 1 1] (white)]
% minval        lower limit of color axis 
%               [default: min(im.CData, [], 'all')]
% maxval        upper limit of color axis
%               [default: max(im.CData, [], 'all')]
%
% Last modified by sirawich-at-princeton.edu, 06/07/2022

defval('c', [1 1 1])

% get the current colormap
cmap = colormap(ax);

% add the target color
cmap = [c; cmap];
colormap(ax, cmap);

% set NaN to the target color
defval('minval', min(min(im.CData, [], 'omitnan'), [], 'omitnan'))
defval('maxval', max(max(im.CData, [], 'omitnan'), [], 'omitnan'))
im.CData(isnan(im.CData)) = minval - 2 * (maxval - minval) / (size(cmap, 1) - 1);
end