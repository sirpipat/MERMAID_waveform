function [bhan, than] = textbox(ax, x, y, str, varargintext, vararginpatch)
% [bhan, than] = TEXTBOX(ax, x, y, str, varargintext, vararginpatch)
%
% Puts a smart text box on a plot wherever you want to put the text. The
% text and the textbox are added to the axes. No new axes is created.
%
% INPUT:
% ax                axes where you want to add the text [default: gca]
% x                 x-position of the text
% y                 y-position of the text
% str               text string
% varargintext      keyword arguments for TEXT
% vararginpatch     keyword arguments for PATCH (aka the textbox)
%
% OUTPUT:
% bhan
% than
%
% SEE ALSO:
% TEXT, PATCH, BOXTEX, TEXTPATCH
%
% Last modified by sirawich-at-princeton.edu: 07/05/2023

defval('ax', gca)
defval('varargintext', [])
defval('vararginpatch', [])

% figure out text extent
if isempty(varargintext)
    tt = text(ax, x, y, str);
else
    tt = text(ax, x, y, str, varargintext{:});
end
te = get(tt, 'Extent');

% delete the original text object
delete(tt);

% margin of the patch
mp = 0.05;

% convert [left bottom width height] to (x,y) of the four corners
x_patch = te(1) + [-mp 1+mp 1+mp -mp] * te(3);
y_patch = te(2) + [-mp -mp 1+mp 1+mp] * te(4);

% create a textbox
if isempty(vararginpatch)
    bhan = patch(ax, x_patch, y_patch, [1 1 1]);
else
    bhan = patch(ax, x_patch, y_patch, [1 1 1], vararginpatch{:});
end

% put the text inside
if isempty(varargintext)
    than = text(ax, x, y, str);
else
    than = text(ax, x, y, str, varargintext{:});
end
end