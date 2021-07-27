function [ax, hs] = drawbackground(fname, ax)
% [ax, hs] = DRAWBACKGROUND(fname, ax)
%
% Draws the background image of the simulation, including the boundaries
% and layers.
%
% INPUT:
% fname         name of the interface file
% ax            axes to plot
%
% OUTPUT:
% ax            axes handle of the plot
% hs            polygon graphic objects
%
% Last modified by sirawich@princeton.edu, 07/26/2021

% list of layer colors
c = {[0.3 0.3 0.3], [0.45 0.45 0.45], [0.6 0.6 0.6], [0.75 0.75 0.75], ...
    [0.45 0.85 1]};
% read the interfacefile
[itfs, ~] = loadinterfacefile(fname);

hold on
for ii = 1:(length(itfs) - 1)
    x = [itfs{ii}.pts(:,1); itfs{ii+1}.pts(end:-1:1,1)];
    z = [itfs{ii}.pts(:,2); itfs{ii+1}.pts(end:-1:1,2)];
    pgon = polyshape(x, z);
    if ii < length(itfs) - 1
        cc = c{mod(ii-1, 4) + 1};
    else
        cc = c{end};
    end
    hs(ii) = plot(ax, pgon, 'FaceColor', cc, 'FaceAlpha', 1);
end

% set axis limit and data aspect ratio
set(ax, 'XLim', [min(itfs{1}.pts(:,1)) max(itfs{1}.pts(:,1))], ...
    'YLim', [min(itfs{1}.pts(:,2)) max(itfs{end}.pts(:,2))], ...
    'DataAspectRatio', [1 1 1]);
end