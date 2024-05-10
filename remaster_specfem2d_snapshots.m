function [imgout, ax, axb] = remaster_specfem2d_snapshots(imgin, interfacefile, method, plt, labeltext)
% [imgout, ax, axb] = REMASTER_SPECFEM2D_SNAPSHOTS(imgin, interfacefile, method, plt, labeltext)
%
% Redraws SPECFEM2D JPEG wavefield snapshots such that the wavefield is
% more readable and save to a JPEC file at the same file directory. See 
% 'method' variable description for how the redrew wavefield is more 
% readable. It can also make a figure and print to $EPS.
%
% INPUT:
% imgin             full filename to a SPECFEM2D JPEG snapshot
% interfacefile     full filename to the interface file for the simulation
% method            method to make a wavefield more readable [default: 0]
%       0 - Set red-blue color scale to KELICOL colormap and the original
%           color intensity dictates transparency
%       1 - Set wavefield transparency to the original color ntensity
%       2 - Make the wavefield a bit brighter but still not readable
%           (This is the old method, not recommended)
% plt               whether to make a figure plot or not [default: true]
% labeltext         text for plot label (e.g. 'a', 'b', 'c', 'd') for
%                   making paper's figures. If this is empty, no boxed
%                   label will be created. Note that it will not affect the
%                   redrew SPECFEM2D snapshot. [default: []]
%
% OUTPUT:
% imgout            full filename to a redred SPECFEM2D snapshot
% ax                axes handle of the snapshot
% axb               axes handle of the boxed label
%
% SEE ALSO:
% BOXEDLABEL, KELICOL
% 
% Last modified by sirawich-at-princeton.edu, 05/09/2024

defval('method', 0)
defval('plt', true)
defval('labeltext', [])

% defined color
COLOR_WATER = [135 206 250];
COLOR_CRUST = [115 115 115];
COLOR_RED = [255 0 0];
COLOR_BLUE = [0 0 255];

% defined colormap
KELICOL = kelicol * 255;
KELICOL_RED = KELICOL(1:81,:);
KELICOL_RED_SCALE = (255:(-255/80):0)';
KELICOL_BLUE = KELICOL(82:162,:);
KELICOL_BLUE_SCALE = (0:(255/80):255)';

%% read image to a [M N P] matrix where
% M = index of Z direction
% N = index of X direction
% P = index of RGB color (1 = RED, 2 = GREEN, 3 = BLUE)
im = imread(imgin);
SIZE_Z = size(im, 1);
SIZE_X = size(im, 2);

% flatten the image to a list of RGB tripets for all pixels
imflat = reshape(im, size(im, 1) * size(im, 2), size(im, 3));

%% read interfacefile
itfs = loadinterfacefile(interfacefile);

% dimension of a physical domain
width = itfs{3}.pts(2, 1);
height = itfs{3}.pts(2, 2);

% create a mesh grid of the physical location for each image pixel
x = (1:size(im, 2)) / size(im,2) * width;
z = (403 - (1:size(im, 1))) / (size(im,1) - 20) * height;
[xx, zz] = meshgrid(x, z);

% interpolated ocean bottom
xob = x;
zob = interp1(itfs{2}.pts(:,1), itfs{2}.pts(:,2), xob, 'linear');

% determine the background type (crust, water air)
is_air = (zz > 9600);
is_crust = (zz < zob);
is_water = and(not(is_air), not(is_crust));

%% determine regions:
% You may try plotting (R,G,B) for all pixels (imflat would work best) to 
% understand how I define the regions as below.

% Positive (blue) and negative (red) regions on unflattened image
id_p = and(im(:,:,2) < 50, im(:,:,1) <= im(:,:,3));
id_n = and(im(:,:,2) < 50, im(:,:,1) > im(:,:,3));

% Positive (blue) and negative (red) regions on   flattened image
% 1: positive perturbation (GREEN < 50 and RED <= BLUE)
id_1 = and(imflat(:,2) < 50, imflat(:,1) <= imflat(:,3));

% 2: negative perturbation (GREEN < 50 and RED > BLUE)
id_2 = and(imflat(:,2) < 50, imflat(:,1) > imflat(:,3));

% 3: air (RED, GREEN, BLUE > 235)
id_3 = and(and(imflat(:,1) > 235, imflat(:,2) > 235),imflat(:,3) > 235);

% 4: water (GREEN > 225 - RED/2, GREEN <= 450 - RED, and BLUE > 185)
id_4 = and(and(imflat(:,2) > 225 - imflat(:,1)/2, imflat(:,2) <= 450 - imflat(:,1)), imflat(:,3) > 185);

% 5: source (GREEN > 225 - RED/2, GREEN <= 450 - RED, and BLUE <= 185)
id_5 = and(and(imflat(:,2) > 225 - imflat(:,1)/2, imflat(:,2) <= 450 - imflat(:,1)), imflat(:,3) <= 185);

% 6: receiver (GREEN >= 50, GREEN < 136, and GREEN <= 225 - RED/2)
id_6 = and(and(imflat(:,2) >= 50, imflat(:,2) < 136), imflat(:,2) <= 225 - imflat(:,1)/2);

% 7: crust (GREEN >= 136 and GREEN <= 225 - RED/2)
id_7 = and(imflat(:,2) >= 136, imflat(:,2) <= 225 - imflat(:,1)/2);

%% rescale the extrema of the waveform to 255
% identify max BLUE in positive perturbation region and scale
max_blue = double(max(imflat(id_1, 3)));

% identify max RED in negative perturbation region and scale
max_red = double(max(imflat(id_2, 1)));

% rescaling
% non-linear scaling (power rule) when n ~= 1
% linear scaling when n == 1
n = 1;
% ignore recoloring if ther is no red color to rescale
if ~isempty(max_red) && max_red > 0
    rescaled_red = 255 * min((double(im(:,:,1)) / max_red) .^ n, 1);
else
    rescaled_red = double(im(:,:,1));
end
% ignore recoloring if ther is no blue color to rescale
if ~isempty(max_blue) && max_blue > 0
    rescaled_blue = 255 * min((double(im(:,:,3)) / max_blue) .^ n, 1);
else
    rescaled_blue = double(im(:,:,3));
end

%% redraw the wavefield image
% For method 0 and 1, the calculation is as follow
% imnew = im(non-wavefield) + im_modified(wavefield)
% where
% im_modified(wavefield) = im(wavefield) .* ...
%   (wavefield_color .* transparency + BG_color .* (1-transparency))
% where
% transparency = color_intensity / 255.
%
% The complexity of the expression comes from the need to broadcast every
% objects to [SIZE_Z SIZE_X 3] arrays. Simply calling REPMAT is not enough
% for some cases.
if method == 0
    % latest/default method - transpraency with kelikol color scheme
    imnew = uint8(double(repmat(not(or(id_p, id_n)), [1 1 3])) .* double(im) + ...
        double(repmat(id_p, [1 1 3])) .* interp1(KELICOL_BLUE_SCALE, KELICOL_BLUE, rescaled_blue, 'linear') .* repmat(rescaled_blue, [1 1 3]) / 255 + ...
        double(repmat(and(id_p, is_crust), [1 1 3])) .* broadcast(COLOR_CRUST, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_blue, [1 1 3]) / 255) + ...
        double(repmat(and(id_p, is_water), [1 1 3])) .* broadcast(COLOR_WATER, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_blue, [1 1 3]) / 255) + ...
        double(repmat(id_n, [1 1 3])) .* interp1(KELICOL_RED_SCALE, KELICOL_RED, rescaled_red, 'linear') .* repmat(rescaled_red, [1 1 3]) / 255 + ...
        double(repmat(and(id_n, is_crust), [1 1 3])) .* broadcast(COLOR_CRUST, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_red, [1 1 3]) / 255) + ...
        double(repmat(and(id_n, is_water), [1 1 3])) .* broadcast(COLOR_WATER, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_red, [1 1 3]) / 255));
elseif method == 1
    % old method - transparency with red/blue color scheme
    imnew = uint8(double(repmat(not(or(id_p, id_n)), [1 1 3])) .* double(im) + ...
        double(repmat(id_p, [1 1 3])) .* broadcast(COLOR_BLUE, [SIZE_Z SIZE_X]) .* double(repmat(im(:,:,3), [1 1 3])) / 255 + ...
        double(repmat(and(id_p, is_crust), [1 1 3])) .* broadcast(COLOR_CRUST, [SIZE_Z SIZE_X]) .* (1 - double(repmat(rescaled_blue, [1 1 3])) / 255) + ...
        double(repmat(and(id_p, is_water), [1 1 3])) .* broadcast(COLOR_WATER, [SIZE_Z SIZE_X]) .* (1 - double(repmat(rescaled_blue, [1 1 3])) / 255) + ...
        double(repmat(id_n, [1 1 3])) .* broadcast(COLOR_RED, [SIZE_Z SIZE_X]) .* double(repmat(im(:,:,1), [1 1 3])) / 255 + ...
        double(repmat(and(id_n, is_crust), [1 1 3])) .* broadcast(COLOR_CRUST, [SIZE_Z SIZE_X]) .* (1 - double(repmat(rescaled_red, [1 1 3])) / 255) + ...
        double(repmat(and(id_n, is_water), [1 1 3])) .* broadcast(COLOR_WATER, [SIZE_Z SIZE_X]) .* (1 - double(repmat(rescaled_red, [1 1 3])) / 255));
elseif method == 2
    % old method - make red/blue brighter but color gradient is bad
    imflat(id_1, 3) = uint8((255^3 * double(imflat(id_1, 3))) .^ 0.25);
    imflat(id_2, 1) = uint8((255^3 * double(imflat(id_2, 1))) .^ 0.25);
    imnew = reshape(imflat, 403, 800, 3);
else
    % DEFAULT METHOD (same as method == 0)
    % transpraency
    imnew = uint8(double(repmat(not(or(id_p, id_n)), [1 1 3])) .* double(im) + ...
        double(repmat(id_p, [1 1 3])) .* interp1(KELICOL_BLUE_SCALE, KELICOL_BLUE, rescaled_blue, 'linear') .* repmat(rescaled_blue, [1 1 3]) / 255 + ...
        double(repmat(and(id_p, is_crust), [1 1 3])) .* broadcast(COLOR_CRUST, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_blue, [1 1 3]) / 255) + ...
        double(repmat(and(id_p, is_water), [1 1 3])) .* broadcast(COLOR_WATER, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_blue, [1 1 3]) / 255) + ...
        double(repmat(id_n, [1 1 3])) .* interp1(KELICOL_RED_SCALE, KELICOL_RED, rescaled_red, 'linear') .* repmat(rescaled_red, [1 1 3]) / 255 + ...
        double(repmat(and(id_n, is_crust), [1 1 3])) .* broadcast(COLOR_CRUST, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_red, [1 1 3]) / 255) + ...
        double(repmat(and(id_n, is_water), [1 1 3])) .* broadcast(COLOR_WATER, [SIZE_Z SIZE_X]) .* (1 - repmat(rescaled_red, [1 1 3]) / 255));
end

if plt
    figure(1)
    clf
    set(gcf, 'Units', 'pixels', 'Position', [1 200 SIZE_X+100 SIZE_Z+150]);
    ax = gca;
    set(ax, 'Units', 'pixels', 'Position', [50 50 SIZE_X SIZE_Z]);
    % draw the remastered image
    imagesc(imnew);
    set(ax, 'Units', 'normalized', 'Box', 'on', 'XTick', [], ...
        'YTick', [], 'YLim', [19.5 SIZE_Z+0.5]);
    nolabels(ax);
    % draw the ocean bottom
    hold on
    plot(itfs{2}.pts(:,1) / width * SIZE_X, ...
        (height - itfs{2}.pts(:,2)) / height * (SIZE_Z - 20) + 20, ...
        'LineWidth', 0.5, 'Color', 'k')
    if ~isempty(labeltext)
        % add boxedlabel
        axb = boxedlabel(ax, 'northwest', 50, 'pixels', labeltext, ...
            'FontSize', 40);
    end
    set(gcf, 'Renderer', 'painters')
    figname = sprintf('%s_%s.eps', mfilename, ...
        cindeks(split(removepath(imgin), '.'), 1));
    figdisp(figname, [], [], 2, [], 'epstopdf')
else
    ax = [];
    axb = [];
end
words = split(imgin, '.');
imgout = [words{1} '_remastered.jpg'];
imwrite(imnew, imgout, 'jpg')
end

% create a 3-dimensional array of size [M N P] by duplicating an input 
% vector V with length P over M x N times
% INPUT:
% v         a vector of length P
% [M N]     how many times in each direction to duplicate V
%
% OUTPUT:
% a         [M N P] array of duplicated vector V
function a = broadcast(v, dim2)
if size(v, 1) > 1
    v = v';
end
a = permute(repmat(v, [dim2(1) 1 dim2(2)]), [1 3 2]);
end