function c = bwrmap(n, order)
% Creates blue-white-red colormap like the one in matplotlib.
%
% INPUT:
% n         number of colors    [default: 255]
% order     'bwr'               [default]
%           'rwb'
%
% OUTPUT:
% c         n x 3 array of RGB triplets
%
% EXAMPLE:
% >> bwrmap(3)
%
% ans = 
%
%       0   0   1       % blue
%       1   1   1       % white
%       1   0   0       % red
%
% Last modified by sirawich-at-princeton.edu, 02/25/2022

defval('n', 255)
defval('order', 'bwr')

% column vector of incrementing numbers starting from zero
nn = (0:(n-1))';

% scale in the denominator
c_factor = (n-1)/2;

c = zeros(n, 3);
if strcmpi(order, 'rwb')
    c(:,1) = min(2 - nn / c_factor, 1);
    c(:,2) = max(1 - abs(nn / c_factor - 1), 0);
    c(:,3) = min(nn / c_factor, 1);
else
    c(:,1) = min(nn / c_factor, 1);
    c(:,2) = max(1 - abs(nn / c_factor - 1), 0);
    c(:,3) = min(2 - nn / c_factor, 1);
end
end

%    |  RED     |   GREEN   |    BLUE
% 1_ |_____     |     _     |     _____
%    |     \    |    / \    |    /
%    |      \   |   /   \   |   /
%    |       \  |  /     \  |  /
% 0- |        \ | /       \ | /