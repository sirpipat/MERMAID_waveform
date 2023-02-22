function [t, p, b, M0] = moment2axes(M)
% [t, p, b, M0] = moment2axes(M)
%
% Computes the principal axes of the moment tensor.
%
% INPUT:
% M             moment tensor (3x3 matrix)
%
% OUTPUT:
% t             tension axis
% p             pressure axis
% b             null axis
% M0            scalar moment
%
% SEE ALSO:
% MOMENT2SLIP
%
% Last modified by sirawich-at-princeton.edu: 02/22/2023

%% Part I : determine the scalar magnitude and then normalize
M0 = sqrt(sum(M.^2, 'all') / 2);
M = M / M0;

%% Part II : determine the principal axes
% compute the principal axes
[V, D] = eig(M);

% sort the eigenvalues
[~, it] = sort(diag(D), 'descend');
V = V(:, it);

% sort into axes
t = V(:, 1);
b = V(:, 2);
p = V(:, 3);
end