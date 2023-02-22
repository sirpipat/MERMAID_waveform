function [Mdc, Miso, Mclvd, pDC] = decomposemoment(M)
% [Mdc, Miso, Mclvd, pDC] = decomposemoment(M)
%
% Computes the double-couple, isotropic, and CLVD elements of a moment
% tensor.
%
% INPUT:
% M             moment tensor (3x3 matrix)
%
% OUTPUT:
% Mdc           double-couple element of a moment tensor
% Miso          isotropic element of a moment tensor
% Mclvd         CLVD element of a moment tensor
% pDC           percentage of the double-couple element
%
% Last modified by sirawich-at-princeton.edu: 02/22/2023

Miso = (trace(M) / 3) * eye(3);

% compute the principal axes
[V, D] = eig(M);

% remove the isotropic element
D_deviatoric = D - (trace(D) / 3) * eye(3);

% sort the eigenvalues by their magnitudes
[~, id] = sort(abs(diag(D_deviatoric)), 'descend');
D_mag = D_deviatoric(id, id);
V_mag = V(:, id);

% determine the double couple percentage
pDC = (1 - 2 * abs(D_mag(3,3) / D_mag(1,1))) * 100;

% diagonalized moment tensor of the double couple element
Ddc = diag([D_mag(1,1)+D_mag(3,3)/2, -D_mag(3,3)/2-D_mag(1,1), 0]);

Mdc = V_mag * Ddc * V_mag';
Mclvd = M - Miso - Mdc;
end