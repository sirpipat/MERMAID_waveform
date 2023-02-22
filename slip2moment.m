function M = slip2moment(slip, normal, M0)
% M = SLIP2MOMENT(slip, normal, M0)
%
% Construct a 3x3 matrix of the moment tensor given a slip on a fault
% planes and the scalar moment
%
% INPUT:
% slip          the slip vector
% normal        the normal vector to the fault plane
% M0            scalar moment
%
% OUTPUT:
% M             moment tensor
%
% SEE ALSO:
% MOMENT2SLIP, MOMENT2AXES
%
% Last modified by sirawich-at-princeton.edu: 02/22/2023

M = M0 * (slip * normal' + normal * slip');
end