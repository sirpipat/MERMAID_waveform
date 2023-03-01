function M = slip2moment(slip, normal, M0)
% M = SLIP2MOMENT(slip, normal, M0)
%
% Construct a 3x3 matrix of the moment tensor given a slip on a fault
% planes and the scalar moment
%
% see Jost and Herrmann, A Student's Guide to and Review of Moment Tensors
% (1989) Equation 6
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
% Last modified by sirawich-at-princeton.edu: 03/01/2023

% Jost and Herrmann, A Student's Guide to and Review of Moment Tensors
% (1989) Equation 6
M = M0 * (slip * normal' + normal * slip');
end