function [slip, normal, M0] = moment2slip(M)
% [slip, normal, M0] = MOMENT2SLIP(M)
%
% Computes the slip vector and normal vector to the fault plane from a
% moment tensor. Note that the slip vector and normal vector can be swapped
% depending on which of the two nodal planes is actually the fault plane.
%
% see Jost and Herrmann, A Student's Guide to and Review of Moment Tensors
% (1989) Equation 11-14
%
% INPUT:
% M             moment tensor (3x3 matrix)
%
% OUTPUT:
% slip          the slip vector
% normal        the normal vector to the fault plane
% M0            scalar moment
%
% SEE ALSO:
% SLIP2MOMENT, MOMENT2AXES
%
% Last modified by sirawich-at-princeton.edu: 03/01/2023

[t, p, ~, M0] = moment2axes(M);

% Jost and Herrmann, A Student's Guide to and Review of Moment Tensors
% (1989) Equation 11-14
slip = (t - p) / sqrt(2);
normal = (t + p) / sqrt(2);
end