function g = diff_central(f, h)
% Compute the derivative of f(x) using first-order central finite 
% difference method on a uniform mesh.
%
% INPUT
% f         f(x)
% h         mesh spacing
% OUTPUT
% g         f'(x)
%
% Last modified by Sirawich Pipatprathanporn, 02/17/2021

g = (circshift(f,-1) - circshift(f,1)) / 2 / h;
g(1) = (f(2) - f(1)) / h;
g(end) = (f(end) - f(end-1)) / h;
end