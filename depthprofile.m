function pf = depthprofile(x, itfs, mt)
% pf = DEPTHPROFILE(x, itfs, mt)
%
% Calculates the depth profile at any position for given interfaces.
%
% INPUT:
% x             x-position of the profile
% itfs          array of interface structs
% mt            interpolation method
%
% OUTPUT:
% pf            depth profile
%
% SEE ALSO:
% MAKEPARAMS, WRITEINTERFACEFILE
%
% Last modified by sirawich@princeton.edu, 08/08/2021

n = length(itfs);
pf = zeros(n, 1);
for ii = 1:n
    pf(ii) = interp1(itfs{ii}.pts(:,1), itfs{ii}.pts(:,2), x, mt);
end
end