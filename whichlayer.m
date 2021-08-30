function layer = whichlayer(pt, itfs, mt)
% material = WHICHLAYER([x z], itfs, mt)
%
% Determines which layer at a specific point (x,z) given the interfaces.
% Interfaces must be somewhat horizontal i.e. z_interface = z(x) and sorted
% from the lowest to the highest, and they do not cross each other.
%
% INPUT:
% [x z]         location of a point
% itfs          array of interface structs
% mt            interpolation method
% 
% OUTPUT:
% layer         the layer at (x, z). It is 0 if the point is below the 
%               lowest interface
%
% SEE ALSO:
% MAKEPARAMS, WRITEINTERFACEFILE
%
% Last modified by sirawich@princeton.edu, 08/08/2021

n = length(itfs);
layer = 0;
for ii = 1:n
    z_interface = interp1(itfs{ii}.pts(:,1), itfs{ii}.pts(:,2), pt(1), mt);
    if z_interface > pt(2)
        return
    end
    layer = layer + 1;
end
end