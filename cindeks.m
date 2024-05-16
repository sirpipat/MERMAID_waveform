function x = cindeks(C, i, j)
% x = CINDEKS(C, i, j)
%
% Extracts linearly indexed position(s) i (1D array) or (i,j) (2D array) 
% out of a cell array C.
%
% INPUT:
%
% C         The input cell array
% i         The requested set of first running linear indices [default: 1]
% j         The requested set of second running linear indices
%
% EXAMPLES:
% % 1D array
% str = 'Mary has a little lamb"
% cindeks(split(str, ' '), 2)
% cindeks(split(str, ' '),'end')
%
% % 2D array
% C = {{'first', 'second'}, {'third', 'forth'}};
% cindeks(C, 1, 2)
% cindeks(C,'end', 1)
%
% See also INDEKS, SINDEKS, RINDEKS, KINDEKS, TINDEKS, DINDEKS, SQUEEZE
%
% Last modified by sirawich-at-princeton.edu, 05/16/2024

defval('i', 1)

switch nargin
    case 2
        if ~ischar(i)
            x = C{i};
        else
            eval([ 'x = C{' i '};'])
        end
    case 3
        if ~ischar(i)
            i = string(i);
        end
        if ~ischar(j)
            j = string(j);
        end
        eval(sprintf('x = C{%s,%s};', i, j));
end