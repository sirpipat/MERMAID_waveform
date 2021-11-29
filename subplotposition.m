function pos = subplotposition(m, n, p, axmargin, figmargin)
% pos = SUBPLOTPOSITION(m, n, p, axmargin, figmargin)
%
% Calculate a subplot's position in m-by-n matrix of small axes, selects
% the p-th axes for the current plot with more adjustable axes and figure 
% margins than SUBPLOT(m,n,p).
%
% INPUT:
% m             number of rows
% n             number of columns
% p             indexing as in SUBPLOT(m,n,p), must be a scalar integer
% axmargin      axes margin     [left bottom right top] in relative size to
%               the whole subplot area (default: [0.15 0.15 0.15 0.15])
% figmargin     figure margin   [left bottom right top] in relative size to
%               the whole figure area (default: [0.03 0.03 0.03 0.03])
%
% OUTPUT:
% pos           position of the axes    [left bottom width height]
%
% SEE ALSO:
% SUBPLOT
%
% Last modified by sirawich-at-princeton.edu, 11/29/2021

defval('axmargin', [0.15 0.15 0.15 0.15])
defval('figmargin', [0.03 0.03 0.03 0.03])

gridpos = [mod(p-1,n) n-ceil(p/n)];
widthfull = (1-figmargin(1)-figmargin(3))/m;
heightfull = (1-figmargin(2)-figmargin(4))/n;
pos = [figmargin(1) + gridpos(1) * widthfull ...
       figmargin(2) + gridpos(2) * heightfull ...
       widthfull ...
       heightfull];
axmargin = [axmargin(1) ...
            axmargin(2) ...
            -(axmargin(1)+axmargin(3)) ...
            -(axmargin(2)+axmargin(4))];
axmargin = axmargin .* [widthfull heightfull widthfull heightfull];
pos = pos + axmargin;
end