function getdatapoint(src, ax, data, name)
% GETDATAPOINT(src, ~ ,ax, data, name)
%
% Display a data point on the scatter plot. Call this function to the 
% scatter plot's button down function. Then, click one of the point to get
% the (x,y) coordinate as well as the index of the data point. If DATA is
% specified, the value at the index will be displayed as well.
%
% INPUT:
% src       selected object
% ax        parent axes of the selected object
% data      data structure containing the plotted data [optional]
%           it could be a matrix, cell arry, or struct.
% name      name of the data [optional]
%
% EXAMPLE:
%
% data = load('patients.mat');
% obj = scatter(data.Age, data.Diastolic);
% obj.ButtonDownFcn = @(src,event)getdatapoint(src, obj.Parent, data, ...
%    'patient');
%
% % Then click at one data point. Below is an example
%
% x : 44.8013
% y : 92.9234
% ii: 99
% - - - - - - - - - - - - - - - - - - - - 
% patient
% -------
%   Gender: Male
%   LastName: Diaz
%   Age: 45
%   Weight: 172
%   Smoker: TRUE
%   Systolic: 136
%   Diastolic: 93
%   Height: 68
%   Location: County General Hospital
%   SelfAssessedHealthStatus: Good
% ----------------------------------------
%
% Last modified by sirawich-at-princeton.edu: 07/24/2023

defval('data', [])
defval('name', 'value')

pt = ax.CurrentPoint(1,1:2);
ii = knnsearch([src.XData' src.YData'], pt);
fprintf('x : %g\n', pt(1));
fprintf('y : %g\n', pt(2));
fprintf('ii: %g\n', ii);
if ~isempty(data)
    fprintf('- - - - - - - - - - - - - - - - - - - - \n');
    displaydata(data, name, ii);
end
fprintf('----------------------------------------\n');
end

function displaydata(data, fieldname, dndex)
defval('fieldname', 'value')
if ~contains(fieldname, ' ')
    nspace = 0;
else
    nspace = indeks(strfind(fieldname, ' '), 'end');
end
if isempty(data)
    return
end
if isa(data, 'struct')
    fprintf('%s\n', fieldname);
    fprintf([repmat(' ', 1, nspace) repmat('-', 1, length(fieldname) - nspace) '\n']);
    names = fieldnames(data);
    for ii = 1:length(names)
        displaydata(data.(names{ii}), ...
            [repmat(' ', 1, nspace+2) names{ii}], dndex);
    end
else
    if dndex <= length(data)
        if isa(data, 'cell')
            value = data{dndex};
        else
            value = data(dndex);
        end
        if isa(value, 'double')
            if strcmp(replace(fieldname, ' ', ''), 'USER7')
                fprintf('%s: %i\n', fieldname, value);
            else
                fprintf('%s: %g\n', fieldname, value);
            end
        elseif isa(value, 'char') || isa(value, 'string')
            fprintf('%s: %s\n', fieldname, value);
        elseif isa(value, 'logical')
            if value
                fprintf('%s: TRUE\n', fieldname);
            else
                fprintf('%s: FALSE\n', fieldname);
            end
        else
            fprintf('%s: CANNOT PRINT\n', fieldname);
        end
    end
end
end