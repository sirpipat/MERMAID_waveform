function struct_arr = array2struct(arr_struct)
% struct_arr = ARRAY2STRUCT(arr_struct)
%
% Converts an array of structs to a struct with field varaiables as arrays.
% Note that the type or arrays in the output depends on the variable types
% of the struct's fields.
% 
% Array type
% 1. cell array      -- strings, character arrays, arrays with length > 1
% 2. matrix          -- int, single, double, datetime (anything that can be
%                       nan)
% 3. struct of array -- struct (The function will be called recursively.)
%
% INPUT:
% arr_struct        array of structs
%
% OUTPUT:
% struct_arr        struct with the same fields as struct in arr_struct
%
% Last modified by sirawich-at-princeton.edu, 03/01/2024

obj = arr_struct(1);
objfields = fieldnames(obj);

% construct the header array
for ii = 1:length(objfields)
    if or(or(isstring(obj.(objfields{ii})), ...
            ischar(obj.(objfields{ii}))), length(obj.(objfields{ii})) > 1)
        struct_arr.(objfields{ii}) = cell(size(arr_struct));
    elseif isstruct(obj.(objfields{ii}))
        struct_arr.(objfields{ii}) = repmat(obj.(objfields{ii}), ...
            size(arr_struct));
    else
        struct_arr.(objfields{ii}) = nan(size(arr_struct));
    end
end

% list of field names
for ii = 1:length(arr_struct)
    obj = arr_struct(ii);
    for jj = 1:length(objfields)
        fieldval = obj.(objfields{jj});
        if or(or(ischar(fieldval), isstring(fieldval)), ...
                length(fieldval) > 1)
            eval(sprintf('struct_arr.%s{%d} = fieldval;', ...
                objfields{jj}, ii));
        else
            eval(sprintf('struct_arr.%s(%d) = fieldval;', ...
                objfields{jj}, ii));
        end
    end
end

% repeat if any field is a struct
for ii = 1:length(objfields)
    if isstruct(obj.(objfields{ii}))
        struct_arr.(objfields{ii}) = array2struct(struct_arr.(objfields{ii}));
    end
end
end