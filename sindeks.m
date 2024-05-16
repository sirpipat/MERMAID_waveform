function x = sindeks(s, i, n, c)
% x = SINDEKS(s, i, n, c)
%
% Extract linearly indexed position(s) i out of matrix/cell members of a 
% struct s. Works for logical and numeric indices.
%
% INPUT:
% s             The struct containing matricies
% i             The requested set of running linear indices [default: 1]
% n             Expected number of elements of the matrices. This function
%               only extract from matrices of length n.
%               [default: [] -- try to extract from every matrix. May lead
%               to unexpected results.]
% c             Whether to extract char matrix of length n [default: false]
%
% OUTPUT:
% x             The struct with the same structure as s but the matrix
%               members are sliced
%
% EXAMPLES:
% s.a = [1 2 3 4 5];
% s.b = 12;
% s.c = 'apple';
% s.d.aa = {1 2 3 4 5};
% s.d.bb = [6 7 8 9 0];
% 
% sindeks(s, 2:4, 5)
% sindeks(s, 'end-1:end', 5, true)  % See what happens for s.c
% sindeks(s, mod((1:5), 3) == 1, 5)
%
% SEE ALSO:
% INDEKS, CINDEKS, RINDEKS, KINDEKS, TINDEKS, DINDEKS, SQUEEZE
%
% Last modifid by sirawich-at-princeton.edu, 05/15/2024

defval('i', 1)
defval('n', -1)
defval('c', false)

if ~isstruct(s)
    warning('s must be a struct. x = INDEKS(s, i) is called instead.')
    x = indeks(s, i);
end

names = fieldnames(s);

for ii = 1:length(names)
    if n > 0
        % extract the matrices/cell array or char array if allowed of
        % length n
        if (c || ~ischar(s.(names{ii}))) && length(s.(names{ii})) == n
            x.(names{ii}) = rindeks(s.(names{ii}), i);
        % if the length is not matches
        else
            % reiterate over a member struct
            if isstruct(s.(names{ii}))
                x.(names{ii}) = sindeks(s.(names{ii}), i, n, c);
            % copy the variable with no change
            else
                x.(names{ii}) = s.(names{ii});
            end
        end
    else
        % try to extract every variable until getting bad subscript error
        if c || ~ischar(s.(names{ii}))
            try
                x.(names{ii}) = rindeks(s.(names{ii}), i);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:badsubscript')
                    if isstruct(s.(names{ii}))
                        x.(names{ii}) = sindeks(s.(names{ii}), i, n, c);
                    else
                        x.(names{ii}) = s.(names{ii});
                    end
                else
                    rethrow(ME)
                end
            end
        else
            x.(names{ii}) = s.(names{ii});
        end
    end
end

end