function OUT = iscellnum(IN)
% ISCELLNUM(S) returns 1 if IN is a cell array of numerics and 0
%   otherwise.

    OUT = all(cellfun(@isnumeric,IN(:)));

end