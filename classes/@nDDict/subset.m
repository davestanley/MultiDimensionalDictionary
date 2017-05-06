function [obj2, ro] = subset(obj,varargin)
% Define variables and check that all dimensions are consistent
% ro - if regular expressions are used, returns the index
% values discovered by the regular expression.

% Verify that size of obj is correct
checkDims(obj);

% Validate that "selection" is the right size
selection = varargin(:);
Na = length(obj.axis_pr);
Nd = ndims(obj.data_pr);

% Fill out selection with empties if too short. For example, say
% size(obj.data_pr) is MxNx1x1. The user might enter a selection with
% length(selection) = 2, such as selection = {[A],[B]}. In this
% case, we convert it from to selection = {[A],[B],[],[]}.
Ns = length(selection);
if Ns < Na && Ns >= Nd
    selection2 = repmat({[]},1,Na);
    selection2(1:Ns) = selection;
    selection = selection2;
    clear selection2
end

% Trim back selection if too long. If size(obj.data_pr) is MxN,
% it is okay for selection to be {5,3,Y,Z} as long as Y and Z
% are either empty (:) or ones. If they are anything else,
% return error!
Ns = length(selection);
if Ns > Na
    % Make sure extra selection entries are either "1" or []
    selection_extra = selection(Na+1:end);
    are_empties = cellfun(@isempty,selection_extra);
    are_ones = cellfun(@(s) s == 1, selection_extra);
    if any(~(are_empties | are_ones))      % If any of the extra selections are NOT either empty or one's...
        error('Index exceeds dimensions of nDDict.data');
    end
    selection = selection(1:Na);
end

% Replace any ':' entries in selection with []. Empty
% entries code for taking all entries along a dimension; ':' is
% an alias for this.
for i = 1:length(selection)
    if ischar(selection{i})
        if strcmp(selection{i},':')
            selection{i} = [];
        end
    end
end

% If Ns is still wrong dimensions, return error
Ns = length(selection);
if Ns ~= Na
    error('Number of inputs must match dimensionality of nDDict.data');
end

% Convert selection to index if using regular expressions
ro = {};
for i = 1:Ns
    if ischar(selection{i})
        ro{i}(1,:) = obj.axis_pr(i).values;
        [selection{i} ro{i}(2,:)] = nDDict.regex_lookup(obj.axis_pr(i).values, selection{i});
        
    end
end

% Make sure that size of selection doesnt exceed size of data
sz = size(obj);
for i = 1:Ns
    if selection{i} > sz(i); error('Selection index exceeds dimensions'); end
end

% Initialize

obj2 = obj;             % Create new class of same type as original
obj2 = obj2.reset;
obj2.meta = obj.meta;

% First update each axis with the corresponding selection and
% convert empty cells to code for full range.
for i = 1:Ns
    
    if isempty(selection{i})
        selection{i} = 1:sz(i);
    end
    
    obj2.axis_pr(i) = obj.axis_pr(i);       % Import axis information
    obj2.axis_pr(i).values = obj.axis_pr(i).values(selection{i});   % Overwrite values field; leave everything else the same.
end

% Update the data
obj2.data_pr = obj.data_pr(selection{:});

% Corrects number of axes. The above code automatically
% converts obj.data_pr from MxNx1x1 to MxN, whereas axis will stay
% as MxNx1x1 (e.g. length of 4). Thus, fixAxes corrects this.
% This should no longer be necessary!!!
% obj2 = fixAxes(obj2);

end