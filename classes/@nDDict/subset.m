function [obj2, ro] = subset(obj,varargin)
% Purpose: get subset of data based on indicies for numeric or regular
%          expression for cellstring. also gets used by valSubset to index numerics 
%          based on values instead of indicies.
%
% Inputs:
%   varargin: each argument corresponds to each axis.
%       1) [] or ':' for all indicies
%       2) regular expression string for cellstring axis values
%       3) numeric array or cellnum of indicies for numeric or cellnum axis
%          values
%
% Outputs:
%   obj2: object with subset of data
%   ro:  if regular expressions are used, contains regular expressions and 
%        results of regexp 'start' indicies.

% Verify that size of obj is correct
checkDims(obj);

% Check for numericsAsValuesFlag
if ischar(varargin{end}) && strcmp(varargin{end}, 'numericsAsValuesFlag')
    numericsAsValuesFlag = true; % tells subset to use numerics as values
    varargin(end) = []; % remove string
%     nargin = nargin - 1; % reduce num input arguments
else
    numericsAsValuesFlag = false;
end

% Get params to validate that "selection" is the right size
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

axClasses = getclass_obj_axis_values(obj);

% Convert selection to index if using numeric values
if numericsAsValuesFlag
    for i = 1:Ns
        if isnumeric(selection{i})
            thisAxVals = obj.axis_pr(i).values;
            if strcmp('cellnum', axClasses{i})
                thisAxVals = [thisAxVals{:}]; % convert to numeric array
            end
            [~, selection{i}] = intersect(thisAxVals, selection{i});
        end
    end
end

ro = {};
re = '([\d.]*\s*[<>=]{0,2})\s*[a-z_A-Z]?\s*([<>=]{0,2}\s*[\d.]*)';
for i = 1:Ns
    if ischar(selection{i})
        if strcmp(axClasses{i}, 'cellstr')
            % Convert selection to index if using regular expressions
            ro{i}(1,:) = obj.axis_pr(i).values;
            [selection{i} ro{i}(2,:)] = nDDict.regex_lookup(obj.axis_pr(i).values, selection{i});
        elseif numericsAsValuesFlag
            % Convert expression on vals for numerics to indicies
            
            thisAxVals = obj.axis_pr(i).values;
            if strcmp('cellnum', axClasses{i})
                thisAxVals = [thisAxVals{:}]; % convert to numeric array
            end
            
            % Parse expression
            expr = selection{i};
            tokens = regexpi(expr, re, 'tokens');
            tokens = tokens{1}; % enter outer cell
            tokens = regexprep(tokens, '\s',''); % remove whitespace
            tokens = tokens(~cellfun(@isempty, tokens)); % remove empty tokens
            if ~all(cellfun(@isempty, cellfunu(@str2num, tokens)))
                tokens = {[tokens{:}]}; % cat since single expression got split up
            end
            lhsBool = isstrprop(cellfun(@(x) x(1), tokens, 'uniform', false), 'digit');
            lhsBool = [lhsBool{:}];
            
            % check LHS and RHS for up to 2 tokens
            parsedSelection = true(size(thisAxVals));
            for k=1:length(tokens)
                thisExpr = tokens{k};
                if lhsBool(k)
                    parsedSelection = parsedSelection & eval(['(' thisExpr 'thisAxVals)']);
                else
                    parsedSelection = parsedSelection & eval(['(thisAxVals' thisExpr ')']);
                end
            end
            [~, selection{i}] = intersect(thisAxVals, thisAxVals(parsedSelection));
        end
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