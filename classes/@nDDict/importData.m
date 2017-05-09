function obj = importData(varargin)
%% importData - import multidimensional data

% Note: can be called as a static (ie class) or object method
%   obj = importData(data,axis_vals,axis_names)
%   obj = importData(obj, data,axis_vals,axis_names)

%% parse arguments
[varargin, objClass] = nDDict.nonObjArgs(varargin{:}); % remove possible object arg
nargin = length(varargin); % redefine nargin in case first arg was object
if nargin < 3
    varargin{3} = []; % fill in missing args with []
end
[data, axis_vals, axis_names] = deal(varargin{:});

%% instantiate object
if isempty(objClass) % didn't use object.method syntax
    obj = nDDict();
else % used object.method syntax
    obj = feval(str2func(objClass)); % support subclass
end

obj.data_pr = data;
obj = obj.fixAxes;

%% import data
if nargin > 1 && ~isempty(axis_vals)
    if ~iscell(axis_vals); error('axis_vals must be a cell array.'); end
    
    for i = 1:length(axis_vals)
        obj.axis_pr(i).values = axis_vals{i};
    end
    
    obj.checkDims;
end

if nargin > 2 && ~isempty(axis_names)
    obj = obj.importAxisNames(axis_names);
end

obj.fixAxes(1);     % Convert any axis vallues that are cellnums to numeric matrices

end