function obj = importData(obj, data, axis_vals, axis_names)
%% importData - overwrite object data with multidimensional data from variable

% Note: functionality can be called from a static (ie class) or object method
%   obj = ImportData(data,axis_vals,axis_names) % uppercase method
%   obj = importData(obj, data,axis_vals,axis_names) % lowercase method

obj.data_pr = data;
obj = obj.fixAxes;

%% import data
if nargin > 2 && ~isempty(axis_vals)
    if ~iscell(axis_vals); error('axis_vals must be a cell array.'); end
    
    % Handle vector data
    if length(axis_vals) == 1 && ismatrix(obj)% #checkthis
        vecDim = cellfun(@length,obj.exportAxisVals) == length(axis_vals{1});
        axis_vals{vecDim} = axis_vals{1}; % move axis_vals to right dim
        axis_vals{~vecDim} = 1; % set other dim to 1
    end
    
    for i = 1:length(axis_vals)
        obj.axis_pr(i).values = axis_vals{i};
    end
    
    obj.checkDims;
end

if nargin > 3 && ~isempty(axis_names)
    obj = obj.importAxisNames(axis_names);
end

obj.fixAxes(1);     % Convert any axis vallues that are cellnums to numeric matrices

end
