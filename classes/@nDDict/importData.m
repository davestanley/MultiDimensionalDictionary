function obj = importData(obj,data,axis_vals,axis_names)
            
obj.data_pr = data;
obj = obj.fixAxes;

if nargin > 2
    if ~isempty(axis_vals)
        if ~iscell(axis_vals); error('axis_vals must be a cell array.'); end
        for i = 1:length(axis_vals)
            obj.axis_pr(i).values = axis_vals{i};
        end
        obj.checkDims;
    end
end

if nargin > 3
    if ~isempty(axis_names)
        obj = obj.importAxisNames(axis_names);
    end
end

obj.fixAxes(1);     % Convert any axis vallues that are cellnums to numeric matrices

end