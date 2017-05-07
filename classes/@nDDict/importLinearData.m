function obj = importLinearData(obj,X,varargin)
%% obj = importLinearData(data,varargin)
%     Purpose:
%     Imports a linear array of data, converts it into a matrix
%     based on the supplied axislabels, and stores it in xp.data.
%     Also populates the xp.axis.values appropriately.
% 
%     Forms:
%     xp = importLinearData(data,axis_values1,...,axis_valuesN)
% 
%     Inputs:
%     data - vector containing input data. Can be numeric or cell array.
%     axis_values1 - vector containing axis values/labels for dimension1 in data
%     ...
%     axis_valuesN - vector containing containing axis values/labels for dimensionN in data
%
%     NOTE: axis_values 1:N must be either numeric array, cell arrays of numerics,
%           or cell arrays of character vectors

    % Initialize
    axlinear = varargin;
    lenX = length(X);
    Ndims = length(axlinear);

% %     Don't need this option for now.
%     % Check if final argument in varargin is a name/value pair
%     if ischar(varargin{end-1})
%         if strcmp(varargin{end-1},'format')
%             outputformat = varargin{end};
%             axeslinear = axeslinear(1:end-2);       % Remove the name-value pair from axeslinear
%         end
%     end

    % Error checking - X must be linear
    if ~isvector(X); error('data must be linear'); end

    % Error checking - X must be cell or numeric
    [~, XsimpleFormat] = nDDict.calcClasses(X,'data');
    if strcmp(XsimpleFormat,'unknown'); error('data must be a numeric or cell array'); end

    % Error checking - each entry in axislinear must be either numeric or
    % cell. If it's a cell, all entries must char.
    axLinearFormat = cell(1, Ndims);
    for k = 1:Ndims
        axLinearFormat{k} = nDDict.calcClasses(axlinear{k}, 'axis_values');
    end
    if any(strcmp(axLinearFormat,'unknown')) || isempty(axlinear); error('axis_values must be a numeric array, cell array of numerics, or cell array of chars'); end

    % Set up xp.axis_pr
    sz = zeros(1, Ndims);
    for iDim = 1:Ndims
        if strcmp(axLinearFormat{iDim}, 'cellnum')
            axlinear{iDim} = [axlinear{iDim}{:}]; % convert cellnum to numeric array
            fprintf('  Note: Converting dim %i axis_values to numeric array from cell array of numerics\n', iDim)
        end
        
        obj.axis_pr(iDim).values = unique(axlinear{iDim},'stable');
        sz(iDim) = length(obj.axis_pr(iDim).values);

        if isnumeric(axlinear{iDim}(1))
            if any(isnan(axlinear{iDim})) || any(isinf(axlinear{iDim}))
                error('Axis cannot contain NaNs or Infs');
            end
        end
    end
    
    if length(sz) == 1; sz(2) = 1; end

    % Set up target matrix
    switch XsimpleFormat
        case 'cell'
            obj.data_pr = cell(sz);
%         case 'string'
%             xp.data_pr = repmat(string(''),sz);
        case 'numeric'
            obj.data_pr = zeros(sz);
        otherwise
            error('Case not implemented');
    end

    % Set up xp.data_pr -> Convert linear data into a multi dimensional matrix
    for indLinear = 1:lenX
        % Get subscripts
        subs = cell(1,Ndims);
        for iDim = 1:Ndims
            if iscellstr(axlinear{iDim})
                subs{iDim} = find(strcmp(axlinear{iDim}{indLinear},obj.axis_pr(iDim).values));
            else
                subs{iDim} = find(axlinear{iDim}(indLinear) == obj.axis_pr(iDim).values);
            end
        end

        % Add data to sparse cell array or matrix based on subscripts
            % Note: Need to find a good way for dealing with duplicate
            % rows. Right now, default behavior is to overwrite
        obj.data_pr(subs{:}) = X(indLinear);

    end
    
    obj = obj.fixAxes(1);

end
