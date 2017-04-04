function obj_out = mean_over_axis(obj, axis_name, function_handle, varargin)

    if nargin < 3, function_handle = []; end
    
    if isempty(function_handle), function_handle = @mean; end
    
    if nargin < 4, varargin = {}; end

    %% Finding first singleton dimension.

    data_dims = cellfun(@(x) length(size(x)), obj.data);
    
    max_dim = max(data_dims(:));
    
    for d = 1:max_dim
        
        internal_sz_d = cellfun(@(x) size(x, d), obj.data);
        
        internal_sz(:, d) = internal_sz_d(:);
    
    end
    
    number_non_singleton_cells = sum(internal_sz > 1);
    
    number_non_singleton_cells(end + 1) = 0;
    
    first_all_singleton = find(number_non_singleton_cells == 0, 1, 'first');
   
    %% Packing axis to take mean over.
    
    obj_out = obj.packDim(axis_name, first_all_singleton);
    
    %% Taking mean.
    
    varargin{end + 1} = first_all_singleton;
    
    obj_out.data = cellfun(@(x) feval(function_handle, x, varargin{:}), obj_out.data, 'UniformOutput', 0);
    
    obj_out = squeeze(obj_out);
    
end