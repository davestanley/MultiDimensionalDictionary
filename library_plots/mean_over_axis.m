function obj_out = mean_over_axis(obj, axis_name)

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
    
    obj_out.data = cellfun(@(x) nanmean(x, first_all_singleton), obj_out.data, 'UniformOutput', 0);
    
    obj_out = squeeze(obj_out);
    
end