

function xp_matrix_imagesc (xp)
    % xp must be 1x1 (e.g. zero dimensional)
    
    meta = xp.meta;
    
    for d = 1:2
        dim_name = ['matrix_dim_' num2str(d)];
        if isfield(meta, dim_name)
            axis_labels{d} = meta.(dim_name).name;
            axis_values{d} = meta.(dim_name).values;
        else
            axis_labels{d} = '';
            axis_values{d} = 1:size(xp.data{1}, d);
        end
    end
    
    imagesc(axis_values{1}, axis_values{2}, xp.data{1}');
    
    axis xy
    %colorbar
    
end


