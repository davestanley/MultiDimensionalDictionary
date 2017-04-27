

function xp_matrix_imagesc (xp, transpose_on, colorbar_on)
    % xp must be 1x1 (e.g. zero dimensional)
    xp_dims = sort(size(xp), 2, 'descend');
    if xp_dims(1) ~= 1
        error('xp_matrix_imagesc can only be used with a scalar xp object.')
    end
    
    if nargin < 2, transpose_on = []; end
    
    if isempty(transpose_on), transpose_on = 0; end
    
    if transpose_on, xp = xp_matrix_transpose(xp); end
    
    if nargin < 3, colorbar_on = []; end
    
    if isempty(colorbar_on), colorbar_on = 0; end
    
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
        
        if iscell(axis_values{d})
            axis_values{d} = 1:size(xp.data{1}, d);
        end
    end
    
    imagesc(axis_values{1}, axis_values{2}, xp.data{1}')
    
    axis xy
    
    if colorbar_on, colorbar; end
    
end


