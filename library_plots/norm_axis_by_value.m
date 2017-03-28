function obj = norm_axis_by_value(obj, axis_dim, axis_value)

if ischar(axis_dim)
    dim_string = axis_dim;
    axis_dim = obj.findaxis(dim_string);
    if ~isscalar(axis_dim) || isempty(axis_dim)
        error('Multiple or zero dimensions matching %s.', dim_string)
    end
end
            
if ~isscalar(axis_dim) || isempty(axis_dim) || axis_dim == 0
    error('Dimension to normalize must be a nonempty, nonzero scalar.')
end

obj_denom = squeeze(obj.axissubset(axis_dim, axis_value));

obj_denom = obj_denom.repmat(obj.axis(axis_dim).values, obj.axis(axis_dim).name, axis_dim);

obj.data = cellfun(@(x,y) 100*(x./y - 1), obj.data, obj_denom.data, 'UniformOutput', 0);

end