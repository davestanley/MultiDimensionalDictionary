

function xp_PlotData (xp, op)
    % xp must be 1x1 (e.g. 0 dimensional)
    if nargin < 2
        op = struct;
    end
    
    if isempty(op); op = struct; end;
    
    op = struct_addDef(op,'args',{});
    
    xlims = op.xlims;
    ylims = op.ylims;
    
    % If doing rastergrams or raster, re-add axis for populations
    plot_type = [];
    ind = find(strcmp(op.args,'plot_type'));
    if ~isempty(ind); plot_type = op.args{ind+1}; end
    if any(strcmp(plot_type,{'rastergram','raster'}))
        if isfield(xp.meta,'matrix_dim_3')
            if strcmp(xp.meta.matrix_dim_3.name,'populations')
                %Nd = ndims(xp.data{1});
                Nd = 3;
                xp = xp.unpackDim(Nd);
            end
        end
    end

    
    % Add dummy axis for variables
    Na = length(xp.axis);
    xp.axis(Na+1).name = 'variables';
    xp.axis(Na+1).values = {'v'};
    
    % Convert xp to DynaSim data struct
    data = xPlt2DynaSim(xp);
    
    % Feed into original PlotData command, making sure it doesn't generate
    % new figures (rather, should produce it in the current subplot)
    PlotData(data,op.args{:},'lock_gca',true);
    
    if ~isempty(xlims); xlim(xlims); end
    if ~isempty(ylims); ylim(ylims); end

end


