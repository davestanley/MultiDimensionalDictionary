

function xp_matrix_imagesc (xp,options)
    % xp must be 1x1 (e.g. zero dimensional)
    
    
    if nargin < 2
        options = struct;
    end
    
    if isempty(options); options = struct; end;
    
    if ~isfield(options,'transpose_on'); options.transpose_on= 0; end
    if ~isfield(options,'xlims'); options.xlims = []; end
    if ~isfield(options,'ylims'); options.ylims = []; end
    if ~isfield(options,'xdat'); options.xdat= []; end
    if ~isfield(options,'ydat'); options.ydat = []; end
    if ~isfield(options,'clims'); options.clims = []; end
    if ~isfield(options,'show_colorbar'); options.show_colorbar = 0; end
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
                          
    transpose_on = options.transpose_on;
    xlims = options.xlims;
    ylims = options.ylims;
    xdat = options.xdat;
    ydat = options.ydat;
    clims = options.clims;
    show_colorbar = options.show_colorbar;
    
    
    if transpose_on;
        d = xp.data{1}';
    else
        d = xp.data{1}';
    end
    
    if ~isempty(clims)
        imagesc(xdat,ydat,d,clims);
    else
        imagesc(xdat,ydat,d);
    end
    
    if ~isempty(xlims); xlim(xlims); end
    if ~isempty(ylims); ylim(ylims); end
    
    if show_colorbar
        colorbar;
    end
    
    
end


