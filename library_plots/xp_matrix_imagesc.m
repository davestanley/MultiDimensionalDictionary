

function hxp = xp_matrix_imagesc (xp,options)
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
    if ~isfield(options,'zlims'); options.zlims = []; end
    if ~isfield(options,'do_colorbar'); options.do_colorbar = false; end
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
                          
    transpose_on = options.transpose_on;
    xlims = options.xlims;
    ylims = options.ylims;
    xdat = options.xdat;
    ydat = options.ydat;
    zlims = options.zlims;
    do_colorbar = options.do_colorbar;
    
    
    if transpose_on
        d = xp.data{1}';
    else
        d = xp.data{1}';
    end
    
    if ~isempty(zlims)
        hxp.hcurr = imagesc(xdat,ydat,d,zlims);
    else
        hxp.hcurr = imagesc(xdat,ydat,d);
    end
    set(gca,'YDir','normal');
    
    if ~isempty(xlims); xlim(xlims); end
    if ~isempty(ylims); ylim(ylims); end
    
    if do_colorbar
        colorbar;
    end
    
    
end


