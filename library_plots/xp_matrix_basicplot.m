

function xp_matrix_basicplot (xp, options)
    % xp must be 1x1 (e.g. 0 dimensional)
    
    
    if nargin < 2
        options = struct;
    end
    
    if isempty(options); options = struct; end;
    
    if ~isfield(options,'xlims'); options.xlims = []; end
    if ~isfield(options,'ylims'); options.ylims = []; end
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
                          
    xlims = options.xlims;
    ylims = options.ylims;
    
    plot(xp.data{1});
    if ~isempty(xlims); xlim(xlims); end
    if ~isempty(ylims); ylim(ylims); end
    
end


