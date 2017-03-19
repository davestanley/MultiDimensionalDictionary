

function hsg = xp_subplot_grid (xp, options)
	% This handles 1D or 2D xp data. For 3D data see xp_subplot_grid3D.
    
    if nargin < 2
        options = struct;
    end
    
    if isempty(options); options = struct; end;
    
    if ~isfield(options,'transpose_on'); options.transpose_on = 0; end
    if ~isfield(options,'display_mode'); options.display_mode = 0; end
    if ~isfield(options,'xlims'); options.xlims = []; end
    if ~isfield(options,'ylims'); options.ylims = []; end
    if ~isfield(options,'subplotzoom_enabled'); subplotzoom_enabled = []; end
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
                          
    transpose_on = options.transpose_on;
    display_mode = options.display_mode;
    xlims = options.xlims;
    ylims = options.ylims;
    subplotzoom_enabled = options.subplotzoom_enabled;
    
    if verLessThan('matlab','8.4') && display_mode == 1; warning('Display_mode==1 might not work with earlier versions of MATLAB.'); end
    if transpose_on && ismatrix(xp)
        xp = xp.transpose;
    elseif transpose_on && ~ismatrix(xp.data)
        error('xp must be a matrix (e.g. ndims < 3) in order to use transpose');
    end
    
    % Parameters
    %subplot_grid_options = {'no_zoom'};
    subplot_grid_options = {};
    
    sz = size(xp);
    
    if ndims(xp.data) <= 2
        N1 = sz(1);
        N2 = sz(2);
        
        
            if display_mode == 1 
                h0 = gcf; ha0 = gca;
                h = figure('visible','off');
            else
                %figure;
            end
            
            if subplotzoom_enabled
                hsg = subplot_grid(N1,N2,subplot_grid_options{:});
            else
                hsg = subplot_grid(N1,N2,'no_zoom',subplot_grid_options{:});
            end
            c=0;
            for i = 1:N1
                for j = 1:N2
                    c=c+1;
                    hsg.set_gca(c);
                    xp.data{i,j}(); 
                end
            end
            

            
            % Do labels for rows
            rowstr = setup_axis_labels(xp.axis(1));
            hsg.rowtitles(rowstr);
            
            % Do labels for columns
            colstr = setup_axis_labels(xp.axis(2));
            hsg.coltitles(colstr);
            
            
            if display_mode == 1
                
                cdata = print(h,'-RGBImage');
                close(h);

                % Restore original axes and display image
                figure(h0); axes(ha0);
                imshow(cdata);
                
            end
        
        
    elseif ndims(xp.data) == 3
        error('For 3D xp data, use instead xp_subplot_grid3D');
        
    end
    
    
    
    
end

function outstr = setup_axis_labels(xpa)
    vals = xpa.getvaluescellstring;
    vals = strrep(vals,'_',' ');
    outstr = cell(size(vals));
    for j = 1:length(outstr)
        outstr{j} = {'',vals{j}};
    end
    outstr{round(end/2)}{1} = strrep(xpa.name,'_',' ');
end
