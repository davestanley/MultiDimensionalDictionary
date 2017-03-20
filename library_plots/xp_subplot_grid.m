

function hsg = xp_subplot_grid (xp, op)
	% This handles 1D or 2D xp data. For 3D data see xp_subplot_grid3D.
    
    if nargin < 2
        op = struct;
    end
    
    if isempty(op); op = struct; end;
    
    op = struct_addDef(op,'transpose_on',0);
    op = struct_addDef(op,'display_mode',0);
    op = struct_addDef(op,'subplotzoom_enabled',1);
    op = struct_addDef(op,'legend1',[]);
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
                          
    transpose_on = op.transpose_on;
    display_mode = op.display_mode;
    subplotzoom_enabled = op.subplotzoom_enabled;
    legend1 = op.legend1;
    
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
                    if i == 1 && j == 1 && ~isempty(legend1)
                        % Place a legend in the 1st subplot
                        legend(legend1{:});
                    end
                end
            end
            
            % Do labels for rows
            if ~strcmp(xp.axis(1).name(1:3),'Dim')          % Only display if its not an empty axis
                rowstr = setup_axis_labels(xp.axis(1));
                hsg.rowtitles(rowstr);
            end
            
            % Do labels for columns
            if ~strcmp(xp.axis(2).name(1:3),'Dim')          % Only display if its not an empty axis
                colstr = setup_axis_labels(xp.axis(2));
                hsg.coltitles(colstr);
            end
            
            
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
