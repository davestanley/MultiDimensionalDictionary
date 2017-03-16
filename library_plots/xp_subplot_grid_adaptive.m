

function hsg = xp_subplot_grid_adaptive (xp, dim_order, max_subplot_side, display_mode, transpose_on)
	% This handles 1D or 2D xp data. For 3D data see xp_subplot_grid3D.
    
    if nargin < 5
        transpose_on = [];
    end
    
    if nargin < 4
        display_mode = [];
    end
    
    if nargin < 3, max_subplot_side = []; end
    
    if nargin < 2, dim_order = []; end
    
    if isempty(transpose_on), transpose_on = 0; end
    
    if isempty(display_mode), display_mode = 0; end
            % Display_mode: 0-Just plot directly
                          % 1-Plot as an image (cdata)
                          % 2-Save to a figure file 
    
    if verLessThan('matlab','8.4') && display_mode == 1 
        warning('display_mode==1 might not work with earlier versions of MATLAB.')
    end
    
    if isempty(max_subplot_side), max_subplot_side = 15; end
    
    % Parameters
    %subplot_grid_options = {'no_zoom'};
    subplot_grid_options = {};
    
    sz = size(xp);
    
    if isempty(dim_order), [~, dim_order] = sort(sz, 2, 'descend'); end
    
    if iscellstr(dim_order)
        dim_order_cell = dim_order;
        dim_order = nan(size(dim_order_cell));
        for d = 1:length(dim_order)
            dim_order(d) = xp.axis.findAxes(dim_order_cell{d});
        end
    end
        
    xp = permute(xp, dim_order);
    
    sz = size(xp);
    
    [dim_indices, subplot_indices, figure_indices, no_rows, no_cols, figs_through] = adaptive_indices(sz, max_subplot_side);
    
    new_fig_indices = diff([0; figure_indices]);
    
    open_figures = findall(0, 'Type', 'figure');
    
    if ~isempty(open_figures)
        last_figure = max([open_figures(:).Number]);
    else
        last_figure = 0;
    end
    
    for plot = 1:size(dim_indices, 1)
        
        fig_for_plot = figure_indices(plot);

        if new_fig_indices(plot)
            
            figure(last_figure + fig_for_plot)
            
            % if display_mode == 1
            %     h0 = figure(figure_indices(plot)); ha0 = gca;
            %     h = figure(figure_indices(plot), 'visible', 'off');
            % else
            %     h0 = figure(figure_indices(plot));
            % end
            
            hsg(fig_for_plot) = subplot_grid(no_rows, no_cols, subplot_grid_options{:});
            
        end
        
        hsg(fig_for_plot).set_gca(subplot_indices(plot));
        xp.data{dim_indices{plot, :}}();
        
        % Titles for subplots or rows/columns.
        
        if figs_through(2) > 1
            
            title([figformat_str(xp.axis(1).name) ': ' figformat_str(xp.axis(1).getvaluestring(plot))])
            
        else
            
            if subplot_indices(plot) == 1
                
                % Do labels for rows
                rowstr = setup_axis_labels(xp.axis(1));
                hsg(fig_for_plot).rowtitles(rowstr);
                
                % Do labels for columns
                colstr = setup_axis_labels(xp.axis(2));
                hsg(fig_for_plot).coltitles(colstr);
                
            end
            
        end
        
        if subplot_indices(plot) == 1
            
            % Overall figure title.
            
            start_title_axes = find(cumprod(figs_through) > 1, 1, 'first');
            
            title_axis = xp.axis(start_title_axes:end);
            mytitle = '';
            for a = 1:length(title_axis)
                mytitle = [mytitle, figformat_str(title_axis(a).name) ': '...
                    figformat_str(title_axis(a).getvaluestring(dim_indices{plot, start_title_axes + a - 1})), ' '];
            end
            hsg(fig_for_plot).figtitle(mytitle);
            
        end
        
    end
    
    % if display_mode == 1
    % 
    %     cdata = print(h,'-RGBImage');
    %     close(h);
    % 
    %     % Restore original axes and display image
    %     figure(h0); axes(ha0);
    %     imshow(cdata);
    % 
    % end
        
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
