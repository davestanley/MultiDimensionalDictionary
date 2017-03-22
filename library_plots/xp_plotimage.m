

function xp_plotimage (xp,op)
    % xp must be 1D
    
    if nargin < 2
        op = struct;
    end

    if isempty(op); op = struct; end;
    
    op = struct_addDef(op,'scale',[]);
    op = struct_addDef(op,'saved_fignum',1);
    
    
    if iscellstr(xp.data{1})
        rgb = imread(xp.data{1}{op.saved_fignum});
    elseif ischar(xp.data{1})
        rgb = imread(xp.data{1});
    end
    
    if ~isempty(op.scale)
        rgb = imresize(rgb,op.scale);
    end
    
    imshow(rgb);
    
end


