

function xp_matrix_basicplot (xp)
    % xp must be 1x1 (e.g. zero dimensional)
    xp_dims = sort(size(xp), 2, 'descend');
    if xp_dims(1) ~= 1
        error('xp_matrix_imagesc can only be used with a scalar xp object.')
    end
    
    plot(xp.data{1});
    
end


