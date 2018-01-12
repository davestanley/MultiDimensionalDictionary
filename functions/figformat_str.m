function output = figformat_str(input)
    % Formats string or cell array of strings for plotting. (e.g. escapes underscores).

    func1 = @(s) strrep(s,'_','\_'); % escape underscores

    if iscellstr(input)
        output = cellfun(func1,input,'Uniform_Output',0);
    elseif ischar(input)
        output = func1(input);
    else
        error('Unknown input type');
    end
        

end
