classdef xPlt < nDDict
    
    % xPlt class inherets from the nDDict class (AKA multidimensional
    % dictionaries). xPlt adds some plotting functionality (see
    % xPlt.RecursivePlot)
    
    methods
        
        %% Constructor
        function obj = xPlt(varargin)
            obj = obj@nDDict(varargin{:});
            obj.meta.datainfo = nDDictAxis; % #isitoutdated
        end
        
    end
    
end