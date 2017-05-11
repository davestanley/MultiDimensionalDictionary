classdef xPlt < MDDict
    
    % xPlt class inherets from the MDDict class (AKA multidimensional
    % dictionaries). xPlt adds some plotting functionality (see
    % xPlt.RecursivePlot)
    
    methods
        
        %% Constructor
        function obj = xPlt(varargin)
            obj = obj@MDDict(varargin{:});
            obj.meta.datainfo = MDDictAxis; % #isitoutdated
        end
        
    end
    
end