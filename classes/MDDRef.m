%% MDDRef - a copyable handle/reference class wrapper for MDD
%
% Purpose: unlike the base MDD class, this class permits...
%   1) pass-by-value (to avoid copying entire object when passing to functions as arg)
%   2) event-driven callbacks in a subclass of MDDRef
%
% Usage: MDDRef objects have the same interface as MDD objects.
%
% Note: Due to the wrapper, tab completion show the available MDD methods. 
%       However, those methods are still accessible when used.
%
% See also: MDD

classdef MDDRef < matlab.mixin.Copyable
    
    properties (Access=private)
        baseObj
    end
    
    properties (Dependent)
        data
        axis
    end
    
    methods
        %         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %         % % % % % % % % % % % Getter and Setters % % % % % % % % % % % %
        %         % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function obj = MDDRef(varargin)
            obj.baseObj = MDD(varargin{:});
        end
        
%         function set.data(obj, value)
%             obj.baseObj.data = value;
%         end
        
        function value = get.data(obj)
            value = obj.baseObj.data;
        end
        
%         function set.axis(obj,value)
%             obj.baseObj.axis = value;
%         end
        
        function varargout = get.axis(obj)
            [varargout{1:nargout}] = obj.baseObj.axis;
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % OVERLOADED OPERATORS % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function varargout = subsref(obj, S)
            if strcmp(S(1).type, '.') && strcmp(S(1).subs, 'copy')
                [varargout{1:nargout}] = builtin('subsref', obj, S);
            else
                [varargout{1:nargout}] = builtin('subsref', obj.baseObj, S);
                if isa(varargout{1}, 'MDD') && ~strcmp(S(1).type, '()')
                    obj.baseObj = varargout{1};
                end
            end
        end
        
        function obj = subsasgn(obj, S, value)
            obj.baseObj = builtin('subsasgn', obj.baseObj, S, value);
        end
        
    end
    
end