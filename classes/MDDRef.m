%% MDDRef - a copyable handle/reference class wrapper for MDD
%
% Purpose: unlike the base MDD class, this class permits...
%   1) pass-by-reference (to avoid copying entire object when passing to functions as arg)
%   2) event-driven callbacks in a subclass of MDDRef
%
% Usage: MDDRef objects have the same interface as MDD objects.
%
% Note: Due to the wrapper, tab completion show the available MDD methods. 
%       However, those methods are still accessible when used.
%
% Author: Erik Roberts
%
% See also: MDD

classdef MDDRef < matlab.mixin.Copyable
    
    properties (Access=private)
        baseObj
        baseObjClass = MDD
    end
    
    properties (Dependent)
        data
        axis
    end
    
    methods
        function obj = MDDRef(varargin)
          % MDDRef - Default constructor
            %
            % Usage:
            %   obj = MDDRef()
            %   obj = MDDRef(data) % multidimensional data
            %   obj = MDDRef(data, axis_vals, axis_names) % multidimensional or linear data
            %   obj = MDDRef(axis_class, data, axis_vals, axis_names) % for subclassing MDDAxis
            %   obj = MDDRef(baseObjClass, axis_class, data, axis_vals, axis_names) % for subclassing MDD and MDDAxis
            %   obj = MDDRef(baseObjClass, data, axis_vals, axis_names) % for subclassing MDD
            
            if nargin && (isobject(varargin{1}) && any(strcmp(superclasses(varargin{1}), 'MDD')))
                obj.baseObjClass = varargin{1};
                varargin(1) = [];
            end
            
            obj.baseObj = feval(str2func(class(obj.baseObjClass)), varargin{:});
        end
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % Getters % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        function value = get.data(obj)
            value = obj.baseObj.data;
        end
        
        function varargout = get.axis(obj)
            [varargout{1:nargout}] = obj.baseObj.axis;
        end
        
        %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % OVERLOADED OPERATORS % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        function varargout = subsref(obj, S)
            if strcmp(S(1).type, '.') && (strcmp(S(1).subs, 'copy') || strcmp(S(1).subs, 'methods'))
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
        
        function mObj = methods(obj)
            mObj = builtin('methods', obj);
            mBaseObj = methods(obj.baseObj);
            mBaseObj(strcmp(mBaseObj, class(obj.baseObj))) = []; % remove class name
            
            mObj = unique([mObj;mBaseObj]);
        end
        
    end
    
end