classdef clFacet < handle & matlab.mixin.Heterogeneous
    %clFacet As in facets of a gem.
    %   Can be any type of planar polygon figure
    
    properties
        % must be set
        base(3,1) double = [0; 0; 0]
        normal(3,1) double = [1; 0; 0]
        
        % call .init() to calculate
        areabarycenter(3,1) double = [0; 0; 0]
        area double = -1 % m^3
        
        spflag logical = false %solar panel indicator
        
        % facet efficiency (takes into account that a facet might not be to
        % 100% made of one material or other geometry factors that are
        % neglected due to simplification of overall geometry)
        fef double = 1
    end
    
    methods
        function obj = clFacet()
            %clFacet Construct an instance of this class
        end
        
        function moveBase(obj, move_by)
            obj.base = obj.base - move_by;
            obj.areabarycenter = obj.areabarycenter - move_by;
        end
    end
    
    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = clTriangle;
        end
    end
    
    methods (Abstract)
        init(obj)
        rayIntersection(obj, satpos, ray)
        getIntersection(obj, satpos, ray, laserstation)
        getReferenceCube(obj)
        compareReferenceCube(obj, refcube)
        getPointNet(obj, resolution)
    end
end

