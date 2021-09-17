classdef clTarget < cl6dofObject
    %clTarget Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         mass double = -1 % [kg]
        volume double = -1 % [m^3]
        shape eShape
    end
    
    methods
        function obj = clTarget()
            %clTarget Construct an instance of this class
        end
    end
    
    methods (Abstract)
        init(obj, vargin)
%         [ impulse, momentum ] = getLaserAccel(obj, laserobj)
        move(obj, varargin)
        spin(obj, varargin)
        getLongestSide(obj)
    end
end

