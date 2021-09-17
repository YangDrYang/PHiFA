classdef clObject < handle
    %clObject This class provides an interface to the propagator class, so
    %each class that inherites from this can be easily propagated
    
    properties
        name string = 'default'
        
        xv(6,1) double = [0; 0; 0; 0; 0; 0] % m and m/s ECI
        qw(7,1) double = [0; 0; 0; 0; 0; 0; 0] % rad and rad/s
        %attitude is depicted by quaternion
        currentEpoche double = 58414.0 % MJD
        
        mass double = -1 % [kg]
    end
    
    methods
        function obj = clObject()   
        end
        
        function height = getHeight(obj)
            height = norm(obj.xv(1:3))-6378.137e3;
        end
    end
    
    methods (Abstract)
        dY = Accel(obj, t, Y)
    end
end

