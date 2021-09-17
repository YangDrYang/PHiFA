classdef clDrag
    %clDrag This class represents a hit by wind
    
    properties
        % a-priori information for init
        affectedArea double = -1 % [m^2]
        hitpos(3,1) double % Position of Hit in local coordinate system of target object [m]
        angleOfIncidence double = 0 % [rad]
        
        % call init to calculate
        distFromLaser double % used to calculate closest hit - beam method [m]
        las2hitpos(3,1) double % Position of laser in local coordinate system of target object seen from hit [m]
        
        force(3,1) double % force produced by drag [N]
    end
    
    methods
        function obj = clDrag()
            %clHit Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = init(obj, surfAtr, ray, facet, dens)
            obj.las2hitpos = obj.hitpos - ray.origin;
            obj.distFromLaser = norm(obj.las2hitpos);
            
            %drag
            velocity=ray.origin;
            e_n = facet.normal;
            velonorm = norm(velocity);
            e_v = velocity/velonorm;
            
            obj.force = - 2*dens*velonorm^2 ...
                * ( (2-surfAtr.sigma_n-surfAtr.sigma_t)*(dot(e_v,e_n))^2 * e_n ...
                + surfAtr.sigma_t*dot(e_v,e_n) * e_v ) ...
                * obj.affectedArea; % \citep[eq. 7.4]{wheeler1965magnetic}
        end
    end
end

