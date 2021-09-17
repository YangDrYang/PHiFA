classdef clDrag
    %clDrag This class represents a hit by wind
    
    properties
        % a-priori information for init
        affectedArea double = -1 % [m^2]
        projectedArea double = -1 % [m^2]
        hitpos(3,1) double % Position of Hit in local coordinate system of target object [m]
        angleOfIncidence double = 0 % [rad]
        
        % call init to calculate
        distFromLaser double % used to calculate closest hit - beam method [m]
        las2hitpos(3,1) double % Position of laser in local coordinate system of target object seen from hit [m]
        
        force(3,1) double = [0;0;0]; % force produced by drag [N]
        torque(3,1) double = [0;0;0]
    end
    
    methods
        function obj = clDrag()
            %clHit Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = init(obj, ~, seg_obj, ray, facet, dens, w, varargin)
            surfAtr = seg_obj.surfAtr;
            obj.las2hitpos = obj.hitpos - ray.origin;
            obj.distFromLaser = norm(obj.las2hitpos);
            
            %drag
            velocity=-(ray.origin + cross(obj.hitpos,w));
            e_n = facet.normal;
            velonorm = norm(velocity);
            e_v = velocity/velonorm;
            
            obj.force = -dens*velonorm^2 ...
                * ( (2-surfAtr.sigma_n-surfAtr.sigma_t)*(dot(e_v,e_n))^2 * e_n ...
                + surfAtr.sigma_t*dot(e_v,e_n) * e_v ) ...
                * obj.affectedArea; % \citep[eq. 7.4]{wheeler1965magnetic}
            obj.torque = cross( obj.hitpos, obj.force );
            obj.projectedArea = obj.affectedArea*cos(obj.angleOfIncidence);
        end
        
        function r = plus(obj1,obj2)
            r.affectedArea = obj1.affectedArea + obj2.affectedArea;
            r.hitpos = obj1.hitpos + obj2.hitpos;
            r.angleOfIncidence = obj1.angleOfIncidence + obj2.angleOfIncidence;
            r.distFromLaser = obj1.distFromLaser + obj2.distFromLaser;
            r.las2hitpos = obj1.las2hitpos + obj2.las2hitpos;
            r.force = obj1.force + obj2.force;
            r.torque = obj1.torque + obj2.torque;
        end
    end
end

