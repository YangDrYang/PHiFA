classdef clHit
    %clHit This class represents a hit by the laser beam on a facet
    
    properties
        % a-priori information for init
        affectedArea double = -1 % [m^2]
        projectedArea double = -1 % [m^2]
        hitpos(3,1) double % Position of Hit in local coordinate system of target object [m]
        angleOfIncidence double = 0 % [rad]
        
        % call init to calculate
        distFromLaser double % used to calculate closest hit - beam method [m]
        las2hitpos(3,1) double % Position of laser in local coordinate system of target object seen from hit [m]
        couplingCoef double = -1
        localFluence double = -1 % Joule [m^2]
        intensity double = -1 % Joule [m^2/s]
        
        impuls(3,1) double = [0; 0; 0] % impuls caused by ablation [N*s]
        moment(3,1) double = [0; 0; 0] % ang momentum caused by ablation [N*m*s]
        
        cwforce(3,1) double = [0; 0; 0] % continious force caused by photon pressure [N]
        cwtorque(3,1) double = [0; 0; 0] % torque from cwforce [Nm]
    end
    
    methods
        function obj = clHit()
            %clHit Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function obj = init(obj, surfAtr, ray, facet, laserstation)
            obj.las2hitpos = obj.hitpos - ray.origin;
            obj.distFromLaser = norm(obj.las2hitpos);
            
            %pulsed
            [obj.localFluence, obj.intensity] = laserstation.getBeamPropertiesAtTarget(obj.hitpos);
            targetFluence = obj.localFluence * cos(obj.angleOfIncidence);
            obj.couplingCoef = surfAtr.getCouplingCoeffientAblation(targetFluence);
            obj.impuls = - facet.fef * facet.normal * obj.couplingCoef * targetFluence * obj.affectedArea;
            obj.moment = cross( obj.hitpos, obj.impuls );
            
            %\gls{cw}
            c = 299792458; % [m/s]
%             kstrich = ray.direction - 2*dot(ray.direction,facet.normal)...
%                 *facet.normal; % \citep[eq. 46]{liedahl2013}
%             obj.cwforce = facet.fef * obj.intensity*obj.affectedArea/c * ...
%                 abs(dot(ray.direction, facet.normal)) * ...
%                 ( ray.direction - surfAtr.alpha * surfAtr.beta * kstrich - ...
%                 0.5*surfAtr.alpha * (1-surfAtr.beta) * facet.normal ); % \citep[eq. 47]{liedahl2013}
            
%             obj.cwforce = (ca+crd)*obj.affectedArea*obj.intensity*ray.direction/c + ...
%                 (2*crd/3*obj.affectedArea + ...
%                 2*crs*obj.affectedArea*cos(obj.angleOfIncidence)) *(-facet.normal);
%             obj.cwtorque = cross( obj.hitpos, obj.cwforce ); %wrong codes to calculate srp force, commented by Yang
%             
            obj.projectedArea = obj.affectedArea*cos(obj.angleOfIncidence);
            
            if ~facet.spflag %no solar panel
                ca = (1-surfAtr.alpha);% absorption coef
                crd = surfAtr.alpha*(1-surfAtr.beta); % diffuse reflection coef
                crs = surfAtr.alpha*surfAtr.beta; % specular reflection coef
                obj.cwforce = facet.fef * obj.intensity/c*obj.affectedArea*( ...
                    (ca+crd)*ray.direction + ...
                    (2*crd/3 + 2*crs*cos(obj.angleOfIncidence)) ...
                    *(-facet.normal));
                obj.cwtorque = cross( obj.hitpos, obj.cwforce);
            else
            
                %%Ref.: Solar Sailing Technology, Dynamics and Mission Applications 1999
                vec_n = -facet.normal;
                vec_t = cross(cross(vec_n,ray.direction),vec_n)/norm(cross(cross(vec_n,ray.direction),vec_n));            
                ang = vectors2angle(vec_n,ray.direction);
                fn = obj.intensity/c*obj.affectedArea*((1+surfAtr.alpha*surfAtr.beta)*cos(ang)^2+...
                    surfAtr.Bf*(1-surfAtr.beta)*surfAtr.alpha*cos(ang)+(1-surfAtr.alpha)*...
                    ((surfAtr.epsf*surfAtr.Bf-surfAtr.epsb*surfAtr.Bb)/(surfAtr.epsf+surfAtr.epsb))*cos(ang)).*vec_n;
                ft = obj.intensity/c*obj.affectedArea*(1-surfAtr.alpha*surfAtr.beta)*cos(ang)*sin(ang).*vec_t;
                obj.cwforce = facet.fef * (fn + ft);
                obj.cwtorque = cross( obj.hitpos, obj.cwforce);
            end
        end
    end
end

