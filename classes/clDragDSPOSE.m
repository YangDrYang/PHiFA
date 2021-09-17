classdef clDragDSPOSE < clDrag
    
    methods
        function obj = init(obj, tar_obj, seg_obj, ray, facet, dens, w, varargin)
            % INPUT
            % w     - satellite rotational speed
            % ray   - ray symbolizing wind speed
            surfAtr = seg_obj.surfAtr;
            
            if nargin < 8
                init@clDrag(obj, tar_obj, surfAtr, ray, facet, dens, w);
            else
                w_rel_i = varargin{1};
                wind_rel_i = varargin{2};
                w_atmos_b = varargin{3};
                wind_rel_b = QForm(tar_obj.qw(1:4), wind_rel_i);

                Cd = tar_obj.Cd;
                obj.las2hitpos = obj.hitpos - ray.origin;
                obj.distFromLaser = norm(obj.las2hitpos);
                normal_i = QTForm(tar_obj.qw(1:4), -facet.normal);
                c_p = facet.areabarycenter+seg_obj.offset;
                c_p_i = QTForm(tar_obj.qw(1:4), c_p);

                %drag
    %             wind_rel_i=-QTForm(tar_obj.qw(1:4),ray.origin);
                velonorm = norm(wind_rel_i);

                drag_force = 0.5*Cd*obj.projectedArea*dens*velonorm^2;
                f1 = drag_force*wind_rel_i/velonorm;
                cxw = cross(c_p_i, w_rel_i);
                f2 = wind_rel_i*dens*obj.affectedArea*dot(normal_i,cxw)*Cd/2;
                f3 = cxw*dens*obj.projectedArea*velonorm*Cd/2;

                cpxwind = cross(c_p, wind_rel_b);
                rpxwind = cross(obj.hitpos, wind_rel_b);
                rpxwrel = cross(obj.hitpos, w_atmos_b);
                m1 = cpxwind*dens*obj.projectedArea*dot(-facet.normal,wind_rel_b)*Cd/2;
                m2 = dot(-facet.normal, rpxwrel)*...
                    rpxwind*obj.affectedArea*Cd*dens/2;
                m3 = obj.projectedArea*velonorm*cross(obj.hitpos,rpxwrel)*Cd*dens/2;


                obj.force = f1 + f2 + f3;
                obj.force = QForm(tar_obj.qw(1:4), obj.force);
                obj.torque = m1 + m2 + m3;
            
            end
        end
    end
end