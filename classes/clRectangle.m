classdef clRectangle < clFacet
    %clRectangle Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        edge1(3,1) double = [0; 0; 1]
        edge2(3,1) double = [0; 1; 0]
    end
    
    methods
        function obj = clRectangle()
            %clRectangle Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function init(obj)
            obj.area = norm(cross(obj.edge1, obj.edge2));
            obj.areabarycenter = 0.5 .* obj.edge1 + 0.5 .* obj.edge2 + obj.base;
            obj.normal = obj.normal/norm(obj.normal);
        end
        
        function yesno = rayIntersection(obj, satpos, ray)
            if dot(ray.direction, obj.normal)~=0
                A = [obj.edge1, obj.edge2, -ray.direction];
                if rcond(A) < 1e-20
                    sv = pinv(A) * ( ray.origin - ( satpos + obj.base ));
                else
                    sv = A \ ( ray.origin - ( satpos + obj.base ));
                end
                yesno = (sv(3)>=0) && (sv(1)>=0) && ...
                    (sv(1)<=1) && (sv(2)<=1) && (sv(2)>=0);
            else
                yesno = false;
            end
        end
        
        function [yesno, hit] = getIntersection(obj, satpos, ray)
            if dot(ray.direction, obj.normal)~=0
                A = [obj.edge1, obj.edge2, -ray.direction];
                if rcond(A) < 1e-20
                    sv = pinv(A) * ( ray.origin - ( satpos + obj.base ));
                else
                    sv = A \ ( ray.origin - ( satpos + obj.base ));
                end
                yesno = (sv(3)>=0) && (sv(1)>=0) && ...
                    (sv(1)<=1) && (sv(2)<=1) && (sv(2)>=0);
            else
                yesno = false;
            end
            if yesno
                hit = clHit();
                hit.hitpos = obj.base + sv(1)*obj.edge1 + sv(2)*obj.edge2;
                hit.las2hitpos = ray.origin + sv(3)*ray.direction;
                hit.distFromLaser = sv(3);
            else
                hit = [];
            end
        end
        
        function refcube = getReferenceCube(obj)
            refcube(2,3) = [obj.base obj.base];
            for i = 1:3
                if obj.base(i) + obj.edge1(i) > refcube(i,1)
                    refcube(i,1) = obj.base(i) + obj.edge1(i);
                elseif obj.base(i) + obj.edge1(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i) + obj.edge1(i);
                end
                
                if obj.base(i) + obj.edge2(i) > refcube(i,1)
                    refcube(i,1) = obj.base(i) + obj.edge2(i);
                elseif obj.base(i) + obj.edge2(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i) + obj.edge2(i);
                end
                
                if obj.base(i) + obj.edge1(i) + obj.edge2(i) > refcube(i,1)
                    refcube(i,1) = obj.base(i) + obj.edge1(i) + obj.edge2(i);
                elseif obj.base(i) + obj.edge1(i) + obj.edge2(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i) + obj.edge1(i) + obj.edge2(i);
                end
            end
        end
        
        function refcube = compareReferenceCube(obj, refcube)
            for i = 1:3 
                if obj.base(i)> refcube(i,1)
                    refcube(i,1) = obj.base(i);
                elseif obj.base(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i);
                end
                
                if obj.base(i) + obj.edge1(i) > refcube(i,1)
                    refcube(i,1) = obj.base(i) + obj.edge1(i);
                elseif obj.base(i) + obj.edge1(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i) + obj.edge1(i);
                end
                
                if obj.base(i) + obj.edge2(i) > refcube(i,1)
                    refcube(i,1) = obj.base(i) + obj.edge2(i);
                elseif obj.base(i) + obj.edge2(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i) + obj.edge2(i);
                end
                
                if obj.base(i) + obj.edge1(i) + obj.edge2(i) > refcube(i,1)
                    refcube(i,1) = obj.base(i) + obj.edge1(i) + obj.edge2(i);
                elseif obj.base(i) + obj.edge1(i) + obj.edge2(i) < refcube(i,2)
                    refcube(i,2) = obj.base(i) + obj.edge1(i) + obj.edge2(i);
                end
            end
        end
        
        % return net of points lying on the surface of the facet. if
        % 1/resolution is higher than facet size, barycenter and facets
        % area are returned
        function [ pts, areapp ] = getPointNet(obj, resolution)
            nelements1 = ceil( norm(obj.edge1) * resolution );
            increment1 = norm(obj.edge1) / nelements1;
            nelements2 = ceil( norm(obj.edge2) * resolution );
            increment2 = norm(obj.edge2) / nelements2;
            
            pts(1:3,1:nelements1*nelements2) = 0;
            npts = 0;
            incedge1 = increment1*obj.edge1/norm(obj.edge1);
            incedge2 = increment2*obj.edge2/norm(obj.edge2);
            stpt = obj.base - 0.5*(incedge1 + incedge2);
            for i = 1:nelements1
                for j = 1:nelements2
                    npts = npts + 1;
                    pts(:,npts) = stpt + i*incedge1 + j*incedge2;
                end
            end
            pts = pts(:,1:npts);
            areapp = increment1*increment2;
                
            if size(pts, 1)<=1
                pts(1:3,1) = obj.areabarycenter;
                areapp = obj.area;
            end
        end

        %% keep it for later
%         function hit = rayHit(obj, satpos, ray)
%             A = [obj.edge1 obj.edge2 ray.direction];
%             sv = ( ray.origin - ( satpos(1:3) + obj.base ))/A;
%             if (sv(1)<=1) && (sv(1)>=0) && (sv(2)<=1) && (sv(2)>=0)
%                 dotp = dot(obj.normal, laserStation.direction);
%                 hit.satpos = obj.barycenter;
%                 hit.ecipos = satpos(1:3) + obj.barycenter;
%                 hit.angleOfIncidence = acos(-dotp);
% %                 hit.localFluence = laserStation.getBeamFluence(hit.eciPos);
%                 hit.distFromLaser = norm(hit.eciPos - ray.origin);
%                 hit.affectedArea = obj.area;
% %                 targetFluence = hit.localFluence * cos( hit.angleOfIncidence);
% %                 hit.linearMomentum = - obj.normal * targetObject.material.getCouplingCoeffient(targetFluence) * targetFluence * hit.affectedArea;
% %                 hit.angularMomentum = cross( hit.satPos, hit.linaerMomentum );
%             end
%         end
%         
%         function hit = beamHit(obj, satpos, beam)
%             dotp = dot(obj.normal, laserStation.direction);
%             hit = clHit();
%             if ( dotp<0 ) && (point_to_line(satpos(1:3) + obj.barycenter, beam.origin, beam.origion + beam.direction) < beam.referenceRadius)
%                 hit.satpos = obj.barycenter;
%                 hit.ecipos = satpos(1:3) + obj.barycenter;
%                 hit.angleOfIncidence = acos(-dotp);
% %                 hit.localFluence = laserStation.getBeamFluence(hit.eciPos);
%                 hit.distFromLaser = norm(hit.eciPos - laserStation.stateVector(1:3));
%                 hit.affectedArea = obj.area;
% %                 targetFluence = hit.localFluence * cos( hit.angleOfIncidence);
% %                 hit.linearMomentum = - obj.normal * targetObject.material.getCouplingCoeffient(targetFluence) * targetFluence * hit.affectedArea;
% %                 hit.angularMomentum = cross( hit.satPos, hit.linaerMomentum );
%             end
%         end
    end
end

