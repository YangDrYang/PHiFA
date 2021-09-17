classdef clCompoundSegment < handle
    %clCompoundSegment segment is part of the target.
    %   It's a compound of facets which share the same surface attributes.
    %   A segment must be a closed surface.
    
    %% PROPERTIES
    properties
        facets clFacet = clTriangle()
        
        %stldata
        stlFile char = ''
        name char = ''
        bInitialiced logical = false
        scale = 1/1000;
        
        %reference data
        biggestArea double = -1 % biggest single facet surface area
        surfaceArea double = -1 % overall surface area
        segmentReferenceCube(3,2) double % holds information about the spatial expanse of this segment
        offset(3,1) double = [0; 0; 0];
        
        %material
        surfAtr clSurfaceAttributes = clSA_const()
        density double = -1 % kg / m^3
        thickness double = -1 % m
        
        %discretization
        resolution double = 1;
        
        %hitdata
        twosided logical = false
        
        %debug
        
    end
    
    %% METHODS
    methods
        function obj = clCompoundSegment()
            %clCompoundSegment Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        %% init
        function init(obj, datastr)
            tic
            disp(['Setting Surface Attributes: ', char(datastr.surfAtr.name) ,'.']);
            obj.surfAtr = datastr.surfAtr;
            obj.offset = datastr.offset;
            obj.name = datastr.name;
            obj.density = datastr.density;
            obj.scale = datastr.scale;
            obj.resolution = datastr.resolution;
            obj.twosided = datastr.twosided;
            obj.thickness = datastr.thickness;
            
            if datastr.bSTL == true
                disp(['Reading ', char(datastr.stlFile), '.']);
                obj.stlFile = datastr.stlFile;
                [f,v,n] = stlread(datastr.stlFile);
                disp(['Time for stl file: ', num2str(toc), ' seconds.'])
                barea = 0;
                disp('Transforming to facets.');
                tic
                obj.surfaceArea = 0;
                nfacets = size(f,1);
                obj.facets(1,nfacets) = clTriangle();
                obj.segmentReferenceCube = [max(v).' min(v).']*obj.scale;
                dispstat('','init');
                for i = 1:nfacets
                    dispstat(sprintf('Transforming %d%%',round(i/nfacets*100)));
                    obj.facets(i) = clTriangle();
                    obj.facets(i).base = v(f(i,1), :).'*obj.scale;
                    obj.facets(i).edge1 = (v(f(i,2), :) - v(f(i,1), :)).'*obj.scale;
                    obj.facets(i).edge2 = (v(f(i,3), :) - v(f(i,1), :)).'*obj.scale;
                    obj.facets(i).normal = n(i,:).'/norm(n(i,:));
                    obj.facets(i).init();
                    obj.surfaceArea = obj.surfaceArea + obj.facets(i).area;
                    if obj.facets(i).area > barea
                        barea = obj.facets(i).area;
                    end
                end
            else
                barea = 0;                
                obj.surfaceArea = 0;
                nfacets = length(obj.facets);
                for i = 1:nfacets
                    obj.surfaceArea = obj.surfaceArea + obj.facets(i).area;
                    if obj.facets(i).area > barea
                        barea = obj.facets(i).area;
                    end
                end
            end
            obj.biggestArea = barea;
            fprintf('Segment has %i facets\n', length(obj.facets));
            fprintf('Object surface area: %.10f m^2\n', obj.surfaceArea);
            fprintf('Mean facet size: %.10f m^2\n', obj.surfaceArea/size(obj.facets,2));            
            fprintf('Boundaries: %.3fm x %.3fm x %.3fm\n', ... 
                obj.segmentReferenceCube(1,1)-obj.segmentReferenceCube(1,2),...
                obj.segmentReferenceCube(2,1)-obj.segmentReferenceCube(2,2),...
                obj.segmentReferenceCube(3,1)-obj.segmentReferenceCube(3,2));
            disp(['Time for facets conversion: ', num2str(toc), ' seconds.'])
            
        end
        
        %% Solid Inertia
        function [volume, mass, angMass, massBarycenter, varargout] = initSolidInertia(obj)
            volume = 0;
            mass = 0;
            angMass = zeros(3);
            fprintf('Calculating Inertia: Assuming object to be solid.\n');
            if obj.density < 0
                fprintf('Warning: Set density before trying to calculate the inertia\n');
                return;
            elseif obj.biggestArea < 0
                fprintf('Warning: Initialice shape before trying to calculate the inertia\n');
                return;
            end
            
            fprintf('Starting discretization of %s segment\n', obj.name);
            tic;
            increment = 1/obj.resolution;
            if obj.getShortestSide()<increment
                fprintf('Warning: Shortest side of object smaller than resolution. Adapting resolution.\n');
                increment = obj.getShortestSide();
                obj.resolution = 1/increment;
            end
            nCounter = 0;
            massBarycenter = [0; 0; 0];
            massCube = increment^3*obj.density;
            angMassCube = diag([1/6*increment^5*obj.density; ...
                1/6*increment^5*obj.density; 1/6*increment^5*obj.density]);        
            
            pt = obj.segmentReferenceCube(1:3,2) + [increment/2; increment/2; increment/2];
            xsteps = ceil( ( obj.segmentReferenceCube(1,1) - obj.segmentReferenceCube(1,2) ) * obj.resolution );
            ysteps = ceil( ( obj.segmentReferenceCube(2,1) - obj.segmentReferenceCube(2,2) ) * obj.resolution );
            zsteps = ceil( ( obj.segmentReferenceCube(3,1) - obj.segmentReferenceCube(3,2) ) * obj.resolution );
            
            ptpc(1:3, 1:xsteps*ysteps*zsteps) = -99;
            nptsinside = 0;
            dispstat('','init');
            maxCounter = xsteps*ysteps*zsteps;
            for i = 0:xsteps-1
                for j = 0:ysteps-1
                    for k = 0:zsteps-1
                        dispstat(sprintf('Processing %d%%',round(nCounter/maxCounter*100)));
                        ipt = pt + [i*increment; j*increment; k*increment];
                        if obj.isPointInside(ipt)
                            nptsinside = nptsinside + 1;
%                             ptpc(:, k*xsteps*ysteps+j*xsteps+i+1) = ipt;
                            ptpc(:, nptsinside) = ipt;
                            volume = volume + increment^3;
                            massBarycenter = massBarycenter + ipt;
                        end
                        nCounter = nCounter + 1;
                    end
                end
            end
            ptpc = ptpc(:,1:nptsinside);
            nElements = length(ptpc);
            mass = volume * obj.density;
            massBarycenter = massBarycenter / nElements;
            fprintf('Local barycenter = [%.4f %.4f %.4f].\n',massBarycenter(1), massBarycenter(2), massBarycenter(3));
            disp(['Split segment into ', num2str(nElements), ' elements.']);
            disp(['Volume: ', num2str(volume), ' m^3.']);
            disp(['Mass: ', num2str(mass), ' kg.']);
            fprintf('Moving facet bases.\n');
            obj.moveBase(massBarycenter);
            fprintf('Calculating moment of inertia.\n');
            
            for i = 1:size(ptpc,2)
%                 if isequal(ptpc(:,i), [ -99 -99 -99 ])
%                     continue;
%                 end
                ipt = ptpc(:,i) - massBarycenter;
                steiner_share = massCube * ...
                                [ipt(2)^2+ipt(3)^2  -ipt(1)*ipt(2)      -ipt(1)*ipt(3); ...
                                -ipt(1)*ipt(2)      ipt(1)^2+ipt(3)^2   -ipt(2)*ipt(3); ...
                                -ipt(1)*ipt(3)      -ipt(2)*ipt(3)      ipt(1)^2+ipt(2)^2];
                angMass = angMass + angMassCube + steiner_share;
%                 angMass = angMass + steiner_share;
                ptpc(:,i) = ipt;
            end
            
            fprintf('Moment of Inertia of Segment:\n');
            fprintf('%.3f\t%.3f\t%.3f\n', angMass(1,1),angMass(1,2),angMass(1,3));
            fprintf('%.3f\t%.3f\t%.3f\n', angMass(2,1),angMass(2,2),angMass(2,3));
            fprintf('%.3f\t%.3f\t%.3f\n', angMass(3,1),angMass(3,2),angMass(3,3));
            
            disp(['Time for inertia calculation: ', num2str(toc), ' seconds.'])
            obj.bInitialiced = true;
            
            if nargout > 4
%                 for i = 1:length(ptpc)
%                     ptpc(:,i) = ptpc(:,i) - massBarycenter;                     
%                 end
                
                if nargout == 5
                    varargout = cell(1);
                elseif nargout == 6
                    varargout = cell(2);
                    varargout{2}(1:length(ptpc)) = increment^3 * obj.density;
                else
                    varargout = cell(nargout-3);
                end
                varargout{1} = ptpc;
            end
        end
        
        %% Hollow Inertia
        function [volume, mass, angMass, massBarycenter, varargout] = initHollowInertia(obj)
            volume = 0;
            mass = 0;
            angMass = zeros(3);
            fprintf('Calculating Inertia: Assuming object to be hollow.\n');
            if obj.density < 0
                fprintf('Warning: Set density before trying to calculate the inertia\n');
                return;
            elseif obj.biggestArea < 0
                fprintf('Warning: Initialice shape before trying to calculate the inertia\n');
                return;
            end
            
            fprintf('Starting discretization of %s segment\n', obj.name);
            tic;
            increment = 1/obj.resolution;
            if increment > obj.thickness
                incrementn = obj.thickness;
                nelementsn = 1;
            else
                nelementsn = round(obj.thickness * obj.resolution);
                incrementn = obj.thickness/nelementsn;
            end
            nCounter = 0;
            massBarycenter = [0; 0; 0];
            massCube = increment^3*obj.density;
            angMassCube = diag([1/6*increment^5*obj.density; ...
                1/6*increment^5*obj.density; 1/6*increment^5*obj.density]);        
            
%             pt = obj.segmentReferenceCube(2,1:3) + [increment/2 increment/2 increment/2];
%             xsteps = size(obj.facets,2);
            ptcl(1:3,1:ceil(obj.surfaceArea/increment^2)*nelementsn) = 0;
            ptm(1,1:ceil(obj.surfaceArea/increment^2)*nelementsn) = 0;
            nptsinside = 0;
            dispstat('','init');
            maxCounter = length(obj.facets);
            for i = 1:maxCounter
                dispstat(sprintf('Processing %d%%',round(nCounter/maxCounter*100)));
                [tmppts, areapp] = obj.facets(i).getPointNet(obj.resolution);
%                 ysteps = length(tmppts);
                for j = 1:size(tmppts,2)
                    stptelem = tmppts(1:3,j) + 0.5*obj.facets(i).normal*incrementn;
                    for k = 1:nelementsn
                        nptsinside = nptsinside + 1;
                        ptelem = stptelem - k*obj.facets(i).normal*incrementn;
                        volumeelem = areapp*incrementn;
                        volume = volume + volumeelem;
                        massBarycenter = massBarycenter + volumeelem*ptelem;
                        ptcl(:,nptsinside) = ptelem;
                        ptm(1,nptsinside) = volumeelem;
                    end
                end
                nCounter = nCounter + 1;
            end
            ptcl = ptcl(:,1:nptsinside);
            ptm = ptm(:,1:nptsinside);
            nElements = length(ptcl);
            mass = volume * obj.density;
            massBarycenter = massBarycenter / volume;
            fprintf('Calculated new mass barycenter.\n');
            disp(['Split segment into ', num2str(nElements), ' elements.']);
            disp(['Volume: ', num2str(volume), ' m^3.']);
            disp(['Mass: ', num2str(mass), ' kg.']);
            fprintf('Moving facet bases.\n');
            obj.moveBase(massBarycenter);
            fprintf('Calculating moment of inertia.\n');
            
            for i = 1:length(ptcl)
                ipt = ptcl(:,i) - massBarycenter;
                steiner_share = ptm(i)*obj.density * ...
                                [ipt(2)^2+ipt(3)^2  -ipt(1)*ipt(2)      -ipt(1)*ipt(3); ...
                                -ipt(1)*ipt(2)      ipt(1)^2+ipt(3)^2   -ipt(2)*ipt(3); ...
                                -ipt(1)*ipt(3)      -ipt(2)*ipt(3)      ipt(1)^2+ipt(2)^2];
                angMass = angMass + angMassCube + steiner_share;
%                 angMass = angMass + steiner_share;

            end
            
            fprintf('Moment of Inertia of Segment:\n');
            fprintf('%.3f\t%.3f\t%.3f\n', angMass(1,1),angMass(1,2),angMass(1,3));
            fprintf('%.3f\t%.3f\t%.3f\n', angMass(2,1),angMass(2,2),angMass(2,3));
            fprintf('%.3f\t%.3f\t%.3f\n', angMass(3,1),angMass(3,2),angMass(3,3));
            
            disp(['Time for inertia calculation: ', num2str(toc), ' seconds.'])
            obj.bInitialiced = true;
            
            if nargout > 4
                for i = 1:length(ptcl)
                    ptcl(:,i) = ptcl(:,i) - massBarycenter;           
                end
                
                if nargout == 5
                    varargout = cell(1);
                elseif nargout == 6
                    varargout = cell(2);
                    varargout{2} = ptm.*obj.density;
                else
                    varargout = cell(nargout-4);
                end
                varargout{1} = ptcl;
            end
        end
        
        %% HIT METHODS
        %% HITS AREA
        function hits = hitsArea(obj, laserstation, ray)
%             hits(1:length(obj.facets)) = clHit();
            hits = [];
%             nHits = 0;
            hits_cell = arrayfun(@(x) obj.hitFacetArea(x,laserstation,ray), obj.facets,...
                    'UniformOutput', false);
            for i = 1:length(hits_cell)
                    hits = cat(2, hits, hits_cell{i});
            end
%             for i = 1:length(obj.facets)
%                 fprintf('Facet %i of %i.\n',i,length(obj.facets));
%                 dotp = dot(obj.facets(i).normal, ray.direction);
%                 if ( dotp<0 )
% %                 if ( dotp<0 && obj.facets(i).inRangeOfBeam(obj, satpos, beam))
%                     nHits = nHits + 1;
%                     hits(nHits).hitpos = obj.facets(i).areabarycenter+obj.offset;
%                     hits(nHits).affectedArea = obj.facets(i).area;
%                     hits(nHits).angleOfIncidence = acos(-dotp/...
%                         (norm(obj.facets(i).normal)*norm(ray.direction)));
%                     hits(nHits) = hits(nHits).init(obj.surfAtr, ...
%                         ray, obj.facets(i), laserstation);
%                 end
%             end
%             hits = hits(1:nHits);
        end
        
        function hit = hitFacetArea(obj, facet, laserstation, ray)
%             fprintf('Facet %i of %i.\n',i,length(obj.facets));
            dotp = dot(facet.normal, ray.direction);
%             crossp = cross(facet.normal, ray.direction);
            if ( dotp<0 )
%                 if ( dotp<0 && obj.facets(i).inRangeOfBeam(obj, satpos, beam))
%                 nHits = nHits + 1;
                hit = clHit();
                hit.hitpos = facet.areabarycenter+obj.offset;
                hit.affectedArea = facet.area;
%                 if crossp(3) == 0
%                     hit.angleOfIncidence = acos(-dotp/...
%                         (norm(facet.normal)*norm(ray.direction)));
%                 elseif crossp(3) > 0
%                     hit.angleOfIncidence = atan2(norm(crossp),dotp);
%                 else
%                     hit.angleOfIncidence = -atan2(norm(crossp),dotp);
%                 end 
%                 hit.angleOfIncidence = acos(-dotp/...
%                     (norm(facet.normal)*norm(ray.direction)));
                hit.angleOfIncidence = vectors2angle(-ray.direction, facet.normal);
                hit = hit.init(obj.surfAtr, ...
                    ray, facet, laserstation);
            else
                hit = [];
            end
        end
        
        function hits = hitsArea_albedo(obj, Y, r_Sun)
            hits = [];
            hits_cell = arrayfun(@(x) obj.hitFacetArea_albedo(x, Y, r_Sun), obj.facets,...
                'UniformOutput', false);
            for i = 1:length(hits_cell)
                hits = cat(1, hits, hits_cell{i});
            end
        end
        
        function hit = hitFacetArea_albedo(obj, facet, Y, r_Sun)
            
            prpgtr = clPropagator.instance();
            albedo = prpgtr.albedo;
            mJD = prpgtr.rso.currentEpoche;
%             Mi2b = Q2Mat(prpgtr.rso.qw(1:4));
            Mi2b = quat2rotm(prpgtr.rso.qw(1:4)');
            mass = prpgtr.rso.mass;
            Mecef2eci = transMatEci2Ecef(mJD)';
            
            %one facet 
            [a_alb,g_alb,a_ir,g_ir,proA_alb,proA_ir] = AlbedoForceTorque(obj, Y(1:3), mJD, r_Sun, ...
                Mecef2eci, Mi2b, 1, facet, mass, albedo, 1, 1, 1, 1);
            hit = {a_alb,g_alb,a_ir,g_ir,proA_alb,proA_ir};

        end
        %% HITS NET
        function hits = hitsNet(obj, laserstation, ray)
            nApproxHits = ceil(obj.surfaceArea*obj.resolution^2);
            hits(1:nApproxHits) = clHit();
            nHits = 0;
%             fprintf('clCompoundSegment.hitsNet()\n');
            for i = 1:length(obj.facets)
%                 fprintf('Facet %i of %i.',i,length(obj.facets));
                dotp = dot(obj.facets(i).normal, ray.direction);
                crossp = cross(obj.facets(i).normal, ray.direction);
                if ( dotp<0 )
%                 if ( dotp<0 && obj.facets(i).inRangeOfBeam(obj, satpos, beam))
                    [ pts, areapp ] = obj.facets(i).getPointNet(obj.resolution);
%                     th(1:length(pts)) = clHit();
                    for j = 1:size(pts,2)
                        nHits=nHits+1;
                        hits(nHits).hitpos = pts(:,j)+obj.offset;
                        hits(nHits).affectedArea = areapp;
                        if crossp(3) == 0
                            hits(nHits).angleOfIncidence = acos(-dotp/...
                                (norm(obj.facets(i).normal)*norm(ray.direction)));
                        elseif crossp(3) > 0
                            hits(nHits).angleOfIncidence = atan2(norm(crossp),dotp);
                        else
                            hits(nHits).angleOfIncidence = -atan2(norm(crossp),dotp);
                        end 

                        hits(nHits) = hits(j).init(obj.surfAtr, ray, ...
                            obj.facets(i), laserstation);
                    end

%                     th(1:end).ecipos = ray.origin + th(1:end).barycenter;
%                     th(1:end).angleOfIncidence = acos(-dotp/...
%                         (norm(obj.facets(i).normal)*norm(ray.direction)));
%                     th(1:end).localFluence = laserstation.getBeamFluence(th(1:end).eciPos);
%                     th(1:end).distFromLaser = norm(th(1:end).eciPos - laserstation.stateVector(1:3));
%                     
%                     targetFluence = th(1:end).localFluence * cos( th(1:end).angleOfIncidence);
%                     th(1:end).lMomAb = - obj.normal * obj.surfAtr.getCouplingCoeffient(targetFluence) * targetFluence * th(1:end).affectedArea;
%                     th(1:end).aMomAb = cross( th(1:end).satPos, th(1:end).linaerMomentum );
%                     th(1:end).lMomPh = 
%                     th(1:end).aMomPh =
%                     hits = [hits, th];
                end
            end
            hits = hits(1:nHits);
        end
        
        %% Point inside
        function yesno = isPointInside(obj, pt) %in local object coord sys
            ray = clRay();
            ray.origin = pt;
            ray.direction = [1; 0; 0];
            satpos = [0; 0; 0];
            yesno = mod(obj.countIntersections(satpos, ray),2);
        end
        
        %% rayIntersections
        function nIntersections = countIntersections(obj, satpos, ray)
            nIntersections = 0;
            for i = 1:size(obj.facets,2)
                if obj.facets(i).rayIntersection(satpos, ray)
                    nIntersections = nIntersections + 1;
                end
            end
        end
        
        %% getIntersections
        function hits = getIntersections(obj, satpos, ray, laserstation)
%             fprintf('clCompoundSegment.getIntersections()\n');

%             hits = [];
%             hits_cell = arrayfun(@(x) obj.getFacetIntersection(x,satpos,ray,laserstation), ...
%                 obj.facets, 'UniformOutput', false);
%             for i = 1:length(hits_cell)
%                     hits = cat(2, hits, hits_cell{i});
%             end
            
            nHits = 0;            
            hits(1:length(obj.facets),1) = clHit();
            for i = 1:length(obj.facets)
                dotp = dot(obj.facets(i).normal, ray.direction);
                crossp = cross(obj.facets(i).normal, ray.direction);
                if(dotp<0)
                    [yesno, hit] = obj.facets(i).getIntersection(satpos, ray);
                    if yesno
                        nHits = nHits + 1;
                        hits(nHits) = hit;
                        hits(nHits).hitpos = hits(nHits).hitpos + obj.offset;
%                         if crossp(3) == 0
%                             hits(nHits).angleOfIncidence = acos(-dotp/...
%                                 (norm(obj.facets(i).normal)*norm(ray.direction)));
%                         elseif crossp(3) > 0
%                             hits(nHits).angleOfIncidence = atan2(norm(crossp),dotp);
%                         else
%                             hits(nHits).angleOfIncidence = -atan2(norm(crossp),dotp);
%                         end
                        hits(nHits).angleOfIncidence = acos(-dotp/...
                            (norm(obj.facets(i).normal)*norm(ray.direction)));                        
                        hits(nHits).affectedArea = 1/laserstation.beamResolution^2;
                        hits(nHits) = hits(nHits).init(obj.surfAtr, ray, ...
                            obj.facets(i), laserstation);
                    end
                end
            end
            if nHits == 0
                hits = [];
            else
                hits = hits(1:nHits);
            end
        end
        
        function hit = getFacetIntersection(obj, facet, satpos, ray, laserstation)
            dotp = dot(facet.normal, ray.direction);
            crossp = cross(facet.normal, ray.direction);
            hit = [];
            if(dotp<0)
                [yesno, hit] = facet.getIntersection(satpos, ray);
                if yesno
%                         nHits = nHits + 1;
%                         hit = hit;
                    hit.hitpos = hit.hitpos + obj.offset;
                    if crossp(3) == 0
                        hit.angleOfIncidence = acos(-dotp/...
                            (norm(facet.normal)*norm(ray.direction)));
                    elseif crossp(3) > 0
                        hit.angleOfIncidence = atan2(norm(crossp),dotp);
                    else
                        hit.angleOfIncidence = -atan2(norm(crossp),dotp);                   
                    end 
                    hit.affectedArea = 1/laserstation.beamResolution^2;
                    hit = hit.init(obj.surfAtr, ray, ...
                        facet, laserstation);
                end
            end
        end
        
        %% getDragImpact
        function hits = getDragIntersections(obj, satpos, ray, dens, resolution, w)
%             hits = [];
%             hits_cell = arrayfun(@(x) obj.getFacetDragIntersection(x,satpos,ray,dens,resolution), ...
%                 obj.facets, 'UniformOutput', false);
%             for i = 1:length(hits_cell)
%                     hits = cat(2, hits, hits_cell{i});
%             end
            
            nHits = 0;
            hits(length(obj.facets)) = clDrag();
            for i = 1:length(obj.facets)
                dotp = dot(obj.facets(i).normal, ray.direction);
                crossp = cross(obj.facets(i).normal, ray.direction);
                if(dotp<0)
                    [yesno, hit] = obj.facets(i).getIntersection(satpos, ray);
                    if yesno
                        nHits = nHits + 1;
                        hits(nHits) = hit2drag(hit);
%                         if crossp(3) == 0
%                             hits(nHits).angleOfIncidence = acos(-dotp/...
%                                 (norm(obj.facets(i).normal)*norm(ray.direction)));
%                         elseif crossp(3) > 0
%                             hits(nHits).angleOfIncidence = atan2(norm(crossp),dotp);
%                         else
%                             hits(nHits).angleOfIncidence = -atan2(norm(crossp),dotp);
%                         end
                        hits(nHits).angleOfIncidence = acos(-dotp/...
                            (norm(obj.facets(i).normal)*norm(ray.direction)));
                        hits(nHits).affectedArea = 1/resolution^2;
                        hits(nHits) = hits(nHits).init(obj.surfAtr, obj, ray, ...
                            obj.facets(i), dens, w);
                    end
                end
            end
            if nHits == 0
                hits = [];
            else
                hits = hits(1:nHits);
            end
        end
        
        function hit = getFacetDragIntersection(obj, facet, satpos, ray, dens, resolution)
            dotp = dot(facet.normal, ray.direction);
            crossp = cross(facet.normal, ray.direction);
            hit = [];
            if(dotp<0)
                [yesno, hit] = facet.getIntersection(satpos, ray);
                if yesno
                    hit = hit2drag(hit);
%                     hit.angleOfIncidence = acos(-dotp/...
%                         (norm(facet.normal)*norm(ray.direction)));
                    if crossp(3) > 0
                        hit.angleOfIncidence = atan2(norm(crossp),dotp);
                    else
                        hit.angleOfIncidence = -atan2(norm(crossp),dotp);
                    end 
                    hit.affectedArea = 1/resolution^2;
                    hit.projectedArea = hit.affectedArea * cos(hit.angleOfIncidence);
                    hit = hit.init(obj.surfAtr, obj, ray, ...
                        facet, dens, w);
                else
                    hit = [];
                end
            end
        end
        
        %% moveBase
        function moveBase(obj, moveBy)
            for i = 1:length(obj.facets)
                obj.facets(i).moveBase(moveBy);
%                 obj.segmentReferenceCube = obj.facets(i).compareReferenceCube(obj.segmentReferenceCube);
            end
            obj.segmentReferenceCube(:,1) = obj.segmentReferenceCube(:,1) - moveBy;
            obj.segmentReferenceCube(:,2) = obj.segmentReferenceCube(:,2) - moveBy;
        end
        
        %% compare Reference Cubes
        function newRefCube = compareRefCube(obj, refCube)
            newRefCube(1:3, 1:2) = 0;
            for i = 1:3
                newRefCube(i,1) = max([ obj.segmentReferenceCube(i,1)+obj.offset(i), refCube(i,1) ]);
                newRefCube(i,2) = min([ obj.segmentReferenceCube(i,2)+obj.offset(i), refCube(i,2) ]);
            end
        end
        
        %% getlongestside
        function length = getLongestSide(obj)
            length = 0;
            for i = 1:3
                tmplength = obj.segmentReferenceCube(i,1)-obj.segmentReferenceCube(i,2);
                if tmplength > length
                    length = tmplength;
                end
            end
        end
        
        %% getshortestside
        function length = getShortestSide(obj)
            length = 10000; % we gettin' problems if there's an object bigger than 10 km #deathstar
            for i = 1:3
                tmplength = obj.segmentReferenceCube(i,1)-obj.segmentReferenceCube(i,2);
                if tmplength < length
                    length = tmplength;
                end
            end
        end
    end
end

