classdef clCompoundTarget < clTarget
    %clCompoundTarget Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % segments
        segments clCompoundSegment = clCompoundSegment()
        nFacets double = 0
        % parameters concerning whole object
        referenceCube(3,2) double % holds information about the spatial expanse of the entire object
        surfaceArea double = 0;
        Lc double = 0 % characteristic length
        % method to calculate momentum
        hitmethod eHitMethod = eHitMethod.Beam
        % miscellaneous
        bInitialiced logical = false
        % coupled orbit/attitude dynamics
        points % point array
        pointsmass double = -1 % mass of each point
    end
    
    methods
        function obj = clCompoundTarget()
            %clCompoundTarget Construct an instance of this class
            %   Detailed explanation goes here
            obj.shape = eShape.discrete;
        end
        
        %% init
        function init(obj, seg_input, tar_input)
            fprintf('#######Initialicing target %s .#######\n', obj.name);
            
            if tar_input.bUseProvidedParam
                obj.initWithProvidedData(seg_input, tar_input);
                return;
            end
            
            nsegments = length(seg_input);
            obj.segments(1,nsegments) = clCompoundSegment();
            tmpvol(1:nsegments) = 0;
            tmpmas(1:nsegments) = 0;
            tmpbc(1:3,1:nsegments) = 0;
            tmpmoi = zeros(3);
            tmpmoi(:,:,nsegments) = zeros(3);
            obj.points = [];
            obj.pointsmass = [];
%             tmpmoi(1:3,1:3,1:nsegments) = zeros(3);
            for i = 1:nsegments
                fprintf('#######Initialicing segment %i of %i .#######\n', i, nsegments);
                if seg_input(i).bSTL == true
                    obj.segments(i) = clCompoundSegment();
                else
                    obj.segments(i) = seg_input(i).segobj;
                end
                obj.segments(i).init(seg_input(i));
                
                if seg_input(i).solid == true
                    [ tmpvol(i), tmpmas(i), tmpmoi(:,:,i), tmpbc(:,i), ...
                        tmppoints, tmppointmass ] = obj.segments(i).initSolidInertia();
                else
                    [ tmpvol(i), tmpmas(i), tmpmoi(:,:,i), tmpbc(:,i), ...
                        tmppoints, tmppointmass ] = obj.segments(i).initHollowInertia();
                end
                obj.points = cat(2,obj.points,tmppoints);
                obj.pointsmass = cat(2,obj.pointsmass,tmppointmass);
                obj.nFacets = obj.nFacets + length(obj.segments(i).facets);
            end
            fprintf('####### Segments Initialicing done. #######\n');
            fprintf('Calculating barycenter of target.\n');
            %calculate barycenter
            obj.mass = sum(tmpmas);
            barycenter = [0; 0; 0];
            obj.volume = 0;
            for i = 1:nsegments
                barycenter = barycenter + tmpmas(i)*(obj.segments(i).offset+tmpbc(:,i));
            end
            barycenter = barycenter / obj.mass;
            fprintf('Target barycenter before applying offset \n\t= [%.4f %.4f %.4f].\n' ...
                ,barycenter(1), barycenter(2), barycenter(3));
            fprintf('Calculating moment of inertia of target\n');
            for i = 1:nsegments
                obj.segments(i).offset = obj.segments(i).offset + tmpbc(:,i) - barycenter;
%                 obj.segments(i).offset = obj.segments(i).offset - barycenter;
                obj.volume = obj.volume + tmpvol(i);
%                 obj.mass = obj.mass + tmpmas(i);
                os = obj.segments(i).offset;
                obj.moi = obj.moi + tmpmoi(:,:,i) + tmpmas(i) * ...
                                [os(2)^2+os(3)^2  -os(1)*os(2)      -os(1)*os(3); ...
                                -os(1)*os(2)      os(1)^2+os(3)^2   -os(2)*os(3); ...
                                -os(1)*os(3)      -os(2)*os(3)      os(1)^2+os(2)^2];
                
                if i == 1
                    obj.referenceCube = obj.segments(i).segmentReferenceCube+...
                        [obj.segments(i).offset obj.segments(i).offset];
                else
                    obj.referenceCube = obj.segments(i).compareRefCube(obj.referenceCube);
                end
                
                obj.surfaceArea = obj.surfaceArea + obj.segments(i).surfaceArea;
            end
            obj.invmoi = inv(obj.moi);
            
            fprintf('Moment of Inertia:\n');
            fprintf('%.3f\t%.3f\t%.3f\n', obj.moi(1,1),obj.moi(1,2),obj.moi(1,3));
            fprintf('%.3f\t%.3f\t%.3f\n', obj.moi(2,1),obj.moi(2,2),obj.moi(2,3));
            fprintf('%.3f\t%.3f\t%.3f\n', obj.moi(3,1),obj.moi(3,2),obj.moi(3,3));
            
            tmp = obj.referenceCube(:,1) - obj.referenceCube(:,2);
            obj.Lc = sqrt(sum(tmp.^2));
            obj.area_drag = tmp(1)*tmp(2)*tmp(3)/min(tmp);
            obj.area_solar = obj.area_drag;
            obj.bInitialiced = true;
            if obj.volume < 1
                disp_vol = [ num2str(obj.volume*1E6) ' cm^3' ];
            else
                disp_vol = [ num2str(obj.volume) ' m^3' ];
            end
            
            fprintf('Volume of mass elements: %s\nMass: %.6f kg\n',disp_vol, obj.mass);
        end
        
        function initWithProvidedData(obj, seg_input, tar_input)
            nsegments = length(seg_input);
            obj.segments(1,nsegments) = clCompoundSegment();
            for i = 1:nsegments
                fprintf('#######Initialicing segment %i of %i .#######\n', i, nsegments);
                if seg_input(i).bSTL == true
                    obj.segments(i) = clCompoundSegment();
                else
                    obj.segments(i) = seg_input(i).segobj;
                end
                obj.segments(i).init(seg_input(i));
                
                obj.nFacets = obj.nFacets + length(obj.segments(i).facets);
                if i == 1
                    obj.referenceCube = obj.segments(i).segmentReferenceCube+...
                        [obj.segments(i).offset obj.segments(i).offset];
                else
                    obj.referenceCube = obj.segments(i).compareRefCube(obj.referenceCube);
                end
                obj.surfaceArea = obj.surfaceArea + obj.segments(i).surfaceArea;
                obj.segments(i).moveBase(tar_input.barycenter);
            end
            
            obj.mass = tar_input.mass;
            obj.moi = tar_input.inertia;
            obj.invmoi = inv(obj.moi);
            obj.mtm = tar_input.magnetic_tensor;
            tmp = obj.referenceCube(:,1) - obj.referenceCube(:,2);
            obj.Lc = sqrt(sum(tmp.^2));
            obj.area_drag = tar_input.area_drag;
            obj.area_solar = tar_input.area_solar;
            obj.Cr = tar_input.Cr;
            obj.Cd = tar_input.Cd;
            obj.bInitialiced = true;
        end
        
        %% getting accelerations
        %% ablation and \gls{cw}
        function [impulse, momentum, cwforce, cwtorque, projectedArea, varargout] = ...
                getLaserAccel(obj, laserstation)
            ray = laserstation.raySAT;
            % get hits
            nfacets = 0;
            hits = [];
            impulse = [0; 0; 0];
            momentum = [0; 0; 0];
            cwforce = [0; 0; 0];
            cwtorque = [0; 0; 0];
            projectedArea = 0;
            for i = 1:length(obj.segments) 
                nfacets = nfacets + length(obj.segments(i).facets); 
            end
            if obj.hitmethod == eHitMethod.Area
                i_seg = 1:length(obj.segments);
                hits_cell = arrayfun(@(x) obj.segments(x).hitsArea(laserstation, ray) ,i_seg ...
                    ,'UniformOutput', false);
%                 for i = 1:length(hits_cell)
%                     hits = cat(1, hits, hits_cell{i});
%                 end
%                 for i = 1:length(obj.segments)
%                     fprintf('Segment %s\n', obj.segments.name);
%                     tmphits = obj.segments(i).hitsArea(laserstation, ray);
%                     hits = cat(1, hits, tmphits);
%                 end
            elseif obj.hitmethod == eHitMethod.Net
                i_seg = 1:length(obj.segments);
                hits_cell = arrayfun(@(x) obj.segments(x).hitsNet(laserstation, ray) ,i_seg ...
                    ,'UniformOutput', false);
%                 for i = 1:length(hits_cell)
%                     hits = cat(1, hits, hits_cell{i});
%                 end
%                 for i = 1:length(obj.segments)
%                     fprintf('Segment %s\n', obj.segments.name);
%                     tmphits = obj.segments(i).hitsNet(laserstation, ray);
%                     hits = cat(1, hits, tmphits);
%                 end
            elseif obj.hitmethod == eHitMethod.Beam
                % discretizise beam
                [rays, ~] = discretizeBeam(ray,obj.Lc,laserstation.beamResolution);
                % check for each beam how often he hits and choose first
                nSeg = length(obj.segments);
                % (START PARFORLOOP)
                hits_cell = cell(1,length(rays));
                parfor i = 1:length(rays)
%                 for i = 1:length(rays)
                    loc_obj = obj;
                    seghits = [];
                    for j = 1:nSeg
                        tmphits = loc_obj.segments(j).getIntersections(loc_obj.segments(j).offset,...
                            rays(i), laserstation);
                        seghits = cat(1,seghits,tmphits);
                    end
                    if ~isempty(seghits)
                        hits_cell{i} = closestHit(seghits);
                    end
                end
                empt = find(cellfun('isempty', hits_cell));
                hits_cell(empt) = [];
                % (END PARFORLOOP) 
            end
            % hits to momentum and output
            if nargout > 5
                varargout = cell(nargout-4,1);
                nhits_ue0 = 0;
                marea = 0; %affected area
                mcm = 0; %coupling coefficient
                maoi = 0; %angle of incidence
                mint = 0; %intensity
                for i = 1:length(hits)
                    impulse = impulse + hits(i).impuls;
                    momentum = momentum + hits(i).moment;
                    cwforce = cwforce + hits(i).cwforce;
                    cwtorque = cwtorque + hits(i).cwtorque;
                    marea = marea + hits(i).affectedArea;
                    mcm = mcm + hits(i).couplingCoef;
                    maoi = maoi + hits(i).angleOfIncidence;
                    mint = mint + hits(i).intensity;
                    if norm(hits(i).impuls)>0
                        nhits_ue0 = nhits_ue0+1;
                    end
                end
                mimp = norm(impulse)/length(hits);
                mamom = norm(momentum)/length(hits);
                marea = marea/length(hits);
                mcm = mcm/length(hits);
                maoi = maoi/length(hits);
                mint = mint/length(hits);
                
                fprintf('Target: %i hits\n', nhits_ue0);
            else
                for i = 1:length(hits_cell)
                    if ~isempty(hits_cell{i})
                        if length(hits_cell{i})>1
                            for j = 1:length(hits_cell{i})
                                impulse = impulse + hits_cell{i}(j).impuls;
                                momentum = momentum + hits_cell{i}(j).moment;
                                cwforce = cwforce + hits_cell{i}(j).cwforce;
                                cwtorque = cwtorque + hits_cell{i}(j).cwtorque;
                                if laserstation.lsType == 2
                                    if laserstation.bPulsed
                                        if norm(hits_cell{i}(j).impuls)&& norm(hits_cell{i}(j).moment)
                                            projectedArea = projectedArea + hits_cell{i}(j).projectedArea;
                                        end
                                    else
                                        if norm(hits_cell{i}(j).cwforce)&& norm(hits_cell{i}(j).cwtorque)
                                            projectedArea = projectedArea + hits_cell{i}(j).projectedArea;
                                        end
                                    end
                                else
                                    if norm(hits_cell{i}(j).cwforce)&& norm(hits_cell{i}(j).cwtorque)
                                        projectedArea = projectedArea + hits_cell{i}(j).projectedArea;
                                    end
                                end
                            end
                        else
                            impulse = impulse + hits_cell{i}.impuls;
                            momentum = momentum + hits_cell{i}.moment;
                            cwforce = cwforce + hits_cell{i}.cwforce;
                            cwtorque = cwtorque + hits_cell{i}.cwtorque;
                            if laserstation.lsType == 2
                                if laserstation.bPulsed
                                    if norm(hits_cell{i}.impuls)&& norm(hits_cell{i}.moment)
                                        projectedArea = projectedArea + hits_cell{i}.projectedArea;
                                    end
                                else
                                    if norm(hits_cell{i}.cwforce)&& norm(hits_cell{i}.cwtorque)
                                        projectedArea = projectedArea + hits_cell{i}.projectedArea;
                                    end
                                end
                            else
                                if norm(hits_cell{i}.cwforce)&& norm(hits_cell{i}.cwtorque)
                                    projectedArea = projectedArea + hits_cell{i}.projectedArea;
                                end
                            end
                        end
                    end
                end    
            end
            % Turn forces into ECI
            cwforce = QTForm(obj.qw(1:4), cwforce);
            impulse = QTForm(obj.qw(1:4), impulse);
            
            normimp = norm(impulse);
            
            ang = acos(dot(impulse, ray.direction)/norm(ray.direction)/normimp);
            raimp = normimp*sin(ang);
            aximp = normimp*cos(ang);
            fom = raimp/aximp;
            if nargout > 5
                varargout{1} = struct('mimp', mimp, 'mamom', mamom, ...
                        'marea', marea, 'mcm', mcm, 'maoi', maoi, ...
                        'mint', mint, 'nhits', nhits_ue0, ...
                        'fom', fom);
            end
            if nargout > 6
                varargout{2} = hits;
            end
        end
        
        %% gravitation
        function [potential,force, torque] = getGravityAccel(obj,r_Sun,r_Moon,E,UT1_UTC,TT_UTC,x_pole,y_pole)
%             C_i2b = Q2Mat(obj.qw(1:4));
             C_i2b = quat2rotm(obj.qw(1:4)');
%             [a_gravity_inertial, g_gravity_body] = ...
%                 AccelGravityGradientTorque(obj.currentEpoche,r_Sun,r_Moon,obj.xv(1:3),E,...
%                 UT1_UTC,TT_UTC,x_pole,y_pole, C_i2b, obj.moi);
            [gravity_potential, a_gravity_inertial, g_gravity_body] = ...
                Gravity_Potential_Force_Torque(obj.currentEpoche,r_Sun,r_Moon,obj.xv(1:3),E,...
                UT1_UTC,TT_UTC,x_pole,y_pole, C_i2b, obj.moi, obj.mass);
            
            potential = gravity_potential;
            force = a_gravity_inertial.*obj.mass;
            torque = g_gravity_body;
        end
        
        %% drag
        function [force, torque, projectedArea] = getDragAccel(obj, dens, T, Omega)
%             tmpmthd = obj.hitmethod;
%             obj.hitmethod = eHitMethod.Beam;
            p = clPropagator.instance();
            if p.AuxParam.dof == 3
                if obj.area_drag == 0
                    obj.area_drag = obj.Lc^2*pi/4;
                end
                force = obj.mass*AccelDrag(dens,obj.xv(1:3),obj.xv(4:6),T,obj.area_drag,obj.mass,obj.Cd,Omega);
                torque = [0;0;0];
                projectedArea = obj.area_drag;
                return;
            end
            xv_ecef = ECI2ECEF(obj.currentEpoche, obj.xv');
            wind_i = wind_def(obj.currentEpoche, xv_ecef);
            v_rel_i = obj.xv(4:6) - wind_i;
            v_abs = norm(v_rel_i);
            % construct ray from satellite velocity
            wind_body = clRay();
            wind_body.origin = v_rel_i;
            wind_body.direction = -wind_body.origin/v_abs;
            wind_body.origin = QForm(obj.qw(1:4), wind_body.origin);
            wind_body.direction = QForm(obj.qw(1:4), wind_body.direction);
            wind_body.xdir = cross(wind_body.direction, obj.xv(1:3));
            wind_body.ydir = cross(wind_body.direction, wind_body.xdir);
            wind_body.xdir = wind_body.xdir/norm(wind_body.xdir);
            wind_body.ydir = wind_body.ydir/norm(wind_body.ydir);
            % discretizise beam
            resolution = p.station.beamResolution;
            if obj.hitmethod == eHitMethod.Beam
                [rays, ~] = discretizeBeam(wind_body,obj.Lc/2,p.station.beamResolution);
                % check for each beam how often he hits and choose first;
                drag_cell = cell(1,length(rays));
                parfor i = 1:length(rays)
                    loc_obj = obj;
                    seghits = [];
                    for j = 1:length(obj.segments)
                        tmphits = loc_obj.segments(j).getDragIntersections(obj.segments(j).offset, rays(i), ...
                            dens, resolution, obj.qw(5:7));
                        seghits = cat(1,seghits,tmphits);
                    end
                    if ~isempty(seghits)
                        drag_cell{i} = closestHit(seghits);
    %                     nHits = nHits+1;
    %                     beamhits(nHits) = closestHit(seghits);
                    end
                end
                empt = find(cellfun('isempty', drag_cell));
                drag_cell(empt) = [];
            else
                drag_cell = cell(1, obj.nFacets);
                n = 0;
                
                w_atmos_i = [0;0;7.29211514670698*10^-5];
                w_atmos_b = QForm(obj.qw(1:4), w_atmos_i);
                wind_rel_i =  wind_i - obj.xv(4:6);
                wind_rel_b = QForm(obj.qw(1:4), wind_rel_i);
                w_sat_i = QForm(obj.qw(1:4), obj.qw(5:7));
                w_rel_i = w_sat_i - w_atmos_i;
                w_rel_b = obj.qw(5:7) - w_atmos_b; % needed for torque
                
                drag = clDrag();
                for i = 1:length(obj.segments)
                    for j = 1:length(obj.segments(i).facets)
                        dotp = dot(obj.segments(i).facets(j).normal, wind_body.direction);
                        if ( dotp<0 )
                            n=n+1;
                            drag_cell{n} = clDrag();
                            drag_cell{n}.affectedArea = 0;
                            drag_cell{n}.projectedArea = 0;
                            drag_cell{n}.hitpos = [0;0;0];
                            drag_cell{n}.angleOfIncidence = 0;
                            drag_cell{n}.distFromLaser = 0;
                            drag_cell{n}.las2hitpos = [0;0;0];
                            if obj.hitmethod == eHitMethod.Net
                                [ pts, areapp ] = obj.segments(i).facets(j).getPointNet(resolution);
                            else
                                pts = obj.segments(i).facets(j).areabarycenter;
                                areapp = obj.segments(i).facets(j).area;
                            end
                            npts = size(pts,2);
                            for k= 1:npts
                                drag.hitpos = pts(:,k) + obj.segments(i).offset;
                                drag.angleOfIncidence = acos(-dotp/...
                                    (norm(obj.segments(i).facets(j).normal)*norm(wind_body.direction)));
                                drag.affectedArea = areapp;
                                drag.projectedArea = areapp*cos(drag.angleOfIncidence);

                                drag = drag.init(obj, obj.segments(i), wind_body, ...
                                    obj.segments(i).facets(j), dens, obj.qw(5:7), w_rel_i, ...
                                    wind_rel_i, w_rel_b);

%                                 drag_cell{n} = drag_cell{n} + drag;
                                drag_cell{n}.force = drag_cell{n}.force + drag.force;
                                drag_cell{n}.torque = drag_cell{n}.torque + drag.torque;
                                drag_cell{n}.affectedArea = drag_cell{n}.affectedArea + drag.affectedArea;
                                if norm(drag_cell{n}.force) && norm(drag_cell{n}.torque) 
                                    drag_cell{n}.projectedArea = drag_cell{n}.projectedArea + drag.projectedArea;
                                end
                                drag_cell{n}.hitpos = drag_cell{n}.hitpos + drag.hitpos;
                                drag_cell{n}.angleOfIncidence = drag_cell{n}.angleOfIncidence + drag.angleOfIncidence;
                                drag_cell{n}.distFromLaser = drag_cell{n}.distFromLaser + drag.distFromLaser;
                                drag_cell{n}.las2hitpos = drag_cell{n}.las2hitpos + drag.las2hitpos;
                            end
                            drag_cell{n}.hitpos = drag_cell{n}.hitpos/npts;
                            drag_cell{n}.angleOfIncidence = drag_cell{n}.angleOfIncidence/npts;
                            drag_cell{n}.distFromLaser = drag_cell{n}.distFromLaser/npts;
                            drag_cell{n}.las2hitpos = drag_cell{n}.las2hitpos/npts;
                        end
                    end
                end
                empt = find(cellfun('isempty', drag_cell));
                drag_cell(empt) = [];
            end
%             hits = beamhits(1:nHits);
            % calculate force and torque from found hits
            force = [0; 0; 0];
            torque = [0; 0; 0];
            projectedArea = 0;%area to mass ratio
            for i = 1:length(drag_cell)
                force = force + drag_cell{i}.force;
                torque = torque + drag_cell{i}.torque;
                if norm(drag_cell{i}.force) && norm(drag_cell{i}.torque)
                    projectedArea = projectedArea + drag_cell{i}.projectedArea;%added by Yang
                end
            end
            % turn into eci coord sys
            force = QTForm(obj.qw(1:4), force);
            
%             obj.hitmethod = tmpmthd;
        end
        
        function [force, torque, projectedArea] = getDragGradientAccel(obj, T, Omega)
%             tmpmthd = obj.hitmethod;
%             obj.hitmethod = eHitMethod.Beam;
            p = clPropagator.instance();
            
            Oplus_flag = 1;        

            [yy, mon, day, hour, minute, sec] = invjday(obj.currentEpoche+2400000.5);
            days = finddays(yy, mon, day, hour, minute, sec);
            idoy = floor(days);
            ut = hour*3600+minute*60+sec; % seconds in day (UT)

            i = find((yy==p.swdata(1,:)) & (mon==p.swdata(2,:)) & (day==p.swdata(3,:)),1,'first');
            sw = p.swdata(:,i);        
            ap(1) = sw(23); % Arithmetic average of the 8 Ap indices for the day
            ap(2) = sw(15); % 3 hr AP index for current time
            sw_1 = p.swdata(:,i-1);
            % Define Solar Flux Values
            if floor(hour/3) == 0
                ap(3) = sw_1(22); % 3 hr AP index for 3 hrs before current time
            else
                ind = 15+floor(hour/3)-1;
                ap(3) = sw(ind); % 3 hr AP index for 3 hrs before current time
            end

            if floor(hour/6) == 0
                if floor(hour/3) == 0
                    ap(4) = sw_1(21); % 3 hr AP index for 6 hrs before current time
                else    
                    ap(4) = sw_1(22); % 3 hr AP index for 6 hrs before current time
                end
            else
                ind = 15+floor(hour/3)-2;
                ap(4) = sw(ind); % 3 hr AP index for 6 hrs before current time
            end

            if floor(hour/9) == 0
                if floor(hour/6) == 0   
                    if floor(hour/3) == 0
                        ap(5) = sw_1(20); % 3 hr AP index for 9 hrs before current time
                    else
                        ap(5) = sw_1(21); % 3 hr AP index for 9 hrs before current time
                    end
                else
                    ap(5) = sw_1(22); % 3 hr AP index for 9 hrs before current time
                end
            else
                ind = 15+floor(hour/3)-3;
                ap(5) = sw(ind); % 3 hr AP index for 9 hrs before current time
            end

            sw_2 = p.swdata(:,i-2);

            % Average of eight 3 hr AP indicies from 12 to 33
                                     % hrs prior to current time
            if floor(hour/12) == 0
                if floor(hour/9) == 0   
                    if floor(hour/6) == 0
                        if floor(hour/3) == 0
                            ap(6) = (sw_1(19)+sw_1(18)+sw_1(17)+sw_1(16)...
                                +sw_1(15)+sw_2(22)+sw_2(21)+sw_2(20))/8.0; 
                        else
                            ap(6) = (sw_1(20)+sw_1(19)+sw_1(18)+sw_1(17)...
                                +sw_1(16)+sw_1(15)+sw_2(22)+sw_2(21))/8.0;                         
                        end
                    else
                        ap(6) = (sw_1(21)+sw_1(20)+sw_1(19)+sw_1(18)...
                            +sw_1(17)+sw_1(16)+sw_1(15)+sw_2(22))/8.0;
                    end
                else
                    ap(6) = (sw_1(22)+sw_1(21)+sw_1(20)+sw_1(19)...
                         +sw_1(18)+sw_1(17)+sw_1(16)+sw_1(15))/8.0;
                end
            else
                ind = 15+floor(hour/3)-4;
                ap(6) = (sw(ind)+sw_1(22)+sw_1(21)+sw_1(20)...
                     +sw_1(19)+sw_1(18)+sw_1(17)+sw_1(16))/8.0;
            end

            sw_3 = p.swdata(:,i-3);
            % Average of eight 3 hr AP indicies from 36 to 57
                                     % hrs prior to current time 
            if floor(hour/12) == 0
                if floor(hour/9) == 0   
                    if floor(hour/6) == 0
                        if floor(hour/3) == 0
                            ap(7) = (sw_2(19)+sw_2(18)+sw_2(17)+sw_2(16)...
                                +sw_2(15)+sw_3(22)+sw_3(21)+sw_3(20))/8.0; 
                        else
                            ap(7) = (sw_2(20)+sw_2(19)+sw_2(18)+sw_2(17)...
                                +sw_2(16)+sw_2(15)+sw_3(22)+sw_3(21))/8.0;                         
                        end
                    else
                        ap(7) = (sw_2(21)+sw_2(20)+sw_2(19)+sw_2(18)...
                            +sw_2(17)+sw_2(16)+sw_2(15)+sw_3(22))/8.0;
                    end
                else
                    ap(7) = (sw_2(22)+sw_2(21)+sw_2(20)+sw_2(19)...
                         +sw_2(18)+sw_2(17)+sw_2(16)+sw_2(15))/8.0;
                end
            else
                ind = 15+floor(hour/3)-4;
                ap(7) = (sw_1(ind)+sw_2(22)+sw_2(21)+sw_2(20)...
                     +sw_2(19)+sw_2(18)+sw_2(17)+sw_2(16))/8.0;
            end

            f107 = sw(31);     % observed solar radio noise flux (jansky)
            f107a = sw(32);    % observed 81-day average F10 (jansky)

            flag = ones(1,23);
            flag(9) = -1;
           
            
            if p.AuxParam.dof == 3
                if obj.area_drag == 0
                    obj.area_drag = obj.Lc^2*pi/4;
                end
                
                xv_ecef = ECI2ECEF(obj.currentEpoche, obj.xv');
                [lon, lat, height] = Geodetic(xv_ecef(1:3));
                alt = height/1000;
                xlat = lat*p.const.Deg;
                xlong = lon*p.const.Deg;
                
                [den,~] = nrlmsise00_mex(int32(idoy),ut,alt,xlat,xlong,f107a,f107,ap,flag);
                if Oplus_flag
                    dens = 1.66e-27*(16.0*den(9)) + den(6);
                end 
                force = obj.mass*AccelDrag(dens,obj.xv(1:3),obj.xv(4:6),T,obj.area_drag,obj.mass,obj.Cd,Omega);
                torque = [0;0;0];
                projectedArea = obj.area_drag;
                return;
            end
            xv_ecef = ECI2ECEF(obj.currentEpoche, obj.xv');
            wind_i = wind_def(obj.currentEpoche, xv_ecef);
            v_rel_i = obj.xv(4:6) - wind_i;
            v_abs = norm(v_rel_i);
            % construct ray from satellite velocity
            wind_body = clRay();
            wind_body.origin = v_rel_i;
            wind_body.direction = -wind_body.origin/v_abs;
            wind_body.origin = QForm(obj.qw(1:4), wind_body.origin);
            wind_body.direction = QForm(obj.qw(1:4), wind_body.direction);
            wind_body.xdir = cross(wind_body.direction, obj.xv(1:3));
            wind_body.ydir = cross(wind_body.direction, wind_body.xdir);
            wind_body.xdir = wind_body.xdir/norm(wind_body.xdir);
            wind_body.ydir = wind_body.ydir/norm(wind_body.ydir);
            % discretizise beam
            resolution = p.station.beamResolution;
            if obj.hitmethod == eHitMethod.Beam
                [rays, ~] = discretizeBeam(wind_body,obj.Lc/2,p.station.beamResolution);
                % check for each beam how often he hits and choose first;
                drag_cell = cell(1,length(rays));
                parfor i = 1:length(rays)
                    loc_obj = obj;
                    seghits = [];
                    for j = 1:length(obj.segments)
                        tmphits = loc_obj.segments(j).getDragIntersections(obj.segments(j).offset, rays(i), ...
                            dens, resolution, obj.qw(5:7));
                        seghits = cat(1,seghits,tmphits);
                    end
                    if ~isempty(seghits)
                        drag_cell{i} = closestHit(seghits);
    %                     nHits = nHits+1;
    %                     beamhits(nHits) = closestHit(seghits);
                    end
                end
                empt = find(cellfun('isempty', drag_cell));
                drag_cell(empt) = [];
            else
                drag_cell = cell(1, obj.nFacets);
                n = 0;
                
                w_atmos_i = [0;0;7.29211514670698*10^-5];
                w_atmos_b = QForm(obj.qw(1:4), w_atmos_i);
                wind_rel_i =  wind_i - obj.xv(4:6);
                wind_rel_b = QForm(obj.qw(1:4), wind_rel_i);
                w_sat_i = QForm(obj.qw(1:4), obj.qw(5:7));
                w_rel_i = w_sat_i - w_atmos_i;
                w_rel_b = obj.qw(5:7) - w_atmos_b; % needed for torque
                
                drag = clDrag();
                for i = 1:length(obj.segments)
                    for j = 1:length(obj.segments(i).facets)
                        dotp = dot(obj.segments(i).facets(j).normal, wind_body.direction);
                        if ( dotp<0 )
                            n=n+1;
                            drag_cell{n} = clDrag();
                            drag_cell{n}.affectedArea = 0;
                            drag_cell{n}.projectedArea = 0;
                            drag_cell{n}.hitpos = [0;0;0];
                            drag_cell{n}.angleOfIncidence = 0;
                            drag_cell{n}.distFromLaser = 0;
                            drag_cell{n}.las2hitpos = [0;0;0];
                            if obj.hitmethod == eHitMethod.Net
                                [ pts, areapp ] = obj.segments(i).facets(j).getPointNet(resolution);
                            else
                                pts = obj.segments(i).facets(j).areabarycenter;
                                areapp = obj.segments(i).facets(j).area;
                            end
                            npts = size(pts,2);
                            for k= 1:npts
                                drag.hitpos = pts(:,k) + obj.segments(i).offset;
%                                 drag.angleOfIncidence = acos(-dotp/...
%                                     (norm(obj.segments(i).facets(j).normal)*norm(wind_body.direction)));
                                drag.angleOfIncidence = vectors2angle(-wind_body.direction, obj.segments(i).facets(j).normal);    
                                drag.affectedArea = areapp;
                                drag.projectedArea = areapp*cos(drag.angleOfIncidence);

                                eulerAng = SpinCalc('QtoEA313',obj.qw(1:4)',1e-6,0)';%quarternion to euler angles 313
                                hitpos_eci = rigid_transf_b2i(drag.hitpos,eulerAng)+obj.xv(1:3);
                                
                                [U,~] = transMatEci2Ecef(obj.currentEpoche);
                                x_ecef = U*hitpos_eci;
                                [lon, lat, height] = Geodetic(x_ecef(1:3));
                                alt = height/1000;
                                xlat = lat*p.const.Deg;
                                xlong = lon*p.const.Deg;
                                
                                wind_i = wind_def(obj.currentEpoche, x_ecef);
                                v_rel_i = obj.xv(4:6) - wind_i;
                                v_abs = norm(v_rel_i);
                                % construct ray from satellite velocity
                                wind_body = clRay();
                                wind_body.origin = v_rel_i;
                                wind_body.direction = -wind_body.origin/v_abs;
                                wind_body.origin = QForm(obj.qw(1:4), wind_body.origin);
                                wind_body.direction = QForm(obj.qw(1:4), wind_body.direction);
                                wind_body.xdir = cross(wind_body.direction, obj.xv(1:3));
                                wind_body.ydir = cross(wind_body.direction, wind_body.xdir);
                                wind_body.xdir = wind_body.xdir/norm(wind_body.xdir);
                                wind_body.ydir = wind_body.ydir/norm(wind_body.ydir);
                
                                [den,~] = nrlmsise00_mex(int32(idoy),ut,alt,xlat,xlong,f107a,f107,ap,flag);
                                if Oplus_flag
                                    dens = 1.66e-27*(16.0*den(9)) + den(6);
                                end 
                                drag = drag.init(obj, obj.segments(i), wind_body, ...
                                    obj.segments(i).facets(j), dens, obj.qw(5:7), w_rel_i, ...
                                    wind_rel_i, w_rel_b);

%                                 drag_cell{n} = drag_cell{n} + drag;
                                drag_cell{n}.force = drag_cell{n}.force + drag.force;
                                drag_cell{n}.torque = drag_cell{n}.torque + drag.torque;
                                drag_cell{n}.affectedArea = drag_cell{n}.affectedArea + drag.affectedArea;
                                if norm(drag_cell{n}.force) && norm(drag_cell{n}.torque) 
                                    drag_cell{n}.projectedArea = drag_cell{n}.projectedArea + drag.projectedArea;
                                end
                                drag_cell{n}.hitpos = drag_cell{n}.hitpos + drag.hitpos;
                                drag_cell{n}.angleOfIncidence = drag_cell{n}.angleOfIncidence + drag.angleOfIncidence;
                                drag_cell{n}.distFromLaser = drag_cell{n}.distFromLaser + drag.distFromLaser;
                                drag_cell{n}.las2hitpos = drag_cell{n}.las2hitpos + drag.las2hitpos;
                            end
                            drag_cell{n}.hitpos = drag_cell{n}.hitpos/npts;
                            drag_cell{n}.angleOfIncidence = drag_cell{n}.angleOfIncidence/npts;
                            drag_cell{n}.distFromLaser = drag_cell{n}.distFromLaser/npts;
                            drag_cell{n}.las2hitpos = drag_cell{n}.las2hitpos/npts;
                        end
                    end
                end
                empt = find(cellfun('isempty', drag_cell));
                drag_cell(empt) = [];
            end
%             hits = beamhits(1:nHits);
            % calculate force and torque from found hits
            force = [0; 0; 0];
            torque = [0; 0; 0];
            projectedArea = 0;%area to mass ratio
            for i = 1:length(drag_cell)
                force = force + drag_cell{i}.force;
                torque = torque + drag_cell{i}.torque;
                if norm(drag_cell{i}.force) && norm(drag_cell{i}.torque)
                    projectedArea = projectedArea + drag_cell{i}.projectedArea;%added by Yang
                end
            end
            % turn into eci coord sys
            force = QTForm(obj.qw(1:4), force);
            
%             obj.hitmethod = tmpmthd;
        end        
        
        %% solar radiation pressure
        function [force, torque, projectedArea] = getSolarAccel(obj, sun)
            p = clPropagator.instance();
%             global AuxParam const
            if p.AuxParam.dof == 3
                if obj.area_solar == 0
                    obj.area_solar = obj.Lc^2*pi/4;
                end
                force = obj.mass*AccelSolrad(obj.xv(1:3),sun.xv(1:3),obj.area_solar,obj.mass, ...
                    obj.Cr,p.const.P_Sol,p.const.AU);
                torque = [0;0;0];
                projectedArea = obj.area_solar;
                return;
            end
            tmpmthd = obj.hitmethod;
            obj.hitmethod = eHitMethod.Area;
            [~, ~, force, torque, projectedArea] = obj.getLaserAccel(sun);
            obj.hitmethod = tmpmthd;
        end
        
        %% magnetic torque
        function torque = getMagneticTorque(obj,Y)
            LLA = zeros(4,1);
            prpgtr = clPropagator.instance();
            mJD = prpgtr.rso.currentEpoche;
            Meci2ecef = transMatEci2Ecef(mJD);
            [lon,~,~] = Geodetic(Meci2ecef*Y(1:3));
            LLA(3) = lon;
            M = obj.mtm;
            w = [0;0;7.29211514670698*10^-5];
            
            G = prpgtr.Gnm;
            H = prpgtr.Hnm;
            
            [B_field_i,B_field_dot_i] = magnet_field(mJD, Meci2ecef*Y(1:3), Meci2ecef*Y(4:6), LLA, Meci2ecef', G, H);
            B_field_b = QForm(obj.qw(1:4),B_field_i);
            B_field_dot_b = QForm(obj.qw(1:4),B_field_dot_i);
            torque = magnet_torque(B_field_b, B_field_dot_b, w, M);
        end
        
        function [alfor,altor,irfor,irtor,alproA,irproA] = getAlbedoForceTorque(obj,Y,r_Sun)
                        
            alfor = zeros(3,1);
            altor = zeros(3,1);
            irfor = zeros(3,1);
            irtor = zeros(3,1);
            alproA = 0;
            irproA = 0;
            
            tmp_hitmethod = obj.hitmethod;
			obj.hitmethod = eHitMethod.Area;
			if obj.hitmethod == eHitMethod.Area
				i_seg = 1:length(obj.segments);
				hits_cell = arrayfun(@(x) obj.segments(x).hitsArea_albedo(Y, r_Sun),i_seg ...
                    ,'UniformOutput', false);
			end
			obj.hitmethod = tmp_hitmethod;
            for i = 1:length(hits_cell)
                nfacets = length(obj.segments(i).facets);
                for j = 1:nfacets
                    alfor = alfor + hits_cell{i}{j,1};
                    altor = altor + hits_cell{i}{j,2};
                    irfor = irfor + hits_cell{i}{j,3};
                    irtor = irtor + hits_cell{i}{j,4};
                    alproA = alproA + hits_cell{i}{j,5};
                    irproA = irproA + hits_cell{i}{j,6};
                end
            end
        end
        
        %% generic propagator (no gravity)
        %% move
        function move(obj, impulse, deltat)
            obj.xv = [ (obj.xv(1:3) + deltat*obj.xv(4:6)); (obj.xv(4:6) + QTForm(obj.qw(1:4),impulse)/obj.mass) ];
        end
        
        %% spin
        function spin(obj, moment, deltat)
            torque = moment/deltat;
            obj.qw = RK4( 'FRB', obj.qw, deltat, 0, obj.moi, ...
                    obj.invmoi, QTForm(obj.qw(1:4),torque) );
            absq = sqrt(obj.qw(1)^2 + obj.qw(2)^2 + obj.qw(3)^2 + obj.qw(4)^2);
            obj.qw(1:4) = obj.qw(1:4)./absq;
        end
        
        %% getlongestside
        function length = getLongestSide(obj)
            length = 0;
            for i = 1:3
                tmplength = obj.referenceCube(i,1)-obj.referenceCube(i,2);
                if tmplength > length
                    length = tmplength;
                end
            end
        end
        
        %% getshortestside
        function length = getShortestSide(obj)
            length = 100;
            for i = 1:3
                tmplength = obj.referenceCube(i,1)-obj.referenceCube(i,2);
                if tmplength < length
                    length = tmplength;
                end
            end
        end
        
        %% setResolution
        function setResolution(obj, res)
            for i = 1:length(obj.segments)
                obj.segments(i).resolution = res;
            end
        end
        
        function setSurfaceAttributes(obj,sas)
            if length(sas)~=length(obj.segments)
                fprintf('Uneven number of surface attributes and segments.\nCan not assign attributes.\n');
                return;
            end
            
            for i = 1:length(obj.segments)
                fprintf('Setting %s Surface Attributes: %s\n',obj.segments(i).name, sas(i).name);
                obj.segments(i).surfAtr = sas(i);
            end
        end
    end
end

