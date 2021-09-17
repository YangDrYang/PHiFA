classdef cl6dofObject < cl3dofObject
    %cl6dofObject This class provides an interface to the propagator class, so
    %each class that inherites from this can be easily propagated
    
    properties
        %propagation variables       
        moi(3,3) double = zeros(3) % inertial tensor
        invmoi(3,3) double = zeros(3) % inverse inertial tensor
        mtm(3,3) double = zeros(3) % magnetic tensor matrix
    end
    
    methods
        function obj = cl6dofObject()
            
        end
        
        %% 6 DOF
        function dY = Accel(obj, t, Y)
            fprintf('### Accel %s: time %f\n', obj.name, t);
            p = clPropagator.instance();
            
            logdata = clSimulationOutput();
%             logdata.t = t;
            logdata.Y = Y;
            
            laserstations = p.station;
            sun = p.sun;

            MJD_UTC = p.AuxParam.Mjd_UTC+t/86400;
            force = [0;0;0];
            torque = [0;0;0];
            Wf = 0;%work done by the non-conservative forces
            Wt = 0;%work done by the non-conservative torques
            
            for i = 1:length(laserstations)
                laserstations(i).currentEpoche = MJD_UTC;
            end
            sun.currentEpoche = MJD_UTC;
            obj.currentEpoche = MJD_UTC;
            obj.xv = Y(1:6);
            Y(7:10) = norm_quat(Y(7:10)); 
            obj.qw = Y(7:13);
            
            [x_pole,y_pole,UT1_UTC,LOD,~,~,~,~,TAI_UTC] = IERS(p.eopdata,MJD_UTC,'l');
            [~, ~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
            MJD_TT = MJD_UTC+TT_UTC/86400;
            MJD_UT1 = MJD_UTC+UT1_UTC/86400;
            P = PrecMatrix(p.const.MJD_J2000,MJD_TT);
            N = NutMatrix(MJD_TT);
            T = N * P;
            E = PoleMatrix(x_pole,y_pole) * GHAMatrix(MJD_UT1) * T;

            % Difference between ephemeris time and universal time
            % JD = MJD_UTC+2400000.5;
            % [year, month, day, hour, minute, sec] = invjday(JD);
            % days = finddays(year, month, day, hour, minute, sec);
            % ET_UT = ETminUT(year+days/365.25);
            % MJD_ET = MJD_UTC+ET_UT/86400;
            % [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
            %  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_ET);

%             if p.AuxParam.a_sun || p.AuxParam.a_moon || (p.AuxParam.a_srad || p.AuxParam.g_srad) ...
%                     || p.AuxParam.a_planets
                MJD_TDB = Mjday_TDB(MJD_TT);
                [r_Mercury,r_Venus,~,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
                    r_Neptune,r_Pluto,r_Moon,r_Sun,~] = JPL_Eph_DE430(MJD_TDB);
%             end

            % Acceleration due to harmonic gravity field
%             [force, torque] = obj.getGravityAccel_PointMassAnalytic([0;0;0], const.GM_Earth);
            if p.AuxParam.a_grav || p.AuxParam.g_grav
                [potential, force, torque] = obj.getGravityAccel(r_Sun, r_Moon, E, UT1_UTC, TT_UTC, x_pole, y_pole);
                logdata.grav.for = force;
                logdata.grav.tor = torque;
                Ue = potential;
                Wgg = (i)*dot(torque,obj.qw(5:7));
            end
            
            % laser engagements
            if (p.AuxParam.a_lase || p.AuxParam.g_lase)
                for i = 1:length(laserstations)
                    laserstations(i).xv(1:3) = laserstations(i).getPosition(t);
                    laserstations(i).lsType = 2;%laser type
                    [incAng, bLaserSwitch] = laserstations(i).laserSwitchCriterion(Y(1:3), Y(4:6), r_Sun);
                    logdata.incAng = incAng;
                    if ~laserstations(i).bPulsed && ...
                             bLaserSwitch && ...
                            ~laserstations(i).inducePulse(obj)
                        
                        [laserimpulse, lasermomentum, laserforce, lasertorque, projarea] = ...
                            obj.getLaserAccel(laserstations(i));
                        
%                         if laserstations(i).bPulsed && laserstation(i).laserPulsed(MJD_UTC)
%                             v = (obj.xv(4:6) + QTForm(obj.qw(1:4),laserimpulse)/obj.mass);
%                             w = obj.qw(5:7) + invI*lasermomentum;
%                         elseif ~laserstations(i).bPulsed
                            force = force + laserforce;
                            torque = torque + lasertorque;
                            Wf = Wf + t*dot(laserforce,obj.xv(4:6));
                            Wt = Wt + t*dot(lasertorque,obj.qw(5:7));
%                         end
                        logdata.laser.for = laserforce;
                        logdata.laser.tor = lasertorque;
                        logdata.laser.imp = laserimpulse;
                        logdata.laser.mom = lasermomentum;
                        logdata.laser.projarea = projarea;
                        logdata.lt = 1;
                    else
                        logdata.laser.for = [0;0;0];
                        logdata.laser.tor = [0;0;0];
                        logdata.laser.imp = [0;0;0];
                        logdata.laser.mom = [0;0;0];
                        logdata.laser.projarea = 0; 
                    end
                end
            end
            
            % Luni-solar perturbations
            if (p.AuxParam.a_sun)
                a = AccelPointMass(Y(1:3),r_Sun,p.const.GM_Sun);
                force = force + a*obj.mass;
                logdata.sun.for = a*obj.mass;
            end

            if (p.AuxParam.a_moon)
                a = AccelPointMass(Y(1:3),r_Moon,p.const.GM_Moon);
                force = force + a*obj.mass;
                logdata.moon.for = a*obj.mass;
            end

            % Planetary perturbations
            if (p.AuxParam.a_planets)
                a = AccelPointMass(Y(1:3),r_Mercury,p.const.GM_Mercury);
                a = a + AccelPointMass(Y(1:3),r_Venus,p.const.GM_Venus);
                a = a + AccelPointMass(Y(1:3),r_Mars,p.const.GM_Mars);
                a = a + AccelPointMass(Y(1:3),r_Jupiter,p.const.GM_Jupiter);
                a = a + AccelPointMass(Y(1:3),r_Saturn,p.const.GM_Saturn);
                a = a + AccelPointMass(Y(1:3),r_Uranus,p.const.GM_Uranus);    
                a = a + AccelPointMass(Y(1:3),r_Neptune,p.const.GM_Neptune);
                a = a + AccelPointMass(Y(1:3),r_Pluto,p.const.GM_Pluto);
                
                force = force + a*obj.mass;
                logdata.plan.for = a*obj.mass;
            end

            % Solar radiation pressure
            if (p.AuxParam.a_srad || p.AuxParam.g_srad)
                sun.xv(1:3) = r_Sun;
                if sun.laserSwitchCriterion(obj.xv(1:3), obj.xv(4:6))
        %                 a = a + AccelSolrad(Y(1:3),r_Sun,AuxParam.area_solar,AuxParam.mass, ...
        %                     AuxParam.Cr,const.P_Sol,const.AU);
                    sun.inducePulse(obj, r_Sun);
                    sun.lsType = 1; %laser type: sun
                    [sunforce, suntorque, projarea] = obj.getSolarAccel(sun);
                    
%                     force = force + sunforce;
%                     torque = torque + suntorque;
%                     logdata.srad.for = sunforce;
%                     logdata.srad.tor = suntorque;
                    if p.AuxParam.a_srad
                        force = force + sunforce;
                        logdata.srad.for = sunforce;
                        logdata.srad.projarea = projarea;
                        Wf = Wf + t*dot(sunforce,obj.xv(4:6));
                    else
                        logdata.srad.for = [0;0;0];
                        logdata.srad.projarea = 0;
                    end

                    if p.AuxParam.g_srad
                        torque = torque + suntorque;
                        logdata.srad.tor = suntorque;
                        Wt = Wt + t*dot(suntorque,obj.qw(5:7));
                    else
                        logdata.srad.tor = [0;0;0];
                    end
                else
                    logdata.srad.for = [0;0;0];
                    logdata.srad.tor = [0;0;0];
                end
                
            end

            % Atmospheric drag
            if (p.AuxParam.a_drag || p.AuxParam.g_drag)
                
                Omega = p.const.omega_Earth-0.843994809*1e-9*LOD; % IERS [rad/s]
                if ~p.AuxParam.draggradient
                    % Atmospheric density

%                     Omega = p.const.omega_Earth; % IERS [rad/s] % just to compare with DSPOSE
                    dens = nrlmsise00(MJD_UTC,E*Y(1:3),UT1_UTC,TT_UTC);

%                     Oplus_flag = 1;        
% 
%                     [yy, mon, day, hour, minute, sec] = invjday(MJD_UTC+2400000.5);
%                     days = finddays(yy, mon, day, hour, minute, sec);
%                     idoy = floor(days);
%                     ut = hour*3600+minute*60+sec; % seconds in day (UT)
% 
%                     [lon, lat, height] = Geodetic(E*Y(1:3));
%                     alt = height/1000;
%                     xlat = lat*p.const.Deg;
%                     xlong = lon*p.const.Deg;
% 
%                     i = find((yy==p.swdata(1,:)) & (mon==p.swdata(2,:)) & (day==p.swdata(3,:)),1,'first');
%                     sw = p.swdata(:,i);        
%                     ap(1) = sw(23); % Arithmetic average of the 8 Ap indices for the day
%                     ap(2) = sw(15); % 3 hr AP index for current time
%                     sw_1 = p.swdata(:,i-1);
%                     % Define Solar Flux Values
%                     if floor(hour/3) == 0
%                         ap(3) = sw_1(22); % 3 hr AP index for 3 hrs before current time
%                     else
%                         ind = 15+floor(hour/3)-1;
%                         ap(3) = sw(ind); % 3 hr AP index for 3 hrs before current time
%                     end
% 
%                     if floor(hour/6) == 0
%                         if floor(hour/3) == 0
%                             ap(4) = sw_1(21); % 3 hr AP index for 6 hrs before current time
%                         else    
%                             ap(4) = sw_1(22); % 3 hr AP index for 6 hrs before current time
%                         end
%                     else
%                         ind = 15+floor(hour/3)-2;
%                         ap(4) = sw(ind); % 3 hr AP index for 6 hrs before current time
%                     end
% 
%                     if floor(hour/9) == 0
%                         if floor(hour/6) == 0   
%                             if floor(hour/3) == 0
%                                 ap(5) = sw_1(20); % 3 hr AP index for 9 hrs before current time
%                             else
%                                 ap(5) = sw_1(21); % 3 hr AP index for 9 hrs before current time
%                             end
%                         else
%                             ap(5) = sw_1(22); % 3 hr AP index for 9 hrs before current time
%                         end
%                     else
%                         ind = 15+floor(hour/3)-3;
%                         ap(5) = sw(ind); % 3 hr AP index for 9 hrs before current time
%                     end
% 
%                     sw_2 = p.swdata(:,i-2);
% 
%                     % Average of eight 3 hr AP indicies from 12 to 33
%                                              % hrs prior to current time
%                     if floor(hour/12) == 0
%                         if floor(hour/9) == 0   
%                             if floor(hour/6) == 0
%                                 if floor(hour/3) == 0
%                                     ap(6) = (sw_1(19)+sw_1(18)+sw_1(17)+sw_1(16)...
%                                         +sw_1(15)+sw_2(22)+sw_2(21)+sw_2(20))/8.0; 
%                                 else
%                                     ap(6) = (sw_1(20)+sw_1(19)+sw_1(18)+sw_1(17)+sw_1(16)...
%                                         +sw_1(15)+sw_2(22)+sw_2(21)+sw_2(20))/8.0;                         
%                                 end
%                             else
%                                 ap(6) = (sw_1(21)+sw_1(20)+sw_1(19)+sw_1(18)...
%                                     +sw_1(17)+sw_1(16)+sw_1(15)+sw_2(22))/8.0;
%                             end
%                         else
%                             ap(6) = (sw_1(22)+sw_1(21)+sw_1(20)+sw_1(19)...
%                                  +sw_1(18)+sw_1(17)+sw_1(16)+sw_1(15))/8.0;
%                         end
%                     else
%                         ind = 15+floor(hour/3)-4;
%                         ap(6) = (sw(ind)+sw_1(22)+sw_1(21)+sw_1(20)...
%                              +sw_1(19)+sw_1(18)+sw_1(17)+sw_1(16))/8.0;
%                     end
% 
%                     sw_3 = p.swdata(:,i-3);
%                     % Average of eight 3 hr AP indicies from 36 to 57
%                                              % hrs prior to current time 
%                     if floor(hour/12) == 0
%                         if floor(hour/9) == 0   
%                             if floor(hour/6) == 0
%                                 if floor(hour/3) == 0
%                                     ap(7) = (sw_2(19)+sw_2(18)+sw_2(17)+sw_2(16)...
%                                         +sw_2(15)+sw_3(22)+sw_3(21)+sw_3(20))/8.0; 
%                                 else
%                                     ap(7) = (sw_2(20)+sw_2(19)+sw_2(18)+sw_2(17)+sw_2(16)...
%                                         +sw_2(15)+sw_3(22)+sw_3(21)+sw_3(20))/8.0;                         
%                                 end
%                             else
%                                 ap(7) = (sw_2(21)+sw_2(20)+sw_2(19)+sw_2(18)...
%                                     +sw_2(17)+sw_2(16)+sw_2(15)+sw_3(22))/8.0;
%                             end
%                         else
%                             ap(7) = (sw_2(22)+sw_2(21)+sw_2(20)+sw_2(19)...
%                                  +sw_2(18)+sw_2(17)+sw_2(16)+sw_2(15))/8.0;
%                         end
%                     else
%                         ind = 15+floor(hour/3)-4;
%                         ap(7) = (sw_1(ind)+sw_2(22)+sw_2(21)+sw_2(20)...
%                              +sw_2(19)+sw_2(18)+sw_2(17)+sw_2(16))/8.0;
%                     end
% 
%                     f107 = sw(27);     % adjusted solar radio noise flux (jansky)
%                     f107a = sw(29);    % adjusted 81-day average F10 (jansky)
% 
%                     flag = ones(1,23);
%                     flag(9) = -1;
% 
%                     [den,~] = nrlmsise00_mex(int32(idoy),ut,alt,xlat,xlong,f107a,f107,ap,flag);
%                     if Oplus_flag
%                         dens = 1.66e-27*(16.0*den(9)) + den(6);
%                     end
                
                
    %                 a = a + AccelDrag(dens,Y(1:3),Y(4:6),T,AuxParam.area_drag,AuxParam.mass,AuxParam.Cd,Omega);
                    if ~isinf(dens)
                        [dragforce, dragtorque, projarea] = obj.getDragAccel(dens, T, Omega);

                        if p.AuxParam.a_drag
                            force = force + dragforce;
                            logdata.drag.for = dragforce;
                            logdata.drag.projarea = projarea;
                            Wf = Wf + t*dot(dragforce,obj.xv(4:6));
                        else
                            logdata.drag.for = [0;0;0];
                            logdata.drag.projarea = 0;
                        end

                        if p.AuxParam.g_drag
                            torque = torque + dragtorque;
                            logdata.drag.tor = dragtorque;
                            Wt = Wt + t*dot(dragtorque,obj.qw(5:7));                        
                        else
                            logdata.drag.tor = [0;0;0];
                        end
                    else
                        logdata.drag.for = [0;0;0];
                        logdata.drag.tor = [0;0;0];
                    end
                else
                    [dragforce, dragtorque, projarea] = obj.getDragGradientAccel(T, Omega);
                    if p.AuxParam.a_drag
                        force = force + dragforce;
                        logdata.drag.for = dragforce;
                        logdata.drag.projarea = projarea;
                        Wf = Wf + t*dot(dragforce,obj.xv(4:6));
                    else
                        logdata.drag.for = [0;0;0];
                        logdata.drag.projarea = 0;
                    end

                    if p.AuxParam.g_drag
                        torque = torque + dragtorque;
                        logdata.drag.tor = dragtorque;
                        Wt = Wt + t*dot(dragtorque,obj.qw(5:7));                        
                    else
                        logdata.drag.tor = [0;0;0];
                    end
                end
            end

            % Magnetic torque: Eddy-Current Torque
            if (p.AuxParam.g_mag)
                
                magtorque = obj.getMagneticTorque(Y(1:6));
                
                if p.AuxParam.g_mag
                    torque = torque + magtorque;
                    logdata.mag.tor = magtorque;
                    Wt = Wt + t*dot(magtorque,obj.qw(5:7));
                else
                    logdata.mag.tor = [0;0;0];
                end
                
            end
            
            % Earth albedo and infrared radiation force and torque
            if (p.AuxParam.a_erad || p.AuxParam.g_erad)
                [alfor,altor,irfor,irtor,alproA,irproA] = obj.getAlbedoForceTorque(Y(1:6),r_Sun);
                
                if p.AuxParam.a_erad
                    force = force + alfor;
                    logdata.albedo.for = alfor;
                    force = force + irfor;
                    logdata.infrared.for = irfor;
                    logdata.albedo.projarea = alproA;
                    logdata.infrared.projarea = irproA;
                    Wf = Wf + t*dot(alfor,obj.xv(4:6));
                    Wf = Wf + t*dot(irfor,obj.xv(4:6));
                else
                    logdata.albedo.for = [0;0;0];
                    logdata.albedo.projarea = 0;
                    logdata.infrared.for = [0;0;0];
                    logdata.infrared.projarea = 0;                    
                end
                
                if p.AuxParam.g_erad
                    torque = torque + altor;
                    logdata.albedo.tor = altor;
                    torque = torque + irtor;
                    logdata.infrared.tor = irtor;
                    Wt = Wt + t*dot(altor,obj.qw(5:7));
                    Wt = Wt + t*dot(irtor,obj.qw(5:7));
                else
                    logdata.albedo.tor = [0;0;0];
                    logdata.infrared.tor = [0;0;0];                    
                end
            end
            
            % Relativistic Effects
            if (p.AuxParam.a_relativity)
                a = Relativity(Y(1:3),Y(4:6));
                force = force + a*obj.mass;
                logdata.rel.for = force;
                Wf = Wf + t*dot(a*obj.mass,obj.xv(4:6));
            end
            
            v = obj.xv(4:6);
            w = obj.qw(5:7);
            
            logdata.grav.spt = -p.const.GM_Earth*Y(1:3)/norm(Y(1:3))^3.*obj.mass;
            a = force./obj.mass;
            qd = differQuat( obj.qw(1:4), w );
            wd = differAngVel(obj.moi,obj.invmoi,w,torque);
            
            dY = [v;a;qd;wd];
            logdata.dY = dY;
            
            % energy anlysis: work-energy balance 
            % REF.: Modeling and Simulation of Long-term
                  % Rotational Dynamics of Large Space Debris, Sec. 5.3.2
            
            Kt = 0.5*obj.mass*dot(obj.xv(4:6),obj.xv(4:6));
            Kr = 0.5*obj.qw(5:7)'*obj.moi*obj.qw(5:7);
            Ka = Kt + Kr; %kinetic energy
            deltaUe = -obj.mass*p.const.omega_Earth*(obj.xv(1)*obj.xv(5) - obj.xv(2)*obj.xv(4));
            logdata.energy(1) = Ue + Ka - Wgg + deltaUe;
            logdata.energy(2) = Ue + Kt + deltaUe;
            logdata.energy(3) = Kr - Wgg;
            logdata.energy(4) = Kr - Wgg - Wt;
            logdata.energy(5) = Ue + Ka - Wgg - (Wf + Wt) + deltaUe;
            
            if p.used_integrator_fcn == "RK45"
                if abs(mod(t,1.0))<1e-4% to record the last epoch results in each integration step, 1.0 is subject to change according to the integration step.
                   p.prev_logdata = logdata;
                end
            else
                p.prev_logdata = logdata;
            end
            if p.used_integrator_fcn == "RADAU II"
                p.prev_logdata.t = t;
                p.output_table = cat(1,p.output_table,p.prev_logdata);
            end
        end
    end
    
    methods (Abstract)
        [impulse, momentum, force, torque] = getLaserAccel(obj, laserstation)
        [force, torque] = getGravityAccel(obj, varargin)
        [force, torque] = getDragAccel(obj, dens, T, Omega)
        [force, torque] = getDragGradientAccel(obj, T, Omega)
        [force, torque] = getSolarAccel(obj, varargin)
        torque = getMagneticTorque(obj,varagin)
        [alfor,altor,irfor,irtor] = getAlbedoForceTorque(obj,varagin)
    end
end