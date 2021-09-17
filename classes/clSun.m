classdef clSun < clLaserStation
    %clSun Since the radiation pressure of the sun can be similarly treated
    %as a laser station we just inherit from it and modify
    
    properties
        shaMod = 2 % shadow model
        const_int double %[W/m^2]
        satellite_shadowed logical = false % true if earth shadows satellite
    end
    
    properties (Constant)
        solar_radiation_pressure double = 1361; %[W/m^2]
        r_earth double = 6371008.8; %mean earth radius [m]
    end
    
    methods
        function obj = clSun()
            obj.beamResolution = 100;
        end
        
        function varargout = inducePulse(obj, targetobj, r_Sun)          
%             global eopdata const
            p = clPropagator.instance();
            
%             MJD_UTC = obj.currentEpoche;
%             
%             [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
%             [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
%             MJD_TT = MJD_UTC+TT_UTC/86400;
%             MJD_UT1 = MJD_UTC+UT1_UTC/86400;
%             
%             MJD_TDB = Mjday_TDB(MJD_TT);
%             [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
%              r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_TDB);
         
            obj.xv(1:3) = r_Sun;
            
            % calculate new pointing dir of laser station
            obj.pDirECI = obj.xv(1:3) - targetobj.xv(1:3);
            targetDist = norm(obj.pDirECI);
%             fprintf('Target distance to sun: %f m\n', targetDist);
            
%             %check for shadow
%             %assuming perfectly cylindrical shadow cone \citep{frueh2013}
%             refdist = sin(acos(dot(obj.pDirECI, targetobj.xv(1:3))/...
%                 (targetDist*norm(targetobj.xv(1:3)))))...
%                 *norm(targetobj.xv(1:3)); %derived from equation 8
            
            obj.raySAT.origin = QForm(targetobj.qw(1:4), obj.pDirECI);
            obj.pDirECI = obj.pDirECI/norm(obj.pDirECI);
            obj.raySAT.direction = -obj.raySAT.origin;
            obj.raySAT.direction = obj.raySAT.direction/norm(obj.raySAT.direction);

            % calculate x- and ydir arrays
            refvec = [0 0 1];

            obj.raySAT.xdir = cross(refvec, targetobj.xv(1:3)-obj.xv(1:3));
            obj.raySAT.xdir = QForm(targetobj.qw(1:4), obj.raySAT.xdir);
            obj.raySAT.xdir = obj.raySAT.xdir/norm(obj.raySAT.xdir);

            obj.raySAT.ydir = cross(obj.raySAT.xdir, targetobj.xv(1:3)-obj.xv(1:3));
            obj.raySAT.ydir = QForm(targetobj.qw(1:4), obj.raySAT.ydir);
            obj.raySAT.ydir = obj.raySAT.ydir/norm(obj.raySAT.ydir);
            
            obj.const_int = getSolarIntensity(targetDist/p.const.AU, obj.currentEpoche)*...
                getNormalizedSourceIntensity(obj.shaMod, targetobj.xv(1:3), r_Sun);

            if nargout > 0
                varargout = cell(1,1);
                varargout{1} = struct('poe', 0, 'h0', -1, 'za', -1, 'tdist', targetDist...
                    ,'fried',-1, 'rytov',-1, 'maxint',obj.const_int ...
                    ,'turberror',-1, 'nscr',-1);
            end
        end
        
        function [pulseFluence, intensity] = getBeamPropertiesAtTarget(obj, ~)
            pulseFluence = 0;
            intensity = obj.const_int;
        end
        
        function Uout = atmosAbsorption(~,Uin,~)
            Uout = Uin;
        end
        
        function bool = laserSwitchCriterion(obj, r_eci, ~)
            bool = getNormalizedSourceIntensity(obj.shaMod, r_eci, obj.xv(1:3));
        end
    end
end

