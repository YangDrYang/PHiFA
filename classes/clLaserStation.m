classdef clLaserStation < cl3dofObject
    %clLaserStation
    
    properties
        
        pDirECI(3,1) double = [1; 0; 0]       % ECI [m]
        pDirAzEl(2,1) double = [0; 30]       % deg local Horizon - not used
        raySAT clRay = clRay()
        
        % Position
        bPulsed logical = true              % true if pulsed laser, false if \gls{cw}
        ephemerides double                  % position information for laser station
        
        bNoLaser logical = false            %false if laser station is used
        bGroundBased logical = true         % true if station is groundbased 
                                            % false if space-borne
                                            
        intstep double = -1                 % if space-borne, to set for orbit integration step 
        
        lla(1,3) double = [-35.3 149.0 0.8]; %canberra
        % if groundbased: lat [deg], long [deg], alt [km] 
        % if spacebased: raan [deg], inclination [deg], height [km]
        coor(1,6) double = [0,0,0,0,0,0];
        % if spacebased: six orbital state in Cartesian coordinate system
        
        % Laser engagement angle
        engAng double = 0
        
        % Laser Station Parameters
        repetitionRate double = 1           % [Hz]
        Aperture double = 1                 % aperture of telescope [m]
        Wavelength double = 1064*10^-9      % frequency doubled Nd:YAG Laser [m]
        PulseLength double = 1*10^-9        % [s]
        referencePower double = 10000       % [W]
        % laserSwitchCriteria
        minEle double = 0                   % Minimum Elevation [deg]
        minD2Earth double = 200000;         % Minimum Distance of Laser Beam to Earth when space-based[m]
        
        % Pointing Parameters
        trackError double = 2;              % error due to uncertain orbit parameteres [mrad]
        pointError double = 0.2424;         % error due to mechanical and electronical 
                                            % restrictions of laser station [mrad]
%         minHitProbability double = 0.01
        minHitProbability double = 0.5
        bOnlyInDark logical = false
        
        % Efficiencies
        telescopeEff double = 0.95          % average power direct after telescope
        atmosTrans double = 0.82            % default transmission by atmosphere
                                            % replaced by lowtran
        bUseLowtran logical = false         % set false if lowtran is not available
        atmosTransLowtran double = 1
                                            
        % Simulation parameter
        focusBias double = 0;               % focus too far or too short [%]
        focusFluctuation double = 0;        % deviation [%]
        powerFluctuation double = 0;        % deviation [%]
        
        % Tubulence
        tumod = @tumodSLCDay                % SLC Model
        xn                                  % x values of intensity distribution of Uturb
        yn                                  % y values of intensity distribution of Uturb
        Uturb                               % intensity distribution after turbulence
        % not needed
        Uvac                                % intensity distribution after vacuum
        x1                                  % x values of intensity distribution of Uin
        y1                                  % y values of intensity distribution of Uin
        Uin                                 % intensity distribution right after telescope
        
        % Beam parameters
        beamType eLaserBeam = eLaserBeam.Gaussian
        beamResolution double = 100         % [Beams per meter]
        % Gaussian
        FWHM double = 1                     % Full Width of Half Maximum [m]
        % Tophat
        tophatIntensity double = 100000     % J / m^2 meaning Fluence at Object per Pulse
        tophatDiameter double = -1          % [m]
        bUsetophatIntensity logical = false % false if Fluence shall be estimated via referencePower
        
        %probably uneeded
%         gaussianMaximumFluence double = -1 % [J / m^2]
%         gaussianVariance double = -1 % [m^2]
        BeamQualitySquared double = 1.48    % one is ideal [-]
        
        lsType int8 = 2 %laser type: laser 
    end
    
    methods
        function obj = clLaserStation()
            %clLaserStation 
        end
        
        function init(obj)
            transmittedPower = obj.referencePower * obj.telescopeEff;
            if obj.beamType == eLaserBeam.TopHat && ~obj.bUsetophatIntensity
                if obj.tophatDiameter < 0
                    obj.tophatDiameter = obj.Aperture;
                end    
                obj.tophatIntensity = transmittedPower/(pi * (obj.tophatDiameter/2)^2);
            end
            
            if obj.bUseLowtran
                % atmosphere absorbtion and scattering
                [~,~,isloaded] = pyversion;
                if ~isloaded
                    pyversion('C:\ProgramData\Miniconda3\envs\py36\python.exe');
                end
            end
        end
        
        %% INDUCE PULSE
        function [aborted, varargout] = inducePulse(obj, targetobj)
            % calculate new pointing dir of laser station
            aborted=false;
            obj.pDirECI = obj.xv(1:3) - targetobj.xv(1:3);
            targetDist = norm(obj.pDirECI);
%             fprintf('Target distance: %f m\n', targetDist);
            obj.raySAT.origin = QForm(targetobj.qw(1:4), obj.pDirECI);
            obj.pDirECI = obj.pDirECI/norm(obj.pDirECI);
            obj.raySAT.direction = -obj.raySAT.origin;
            % pointing error
            randvec = -ones(size(obj.raySAT.direction)) + 2*rand(size(obj.raySAT.direction));
            randvec = cross(obj.raySAT.direction, randvec);
            randvec = randvec/norm(randvec);
            pointerror = sqrt(obj.trackError^2 + obj.pointError^2) / 1.0e6;
            errorvec = randvec*pointerror*randn()*norm(obj.raySAT.direction);
            obj.raySAT.direction = obj.raySAT.direction+errorvec;
            obj.raySAT.direction = obj.raySAT.direction/norm(obj.raySAT.direction);
            actual_pointinerror = 1.0e6*norm(errorvec)/targetDist;
            fprintf('Pointing error: %f murad\n',actual_pointinerror);
            
            minBeamWidth = obj.Wavelength*targetDist/obj.Aperture;
            targetDiameter = targetobj.getLongestSide();
            hitProbability = 1 - exp(-(minBeamWidth+targetDiameter)^2 ...
                /(8*targetDist^2*pointerror^2));
            fprintf('Hit Probability: %f\n',hitProbability);
            
            if hitProbability>obj.minHitProbability
                % calculate x- and ydir arrays
                if obj.xv(3) >= 0
                    refvec = [0 0 1];
                else
                    refvec = [0 0 -1];
                end

                if norm(cross(refvec,obj.raySAT.direction)) < 1e-10
                    refvec = [1 0 0];
                end

                obj.raySAT.xdir = cross(refvec, targetobj.xv(1:3)-obj.xv(1:3));
                obj.raySAT.xdir = QForm(targetobj.qw(1:4), obj.raySAT.xdir);
                obj.raySAT.xdir = obj.raySAT.xdir/norm(obj.raySAT.xdir);

                obj.raySAT.ydir = cross(obj.raySAT.xdir, targetobj.xv(1:3)-obj.xv(1:3));
                obj.raySAT.ydir = QForm(targetobj.qw(1:4), obj.raySAT.ydir);
                obj.raySAT.ydir = obj.raySAT.ydir/norm(obj.raySAT.ydir);


                % get arriving beam distribution
                D2 = targetDiameter*2; % diameter of area of interest at target plane
                za = acos(abs(dot(obj.pDirECI,obj.xv(1:3))/(norm(obj.pDirECI)*norm(obj.xv(1:3))))); % zenith angle
                fprintf('Zenith angle: %f deg\n', za*180/pi);
                h0 = obj.getHeight(); % height
                fprintf('Laser station height: %f m\n', h0);

                if obj.beamType == eLaserBeam.Gaussian     
                    fprintf('Initializing Gaussian beam\n');
                    rng('shuffle')

                    [aborted, obj.xn,obj.yn,obj.Uturb,propstruct] = ...
                        propBeam(obj,targetDist,D2,obj.tumod,0.05,0.01,za,h0);
                    if aborted
                        fprintf('Bad luck. Try another time.\n');
                        aborted=true;
                    elseif obj.bGroundBased
                        tm = obj.getTransmittance_lowtran(h0, za);
                        obj.atmosTransLowtran = tm;
                        fprintf('Atmosphere Transmission: %.3f\n', tm);
                    end
                else
                    fprintf('Initializing TopHat beam\n');
                    deltan = 0.01;
                    D2 = targetobj.getLongestSide()*2; % diameter of area of interest at target plane
                    N = D2/deltan;
                    N = round((N-2)/2)*2+2;
                    if N>1000
                        N=1000;
                        deltan=D2/N;
                    end
                    % making mesh
                    [obj.xn, obj.yn] = meshgrid((-N/2 : N/2-1) * deltan);
                    obj.x1 = obj.xn; obj.y1 = obj.yn;
                    [~, r1] = cart2pol(obj.x1, obj.y1);
                    if ~obj.bUsetophatIntensity
                        % constructing tophat beam with energy in beam and radius
                        obj.Uin = tophat(r1,obj.tophatDiameter)*obj.referencePower...
                            *obj.telescopeEff/(pi*(obj.tophatDiameter/2)^2);
                        tm = obj.getTransmittance_lowtran(h0, za);
                        obj.Uturb = obj.Uin*tm;
                        fprintf('Atmosphere Transmission: %.3f\n', tm);
                        obj.Uvac = obj.Uin;
                    else
                        obj.Uin = tophat(r1,obj.tophatDiameter)*obj.tophatIntensity;
                        obj.Uvac = obj.Uin;
                        obj.Uturb = obj.Uin;
                    end

                    if nargout > 0
                        propstruct.fried = -1;
                        propstruct.rytov = -1;
                        [propstruct.maxint,~] = findMaxIntensity(obj.Uturb,obj.xn,obj.yn);
                        propstruct.turberror = 0;
                        propstruct.nscr = -1;
                    end
                end
            else
                fprintf('Probablity of hitting too small. Skipping Pulse.\n');
                aborted=true;
                if nargout > 0
                    propstruct.fried = -1;
                    propstruct.rytov = -1;
                    propstruct.maxint = 0;
                    propstruct.turberror = 0;
                    propstruct.nscr = -1;
                end
            end
            
            if nargout > 1
                varargout = cell(1,1);
                varargout{1} = struct('poe', actual_pointinerror, 'h0', h0, 'za', za, 'tdist', targetDist...
                    ,'fried',propstruct.fried, 'rytov',propstruct.rytov, 'maxint',propstruct.maxint...
                    ,'turberror',propstruct.turberror, 'nscr',propstruct.nscr);
            end
        end
        
        function h0 = getHeight(obj)
            if obj.bGroundBased
                h0 = obj.lla(3)*1000;
            else
                h0 = norm(obj.xv(1:3))-6378.137e3;
            end
        end
        
        function [pulseFluence, intensity] = getBeamPropertiesAtTarget(obj, pt)
            % satellite coordinate system
            % find x- and y-coordinates in ray coordinate system
            A = [obj.raySAT.direction(:), obj.raySAT.xdir(:), obj.raySAT.ydir(:)];
            rcond_a = rcond(A);
            if rcond_a < 1e-20 || isnan(rcond_a)
                sv = pinv(A) * ( pt(:) - obj.raySAT.origin(:) );
            else
                sv = A \ ( pt(:) - obj.raySAT.origin(:) );
            end
            xcor = sv(2);
            ycor = sv(3);
            % interpolate beam intensity
%             intensity = interp2(obj.xn,obj.yn,obj.Uturb,xcor,ycor,'linear');
            intensity = 8.8542187817E-12*299792458/2* ...
                (qinterp2(obj.xn,obj.yn,obj.Uturb,xcor,ycor,2))^2 * obj.atmosTransLowtran;
%             intensity = abs( intensity );
            if isnan(intensity)
                intensity = 0;
            end
            % convert intensity to fluence
            pulseFluence = intensity/obj.repetitionRate;
        end
        
        function t = getTransmittance_wnuk2015changes(obj, ~, za)
            el = pi/2-za;
            x = 1/(cos(el)+0.025*exp(-11*cos(el))); %wnuk2015changes eq 7
            t = obj.atmosTrans^x;
        end
        
        function t = getTransmittance_lowtran(obj, h0, za)
            if obj.bUseLowtran
                p.model=5;
                p.h1=h0/1000;
                p.angle=za*180/pi;
                p.wlshort= 200;
                p.wllong= 30000;
                T = lowtran_transmission(p);
                trans = squeeze(xarray2mat(T{'transmission'}));
                wl_nm = xarray2mat(T{'wavelength_nm'});

                t = interp1q(wl_nm, trans, obj.Wavelength*10^9);
            else
                t = interpTransmittance(h0, za, obj.Wavelength);
            end
        end
        
        function bool = laserPulsed(obj, sec)
            if ~clPropagator.initialiced()
                fprintf('Initialice propagator and pass laser station.\n');
                return;
            elseif ~obj.bPulsed
                bool = false;
            else
                pulsetimes = obj.ephemerides(1,1):1/obj.repetitionRate:obj.ephemerides(1,end);
                c = min(abs(pulsetimes-sec));
                if c <= obj.PulseLength/2
                    bool = true;
                else
                    bool = false;
                end
            end
        end
        
        function [incAng,bool] = laserSwitchCriterion(obj, r_eci, v_eci, r_sun)
            dotp0 = dot(obj.xv(1:3),r_sun);
            if obj.bOnlyInDark && ... station in shadow and satellite lit
                    ( ( obj.bGroundBased && dotp0<0 ) ... if station gb, station must be in shadow
                    || getNormalizedSourceIntensity(1, r_eci, r_sun)==0 ) % satellite must be lit
                bool=false;
                incAng = 0;
                return;
            end
            r_las2obj = r_eci-obj.xv(1:3);
            incAng = vectors2angle(r_las2obj,v_eci);
            dotp1 = dot(r_las2obj,v_eci);
            if dotp1 < 0 % deccelerating it
                dotp2 = dot(r_las2obj,obj.xv(1:3));
                if obj.bGroundBased
                    bool = dotp2 > 0;
                else
                    dist2earth = point_to_line([0; 0; 0], obj.xv(1:3), r_eci);
                    bool =  dist2earth > 6378.137e3+obj.minD2Earth ...
                        || dotp2 >= 0;
                end
            else
                bool = false;
            end
        end

        function r_eci = getPosition(obj, sec)
            r_eci = [0;0;0];
            if ~isempty(obj.ephemerides) && ...
                    ( ( sec > obj.ephemerides(1,1) ) || ( sec < obj.ephemerides(1,end) ) )
                r_eci(1) = interp1q(obj.ephemerides(1,:)',obj.ephemerides(2,:)',sec);
                r_eci(2) = interp1q(obj.ephemerides(1,:)',obj.ephemerides(3,:)',sec);
                r_eci(3) = interp1q(obj.ephemerides(1,:)',obj.ephemerides(4,:)',sec);
            else
                fprintf('MJD not in propagation time window or ephemerides not initialiced.\n');
            end
        end
        %% out of date        
        function mjdpulses = getPulseTimes(obj)
            mjdpulses = obj.ephemerides(1,1):1/obj.repetitionRate/86400:obj.ephemerides(1,end);
        end

        function fLocalFluence = getBeamFluence(obj, pt)
            distance = point_to_line(pt, obj.raySAT.origin, obj.raySAT.origin+obj.raySAT.direction);
            if ~obj.bUseRefFluence
                %calc new reference fluence at distance to target from station...
            end
            fLocalFluence = obj.referenceFluence;
            if obj.beamType == eLaserBeam.TopHat
                if distance > obj.referenceDiameter/2
                    fLocalFluence = 0.0;
                end
            else
                if obj.gaussianMaximumFluence < 0
                    beamEnergy = obj.referenceFluence * pi * (obj.referenceDiameter/2)^2;
                    obj.gaussianVariance = obj.FWHM * ( 2 * sqrt( 2 * log( 2 ) ) );
                    obj.gaussianMaximumFluence = beamEnergy / ( 2 * pi * obj.gaussianVariance^2 );
                end
                fLocalFluence = obj.gaussianMaximumFluence * exp( - ( (distance^2) / ( 2 * obj.gaussianVariance^2 ) ) );
            end
        end
        
        function getIntensity(obj, pt)
            
        end
        
        function [x1, y1, Uleave, delta1] = getLeavingFluenceDistribution(obj,tDist)
%             delta1 = 1/obj.beamResolution;
            Dref = obj.Aperture/2;
            N = 512;
            delta1 = Dref/N;
            k = 2*pi / obj.Wavelength; % optical wavenumber [rad/m]
            [x1, y1] = meshgrid((-N/2 : N/2-1) * delta1);
            [~, r1] = cart2pol(x1, y1);
            
            D1 = obj.Wavelength*tDist / Dref;    % width of central lobe [m]

            % Leaving
            if obj.beamType == eLaserBeam.Gaussian
                Uleave = 2*obj.referencePower / (pi * D1) .* exp(-1i*k/(2*tDist) * r1.^2) ...
                    .* exp(-(r1/(4*D1)).^2);
            else
                Uleave = obj.referenceFluence .* exp(-1i*k/(2*tDist) * r1.^2)...
                    .* tophat(r1, obj.Aperture/2);
            end
        end
        
        
%         function fBeamArea = getBeamArea(obj, pt)
%             %calcBeamRadiusAtGivenDistance This gives back the area the
%             %laser beam would had have if diffraction would been take into
%             %account.
%             %   Although modern lasers can via adaptive optics and artifial
%             %   guide stars compensate for that.
% %             fDistance = line_on_point(pt, obj.stateVectorTrans(1:3), obj.stateVectorTrans(1:3)+obj.pointingDirection);
% %             fBeamRadius = (2*obj.fWavelength*obj.fBeamQualitySquared*fDistance)/(Pi*obj.fAperture);
%             fBeamRadius = obj.referenceDiameter;
%             fBeamArea = pi * fBeamRadius^2;
%         end        
    end
end

