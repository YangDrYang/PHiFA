classdef clPropagator < handle
    %clPropagator Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % OBJECTS
        station clLaserStation = clLaserStation()
        stations_param = struct('n_stations', 0, ...
            'bPulsed', false, ...
            'repRate', -1);
        
        rso clTarget = clCompoundTarget()
        sun clSun = clSun();
        
        % SIMULATION
        startmjd double          % starting MJD
        proptime double          % duration of simulation [s]
        t double                 % seconds from simulation start [s]
        Y0 (13,1) double         % initial state vector
            
        i_step                   % step counter of simulation
        
        absTol double = 1e-6;
        relTol double = 1e-11;
        
        min_height double = 100 % if below this, simulation will be aborted [km]
        
        % 'GLOBAL' PARAMETERS
        AuxParam
        eopdata
        const
        PC
        swdata
        Cnm
        Snm
        Gnm
        Hnm
        albedo
        
        % OUTPUT
        outfilename string
        output_table clSimulationOutput
        prev_logdata clSimulationOutput
        used_integrator_fcn string
        timestamps (2,1) datetime
        profiler_info struct
        
        % USER UPDATE
        last_update     % counter
    end
    
    properties (Access=private)
        initialized logical = false
    end
    
    methods (Static)
        function obj = instance(varargin)
            persistent uniqueInstance
            if nargin==1
                uniqueInstance = varargin{1};
                obj = uniqueInstance;
            elseif isempty(uniqueInstance)
                obj = clPropagator();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    
    methods (Access=private)
        function obj = clPropagator()
            fprintf('Reading propagator information...\n');
            % --------------------------------------------------------------------------
            % 
            %                   Satellite Orbit Modeling
            %
            % References:
            % Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
            % Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
            %
            % Montenbruck O., Pfleger T.; Astronomy on the Personal Computer; Springer 
            % Verlag, Heidelberg; 4th edition (2000).
            %
            % Seeber G.; Satellite Geodesy; Walter de Gruyter, Berlin, New York; 2nd
            % completely revised and extended edition (2003).
            %
            % Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
            % , New York; 3rd edition(2007).
            %
            % http://ssd.jpl.nasa.gov/?ephemerides
            %
            % http://celestrak.com/SpaceData/
            %
            % Last modified:   2018/02/11   M. Mahooti
            %
            % --------------------------------------------------------------------------

%             global const Cnm Snm swdata eopdata PC

%             run('propagator/SAT_Const');
            %% CONST VARIABLES
            % Mathematical constants
            obj.const.pi2       = 2*pi;                % 2pi
            obj.const.Rad       = pi/180;              % Radians per degree
            obj.const.Deg       = 180/pi;              % Degrees per radian
            obj.const.Arcs      = 3600*180/pi;         % Arcseconds per radian

            % General
            obj.const.MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
            obj.const.T_B1950   = -0.500002108;        % Epoch B1950
            obj.const.c_light   = 299792458.000000;    % Speed of light  [m/s]; DE430
            obj.const.AU        = 149597870700.000000; % Astronomical unit [m]; DE430

            % Physical parameters of the Earth, Sun and Moon

            % Equatorial radius and flattening
            obj.const.R_Earth   = 6378.137e3;          % Earth's radius [m]; WGS-84
            obj.const.f_Earth   = 1/298.257223563;     % Flattening; WGS-84   
            obj.const.R_Sun     = 696000e3;            % Sun's radius [m]; DE430
            obj.const.R_Moon    = 1738e3;              % Moon's radius [m]; DE430

            % Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
            obj.const.omega_Earth = 15.04106717866910/3600*obj.const.Rad; % [rad/s]; WGS-84

            % Gravitational coefficients
            obj.const.GM_Earth   = 398600.4418e9;               	 % [m^3/s^2]; WGS-84
            obj.const.GM_Sun     = 132712440041.939400e9;      	  	 % [m^3/s^2]; DE430
            obj.const.GM_Moon    = obj.const.GM_Earth/81.30056907419062; % [m^3/s^2]; DE430
            obj.const.GM_Mercury = 22031.780000e9;             	  	 % [m^3/s^2]; DE430
            obj.const.GM_Venus   = 324858.592000e9;            	  	 % [m^3/s^2]; DE430
            obj.const.GM_Mars    = 42828.375214e9;             	  	 % [m^3/s^2]; DE430
            obj.const.GM_Jupiter = 126712764.800000e9;         	  	 % [m^3/s^2]; DE430
            obj.const.GM_Saturn  = 37940585.200000e9;          	  	 % [m^3/s^2]; DE430
            obj.const.GM_Uranus  = 5794548.600000e9;           	  	 % [m^3/s^2]; DE430
            obj.const.GM_Neptune = 6836527.100580e9;           	  	 % [m^3/s^2]; DE430
            obj.const.GM_Pluto   = 977.0000009e9;        	  	 	 % [m^3/s^2]; DE430

            % Solar radiation pressure at 1 AU
            obj.const.I_Sol = 1367; % [W/m^2] (1367 W/m^2); IERS 96
            obj.const.P_Sol = 1367/obj.const.c_light; % [N/m^2] (1367 W/m^2); IERS 96

            obj.loadData();
        end
        
        function loadData(obj)
            %% FURTHER VARIABLES
            load('propagator/DE430Coeff.mat', 'DE430Coeff');
            obj.PC = DE430Coeff;

            % read Earth gravity field coefficients
%             obj.Cnm = zeros(181,181);
%             obj.Snm = zeros(181,181);
%             fid = fopen('propagator/GGM03S.txt','r');
%             for n=0:180
%                 for m=0:n
%                     temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
%                     obj.Cnm(n+1,m+1) = temp(3);
%                     obj.Snm(n+1,m+1) = temp(4);
%                 end
%             end
%             fclose(fid);
            [obj.Cnm, obj.Snm] = loadMahootiGravCoef();
            obj.Cnm(1,1) = 1;
            % read Earth orientation parameters
            fid = fopen('propagator/eop19990101.txt','r');
%             fid = fopen('propagator/eop20130101.txt','r');
            %  ----------------------------------------------------------------------------------------------------
            % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
            % |(0h UTC)           "         "          s          s          "        "          "         "     s 
            %  ----------------------------------------------------------------------------------------------------
            obj.eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
            fclose(fid);

            % read space weather data
            fid = fopen('propagator/sw19990101.txt','r');
%             fid = fopen('propagator/sw20130101.txt','r');
            %  ---------------------------------------------------------------------------------------------------------------------------------
            % |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
            % | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
            %  ---------------------------------------------------------------------------------------------------------------------------------
            obj.swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
            fclose(fid);
        end
    end
    
    methods
        
        function bool = isInitialized(obj)
            bool = obj.initialized;
        end    
        %% init
        function init(obj, initialstate)
            fprintf('Intializing propagator.\n');

            if obj.initialized == true
                obj.loadData();
            end
            
            % model parameters
            obj.AuxParam = initialstate;
            
            if initialstate.sec >= 0
                year = initialstate.year;
                mon = initialstate.mon;
                day = initialstate.day;
                hour = initialstate.hour;
                min = initialstate.min;
                sec = initialstate.sec;
                Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
            else
                Mjd_UTC = initialstate.mjd;
                jd = Mjd_UTC+2400000.5;
                [year,mon,day,hr,min,sec] = invjday(jd);
            end

            obj.AuxParam.Mjd_UTC = Mjd_UTC;
            
            obj.startmjd   = Mjd_UTC;
            obj.proptime   = initialstate.proptime;   % propagation time [s]

            % shorten PC, eopdata, swdata, Cnm, and Snm
            num = fix(obj.proptime/86400)+5;
            JD = Mjd_UTC+2400000.5;
            i = find(obj.PC(:,1)<=JD & JD<=obj.PC(:,2),1,'first');
            obj.PC = obj.PC(i:i+num,:);
            mjd = (floor(Mjd_UTC));
            i = find(mjd==obj.eopdata(4,:),1,'first');
            obj.eopdata = obj.eopdata(:,i-1:i+num);
            i = find((year==obj.swdata(1,:)) & (mon==obj.swdata(2,:)) & (day==obj.swdata(3,:)),1,'first');
            obj.swdata = obj.swdata(:,i-4:i+num);
            obj.Cnm = obj.Cnm(1:obj.AuxParam.n+1,1:obj.AuxParam.n+1);
            obj.Snm = obj.Snm(1:obj.AuxParam.n+1,1:obj.AuxParam.n+1);

            obj.initialized = true;
            
            % read magnetic IGRF-12 Gauss Coefficients
            [obj.Gnm, obj.Hnm] = loadMagCoef();
            
            % read magnetic IGRF-12 Gauss Coefficients
            obj.albedo = loadAlbedoInfraredCoef(obj.AuxParam.albedo_model);
        end
        %% addtarget
        function addTarget(obj, targetobj)
            if ~obj.isInitialized()
                fprintf('Initialice before adding targets.\n');
                return;
            end
            obj.rso = targetobj;
        end
        %% addlaserstation
        function addLaserStation(obj, laserobj, varargin)
            if nargin < 2
                disp('You should be calling the function with at least a laser station to add');
                return;
            elseif nargin == 2
                laserobj.AuxParam = obj.AuxParam;
            elseif nargin > 2
                laserobj.AuxParam = varargin{1};
            end
            
            if ~obj.isInitialized()
                fprintf('Initialice before adding laser stations.\n');
                return;
            end
            obj.stations_param.n_stations = obj.stations_param.n_stations + 1;
            fprintf('Adding Laserstation %s.\n', laserobj.name);
            obj.station(obj.stations_param.n_stations) = laserobj;
%             propsecs = -laserobj.intstep:laserobj.intstep:(obj.proptime+laserobj.intstep);
%             if propsecs(end) < obj.proptime
%                 propsecs(end+1) = propsecs(end)+laserobj.intstep;
%             end
            propsecs = -30:laserobj.intstep:obj.proptime+30;
%             propsecs = 0:laserobj.intstep:obj.proptime;
            if propsecs(end) < obj.proptime
                propsecs(end+1) = propsecs(end)+30;
            end
            itsteps = length(propsecs);
%             laserobj.ephemerides(1:7,1:itsteps) = zeros(7,itsteps);
            if laserobj.bPulsed
                obj.stations_param.bPulsed = true;
                % quick and dirty - the repetition rates of all
                % laserstations is set to the rate of the station that has
                % beed added first
                if obj.stations_param.repRate < 0 
                    obj.stations_param.repRate = laserobj.repetitionRate;
                else
                    laserobj.repetitionRate = obj.stations_param.repRate;
                end
            end

            if laserobj.bGroundBased
%                 wgs84 = wgs84Ellipsoid('meters');
                ecef = [0,0,0,0,0,0];
                rr = laserobj.lla(3)*1000 + obj.const.R_Earth;
                ecef(1)= rr*cos(laserobj.lla(1)*pi/180)*cos(laserobj.lla(2)*pi/180);
                ecef(2)= rr*cos(laserobj.lla(1)*pi/180)*sin(laserobj.lla(2)*pi/180);
                ecef(3)= rr*sin(laserobj.lla(1)*pi/180);
%                 [ecef(1),ecef(2),ecef(3)] = geodetic2ecef(wgs84,laserobj.lla(1),...
%                     laserobj.lla(2),laserobj.lla(3));
                sprintf('Propagating ground-based station.');
                for i = 1:itsteps
                    laserobj.ephemerides(1,i) = propsecs(i);
                    mjd = obj.startmjd+propsecs(i)/86400;
                    laserobj.ephemerides(2:7,i) = ECEF2ECI(mjd,...
                        ecef).'; % \citep{mahooti2018matlab}
                end
            else
%                 Y0 = getInitialState(laserobj.lla(1), laserobj.lla(2), laserobj.lla(3));
%                 lasoptions = odeset('RelTol',1e-11,'AbsTol',1e-6);
%                 [t,yout] = ode23(@laserobj.Accel,propsecs,Y0,lasoptions);
%                 laserobj.ephemerides = cat(1,t',yout');
% %                 laserobj.ephemerides(2:7,:) = yout';
                clPropagator.instance();
                Y0 = laserobj.coor;
%                 lasoptions = odeset('RelTol',1e-11,'AbsTol',1e-6);
%                 [t,yout] = ode23(@laserobj.Accel,propsecs,Y0,lasoptions);
                lasoptions = rdpset('RelTol',1e-13,'AbsTol',1e-16);
                [t,yout] = radau(@laserobj.Accel,propsecs,Y0,lasoptions);
                laserobj.ephemerides = cat(1,t',yout');
                clear global
            end
        end
        
        function resetLaserStations(obj)
            obj.stations_param = struct('n_stations', 0, ...
                'bPulsed', false, ...
                'repRate', -1);
            obj.station = clLaserStation();
        end
         %% propagateTarget
        function ephemerides = propagateTarget(obj, odefcn, Y0, step)
            obj.Y0 = Y0;
            tic;
            obj.used_integrator_fcn = func2str(odefcn);
            % global AuxParam
            fprintf('#->\t->\t##################\t<-\t<-#\n');
            fprintf('Starting simulated Deorbiting of %s.\n', obj.rso.name);
            fprintf('#->\t->\t##################\t<-\t<-#\n');
            if (obj.AuxParam.a_lase || obj.AuxParam.g_lase) && obj.stations_param.bPulsed
                ustep = 1/obj.stations_param.repRate;
                if ustep > step
                    ustep = step;
                end
                % maybe consider choosing ustep to be bigger than f_rep and
                % averaging over several pulses
            else
                ustep = step;
            end
            steps = 0:ustep:obj.proptime;
            if steps(end)<obj.proptime
%                 steps = [steps steps+ustep];
                steps = 0:ustep:obj.proptime+ustep;
            end
            obj.last_update = -5;
            max_output = ceil(steps(end)/obj.AuxParam.n_log_only_every_xth_step+1);
            obj.output_table(1:max_output) = clSimulationOutput();
            ephemerides(1:length(steps), 1:14) = 0;
            ephemerides(1,2:14) = Y0';
            ephemerides(1,1) = 0;
            Y_old = Y0;
            obj.initOutput(length(steps)-1);
%             [accelfcn, odeargin] = obj.getInputArgForODEFCN(odefcn);
            Y_hist = NaN(8,14);%to save historical state for eight epochs
            Y_hist(1,1) = 0;
            Y_hist(1,2:14) = Y0;
            flagLaser = 0;
            if obj.AuxParam.a_lase && obj.AuxParam.g_lase
                flagLaser = 1;
            end
            obj.i_step = 2;
            for i = 1:length(steps)-1
                tspan = [steps(i) steps(i+1)];
                
                
                if strcmp(func2str(odefcn),'ABM8')
                    if i < 8
                        odefcn_ = @RK45;
%                         odefcn_ = @RK4;
                        [accelfcn, odeargin] = obj.getInputArgForODEFCN(odefcn_); % RK45 for the first 7 steps
                        obj.AuxParam.a_lase = 0;
                        obj.AuxParam.g_lase = 0;
                        if isempty(odeargin)
                            [t,Y_new] = odefcn_(accelfcn,tspan,Y_old);
                        else
                            [t,Y_new] = odefcn_(accelfcn,tspan,Y_old,odeargin(:));
                        end
                    else
                        [accelfcn, ~] = obj.getInputArgForODEFCN(odefcn);
                        if flagLaser
                            obj.AuxParam.a_lase = 1;
                            obj.AuxParam.g_lase = 1;                        
                        end
                        Y_tmp = odefcn(accelfcn, tspan(2)-tspan(1), Y_hist);
                        Y_new = Y_tmp';
                        t = tspan(2);
                        Y_hist(1:7,:) = Y_hist(2:8,:);
                        Y_hist(8,1) = tspan(2);
                        Y_hist(8,2:end) = Y_tmp';
                    end
                else 
                    [accelfcn, odeargin] = obj.getInputArgForODEFCN(odefcn);
                    if isempty(odeargin)
                        [t,Y_new] = odefcn(accelfcn,tspan,Y_old);
                    else
                        [t,Y_new] = odefcn(accelfcn,tspan,Y_old,odeargin(:));
                    end
                end

                Y = Y_new(end,:)';
                
%                 obj.output_table(obj.i_step).t = t(end);
                fprintf('Time: %d of %d\n',t(end), obj.proptime);
                
                MJD_UTC = obj.startmjd + t(end)/86400;
                [~,~,UT1_UTC,~,~,~,~,~,TAI_UTC] = IERS(obj.eopdata,MJD_UTC,'l');
                [~, ~, ~, TT_UTC, ~] = timediff(UT1_UTC,TAI_UTC);
                MJD_TT = MJD_UTC+TT_UTC/86400;
                MJD_TDB = Mjday_TDB(MJD_TT);
                [~,~,~,~,~,~,~, ...
                    ~,~,~,r_Sun,~] = JPL_Eph_DE430(MJD_TDB);
                for ii = 1:obj.stations_param.n_stations
                    obj.station(ii).xv(1:3) = obj.station(ii).getPosition(obj.t(end));
                    obj.station(ii).lsType = 2;%laser type: laser
                    [incAng, bLaserSwitch] = obj.station(ii).laserSwitchCriterion(Y(1:3),Y(4:6), r_Sun);
                    obj.prev_logdata.incAng = incAng;
                    if (obj.AuxParam.a_lase || obj.AuxParam.g_lase) && obj.station(ii).bPulsed && ...
                            bLaserSwitch && ...
                            ~obj.station(ii).inducePulse(obj.rso)

                        [laserimpulse, lasermomentum, laserforce, lasertorque, projarea] = ...
                            obj.rso.getLaserAccel(obj.station(ii));

                        obj.prev_logdata.laser.for = laserforce;
                        obj.prev_logdata.laser.tor = lasertorque;
                        obj.prev_logdata.laser.imp = laserimpulse;
                        obj.prev_logdata.laser.mom = lasermomentum;
                        obj.prev_logdata.laser.projarea = projarea;
                        if ~norm(laserimpulse) && ~norm(lasermomentum)
                            obj.prev_logdata.lt = 0;
                        else
                            obj.prev_logdata.lt = 2;
                        end

                        %update velocity and angular velocity
                        if obj.AuxParam.a_lase
                            Y(4:6) = Y(4:6) + laserimpulse/obj.rso.mass;
                            obj.prev_logdata.energy(4:5) = obj.prev_logdata.energy(4:5) -...
                                dot(laserimpulse,laserimpulse) /obj.rso.mass/2.0;
                        end
                        if obj.AuxParam.g_lase
                            Y(11:13) = Y(11:13) + (obj.rso.invmoi*lasermomentum);
                            obj.prev_logdata.energy(4:5) = obj.prev_logdata.energy(4:5) -...
                                lasermomentum'*obj.rso.moi*lasermomentum/2.0;  
                        end
                    end
                end
                Y_old = Y(:);
                if strcmp(func2str(odefcn),'ABM8') && i < 8%for first seven points
                    Y_hist(i+1,1) = tspan(2);
                    Y_hist(i+1,2:14) = Y;
                end
%                 obj.output_table(obj.i_step).Y = Y_old(:);
                ephemerides(i+1,1) = t(end);
                ephemerides(i+1,2:14) = Y';
                fprintf('#->\t->\t##################\t<-\t<-#\n');
%                 obj.updateOutput(i);
                if obj.rso.getHeight()<obj.min_height*1000
                    fprintf('Object reached minimum height. Simulation is being aborted.\n');
                    obj.output_table = obj.output_table(1:obj.i_step);
                    break;
                end
                obj.prev_logdata.t = t(end);
                obj.prev_logdata.simtime = toc;
                obj.prev_logdata.steptime = obj.output_table(2).simtime ...
                    - obj.output_table(1).simtime;
                if obj.AuxParam.b_logging && ...
                        mod(i,obj.AuxParam.n_log_only_every_xth_step)==0
                    obj.output_table(obj.i_step) = obj.prev_logdata;
                    obj.output_table(obj.i_step).Y = Y_old(:);
                    obj.output_table(obj.i_step).t = t(end);
                    obj.i_step = obj.i_step + 1;
                end
            end
            ephemerides = ephemerides(1:i+1,:)';
        end
        %% proptarget_radau
        function ephemerides = propagateTarget_radau(obj, Y0, tepochs)
            if obj.AuxParam.a_lase || obj.AuxParam.g_lase
                fprintf('Pulsed laser wont be able to be considered using radau(), CW laser is not recommended\n');
            end
            obj.Y0 = Y0;
            options = rdpset('RelTol',obj.relTol,'AbsTol',obj.absTol);
%             pulsecs = 0:step:obj.proptime;
%             waitbartitle = sprintf('Propagating target: %s', obj.rso.name);
%             obj.f = waitbar(0,waitbartitle);
%             [t,yout] = radau(@obj.Accel,pulsecs,Y0,options);
            [t,yout] = radau(@obj.Accel,tepochs,Y0,options);
%             close(obj.f);
            ephemerides(:,1) = t;
            ephemerides(:,2:14) = yout;
            ephemerides = ephemerides';
            
%             for i = 1:length(ephemerides)
%                 obj.i_step=i;
%                 obj.Accel(ephemerides(i,1),ephemerides(i,2:14))
%             end
        end
        %% accel
        function dY = Accel(obj,t,Y)
            obj.t = t;
            dY = obj.rso.Accel(t,Y);
        end
        %% accel SCT
        % for some reason the SCT decided to invert order of input
        % arguments compared to \gls{matlab}s ODE functions
        function dY = Accel_sct(obj,Y,t)
            dY = obj.Accel(t,Y);
        end
        %% updateProp
        function updatePropagation(obj, t)
            new_update = floor(t/obj.proptime*1000);
            if mod(new_update,20)==0 && new_update ~= obj.last_update
                fprintf('%i%%.', new_update/10);
                obj.last_update = new_update;
            end
        end
        %% initoutput
        function initOutput(obj, steps)
            obj.output_table(steps) = clSimulationOutput();
        end
        %% initnewsim
        function initNewSim(obj)
            nSim = 1;
            tmpname = sprintf('logfiles%ssim%04d.mat',filesep,nSim);
            while exist(tmpname, 'file')
                nSim=nSim+1;
                tmpname = sprintf('logfiles%ssim%04d.mat',filesep,nSim);
            end
            tmpname = sprintf('logfiles%ssim%04d.out',filesep,nSim);
            if exist(tmpname, 'file')
                delete(tmpname);
            end
            obj.outfilename = sprintf('logfiles%ssim%04d',filesep,nSim);
            
            if obj.AuxParam.b_debugging
                profile on;
                diary2file = sprintf('%s.out', obj.outfilename);
                diary(diary2file);
            end
            obj.timestamps(1) = datetime('now');
        end
        %% called after propagation
        function finishSim(obj)
            s = profile('status');
            if obj.AuxParam.b_debugging && strcmp(s.ProfilerStatus, 'on')
                obj.profiler_info = profile('info');
                profile clear
            end
            obj.timestamps(2) = datetime('now');
            diary off
            
            % make entry
            fid = fopen('logfiles/simulation_diary.txt', 'a');
            inital_state_str = sprintf('%09.3f\t%09.3f\t%09.3f\t%06.3f\t%06.3f\t%06.3f', ...
                obj.Y0(1)/1000,obj.Y0(2)/1000,obj.Y0(3)/1000,obj.Y0(4)/1000,obj.Y0(5)/1000,obj.Y0(6)/1000);
            initial_attitude_str = sprintf('%.3f\t%.3f\t%.3f\t%.3f\t%07.3f\t%07.3f\t%07.3f', ...
                obj.Y0(7),obj.Y0(8),obj.Y0(9),obj.Y0(10),obj.Y0(11),obj.Y0(12),obj.Y0(13));
            fprintf(fid,'%s\t%s\t%s\t%d\t%.3f\t%s\t%s\t%d%d\t%d\t%d\t%d%d\t%d%d\t%s\t%.3f\t%.3f\n', ...
                datestr(obj.timestamps(1)), string(obj.timestamps(2)-obj.timestamps(1)), ...
                obj.outfilename, obj.AuxParam.dof, ...
                obj.proptime/3600, inital_state_str, initial_attitude_str,...
                obj.AuxParam.a_lase, obj.AuxParam.g_lase, obj.station(1).bPulsed, ...
                obj.stations_param.n_stations, obj.AuxParam.a_drag, obj.AuxParam.g_drag, obj.AuxParam.a_srad, ...
                obj.AuxParam.g_srad, obj.rso.name, obj.rso.mass, obj.rso.Lc );
            fclose(fid);
        end
        %% update output
        function updateOutput(obj, step)
            global glout
            obj.output_table(step) = glout;
        end
        
        function [accelfcn, varargout] = getInputArgForODEFCN(obj, odefcn)
            c = func2str(odefcn);
            if strcmp(c,'ABM8')
                varargout = cell(1);
                accelfcn = @obj.Accel;
            elseif strcmp(c,'RK45')
                varargout = cell(1);
                varargout{1} = obj.absTol;
                accelfcn = @obj.Accel_sct;
            elseif strcmp(c, 'RK4')
                varargout = cell(1);
                accelfcn = @obj.Accel_sct;
            elseif strcmp(c(1:3), 'ode')
                varargout = cell(1);
                varargout{1} = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
                accelfcn = @obj.Accel;
            else
                varargout = cell(1);
                accelfcn = @obj.Accel;
            end
        end
    end
end

