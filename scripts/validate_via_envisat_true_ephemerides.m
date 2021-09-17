% clear propagator initstate odefcn station target Y0 sa
clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();

%% Start Diary
% propagator.initNewSim();


%% some initial values    
timestamp = now;
[Y0_dof3, mjd_utc_start, date] = loadEnvisatMahootiInitialState();
% Y0 = [Y0_dof3; q0; w0];
initialstate.year = date.year;
initialstate.mon = date.mon;
initialstate.day = date.day;
initialstate.hour = date.hour;
initialstate.min = date.min;
initialstate.sec = date.sec;

initstate = clPropagatorInitialValues();
initstate.proptime = 1588*60;
initstate.step = 10;
initstate.mjd = mjd_utc_start;
initstate.sec = -1;

initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 0; % using DSPOSE functions       
initstate.n = 80;
initstate.m = 80;

initstate.a_lase = 1; % laser engagements
initstate.g_lase = 1;
initstate.a_drag = 1; % drag
initstate.g_drag = 1;
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 1;
initstate.a_sun = 1; % gravity of sun
initstate.a_moon = 1; % gravity of moon
initstate.a_planets = 1; % gravity of other planets

initstate.a_solidEarthTides = 1; % solid earth tides
initstate.a_oceanTides = 1; % ocean tides
initstate.a_relativity = 1; % relativity effects

propagator.init(initstate);

%% create Target
target = loadDSPOSEEnvisat();
%% create Laserstation
station = createOHigginsLaserstation();

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
target.setSurfaceAttributes(sa);
target.hitmethod = eHitMethod.Beam;
    
%% Add Objects

propagator.addTarget(target);
% propagator.addLaserStation(station);

% target2DSPOSE(target, Y0, propagator.AuxParam, '');

% for i = 1:10
    %% Start Diary
    propagator.initNewSim();
    q0 = angle2quat(rand*360*pi/180, rand*360*pi/180, rand*360*pi/180,'ZYX')';
    w0 = [rand*6-3; rand*6-3; rand*6-3].*pi./180;
    Y0 = [Y0_dof3; q0; w0];
    %% Propagation Matlab's ODE
    tic
    profile on
    % odefcn = @RK45; 
    % odefcn = @ode23;
    propagator.used_integrator_fcn = "RADAU II";
    tepochs = 0:initstate.step:initstate.proptime;
    eph = propagator.propagateTarget_radau(Y0,tepochs);
    toc
    %% compare
%     envisat_eph = loadEnvisatMahootiTrueEphemerides();
%     envisat_eph_eci = ecef2eci_ephemerides(envisat_eph, mjd_utc_start);
%     compareEphemerides(envisat_eph_eci, eph);
%     title(num2str(i));
    %% Write Output
    propagator.finishSim();
    save2file = sprintf('%s.mat', propagator.outfilename);

    % cd dspose
    % !./dspose_exec
    % cd ..
    % dspose_timestamp = getLatestDSPOSETimestamp();
    dspose_timestamp = '';

    % propagator = reconstructPerturbations(propagator, eph);

    if exist('eph', 'var')
        save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'dspose_timestamp');
    else
        save(save2file, 'propagator', 'initstate', 'Y0', 'dspose_timestamp');
    end

% end
