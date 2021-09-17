% clear propagator initstate odefcn station target Y0 sa
clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();

%% Start Diary
% propagator.initNewSim();


%% some initial values    
timestamp = now;
w0 = [0; 0; 2.9].*pi./180;
q0 = eul2quat([0, 50, 0]'.*pi./180);
Y0 = [getInitialState(0,30,600); q0; w0];
% Y0 = [Y0_dof3; q0; w0];

initstate = clPropagatorInitialValues();
% initstate.proptime = 2*24*60*60;
initstate.proptime = 60*60;
initstate.step = 1;

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

initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects

propagator.init(initstate);

%% create Target
% target = loadDSPOSEEnvisat();
target = createCubeSatTarget(1);
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0;
station.referencePower = 20000;
station.bGroundBased = true;
station.lla = [0, 0, 0];
station.beamResolution = 200;
station.bOnlyInDark = true;
station.bPulsed = false;
station.intstep = initstate.step;

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
% sa(2) = sa(1);
target.setSurfaceAttributes(sa);
% target.hitmethod = eHitMethod.Beam;
target.hitmethod = eHitMethod.Area;
    
%% Add Objects
% 
% propagator.addTarget(target);
% propagator.addLaserStation(station);

% target2DSPOSE(target, Y0, propagator.AuxParam, '');

for lt = 1:2
    if lt==1
        station.bGroundBased = true;
        station.lla = [5+9/60+40.3/360, -52-38/60-49/360, 0]; %public swimming pool kourou
    else
        station.bGroundBased = false;
        station.lla = [0 0 400];
    end
    %% Start Diary
    propagator.initNewSim();
    propagator.addTarget(target);
    propagator.addLaserStation(station);
    %% Propagation Matlab's ODE
    tic
    profile on
    % odefcn = @RK45; 
    % odefcn = @ode23;
    eph = propagator.propagateTarget_radau(Y0,initstate.step);
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

%     propagator = reconstructPerturbations(propagator, eph);

    if exist('eph', 'var')
        save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'dspose_timestamp');
    else
        save(save2file, 'propagator', 'initstate', 'Y0', 'dspose_timestamp');
    end

end
