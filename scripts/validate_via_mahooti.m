clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 0].*pi./180;
q0 = angle2quat([0, 0, 0].*pi./180);
[Y0_dof3, mjd_utc_start, date] = loadEnvisatMahootiInitialState();
Y0 = [Y0_dof3; q0; w0];
initstate = clPropagatorInitialValues();
initstate.proptime = 1588*60;
initstate.step = 10;
initstate.mjd = mjd_utc_start;
initialstate.sec = -1;

initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 80;
initstate.m = 80;

initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 0; % drag
initstate.g_drag = 0;
initstate.a_srad = 0; % solar radiation pressure
initstate.g_srad = 0;
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets

initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects



%% create Target
% target = createCubeTarget(0.2, false);
target = loadDSPOSEEnvisat();
% plate = makePlate(0.2, 0.05);
% target = createTarget(plate, 'Solid Plate');
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0.5;
station.referencePower = 10000;
% station.bGroundBased = false;
% station.lla = [0 120 650];

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
target.setSurfaceAttributes(sa);
target.hitmethod = eHitMethod.Beam;

% propagator.init(initstate);
% propagator.addTarget(target);
% propagator.addLaserStation(station);
% 
% target2DSPOSE(target, Y0, propagator.AuxParam, '');

%% loop
for i = 1:1
    
    if i==2
        initstate.a_grav = 0; % harmonic terms of gravity
        initstate.g_grav = 0;
        initstate.a_drag = 1; % drag
        initstate.g_drag = 1;
        initstate.a_srad = 0; % solar radiation pressure
        initstate.g_srad = 0;
    elseif i==3
        initstate.a_grav = 0; % harmonic terms of gravity
        initstate.g_grav = 0;
        initstate.a_drag = 0; % drag
        initstate.g_drag = 0;
        initstate.a_srad = 1; % solar radiation pressure
        initstate.g_srad = 1;
    end
    
    propagator.init(initstate);
    propagator.addTarget(target);
    propagator.addLaserStation(station);

%     target2DSPOSE(target, Y0, propagator.AuxParam, '');    
    
    %% Start Diary
    propagator.initNewSim();

    %% Propagation
    tic
    profile on
    odefcn = @RK45; 
    % odefcn = @ode23;
%     eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
    eph = propagator.propagateTarget_radau(Y0,initstate.step);
    toc
    
    %% compare
    if i==1
        load logfiles/mahooti_envisat_grav.mat Eph
    elseif i==2
        load logfiles/mahooti_envisat_drag.mat Eph
    else
        load logfiles/mahooti_envisat_srp.mat Eph
    end
    Eph = Eph';
    close all;
    compareEphemerides(Eph, eph);
    
    %% Write Output
    propagator.finishSim();
    save2file = sprintf('%s.mat', propagator.outfilename);
    saveAllFigures(propagator.outfilename);
    
    if exist('eph', 'var')
        save(save2file, 'propagator', 'initstate', 'Eph', 'eph', 'Y0');
    else
        save(save2file, 'propagator', 'initstate', 'Eph', 'Y0');
    end

end


