clear
clc

%% Create Propagator
% global propagator
propagator = clPropagator.instance();

%% some initial values
Y0 = [ +7144843.808; +217687.110; -506463.296; +562.650611; -1616.516697; +7358.157263; ...
    1; 0; 0; 0; 0; 0; 0];
initstate = clPropagatorInitialValues();
initstate.proptime = 10*60;
initstate.step = 30;

initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 40;
initstate.m = 40;

initstate.a_lase = 1; % laser engagements
initstate.g_lase = 1;
initstate.a_drag = 1; % drag
initstate.g_drag = 0;
initstate.a_srad = 0; % solar radiation pressure
initstate.g_srad = 0;
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets

initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects

propagator.init(initstate);

%% create Target
target = createCubeTarget(0.2, false);
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 1;
station.pointError = 1;
station.referencePower = 10000;
station.bGroundBased = false;
station.lla = [0 95 650];
    
%% Add Objects

propagator.addTarget(target);
propagator.addLaserStation(station);

%% Function Handles
fcnhandles = {@ode23, @ode45, @RK4, @RK45};

%% Start Diary
propagator.initNewSim();

%% Propagation Matlab's ODE
tic
profile on
eph = propagator.propagateTarget(fcnhandles{3},Y0,initstate.step);
toc

%% Write Output
propagator.finishSim();
save2file = sprintf('%s.mat', propagator.outfilename);
save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'timestamp', 'profiler_info');
