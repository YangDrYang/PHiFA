
%% Create Propagator
% global propagator
propagator = clPropagator.instance();

%% Start Diary
propagator.initNewSim();
diary2file = sprintf('%s.out', propagator.outfilename);
diary(diary2file);


%% some initial values    
timestamp = now;
% Y0 = [ +7144843.808; +217687.110; -506463.296; +562.650611; -1616.516697; +7358.157263; ...
%     1; 0; 0; 0; 0; 0; 0];
% Y0 = [getInitialState(0,-60,400); eul2quat([randn()*pi/2 randn()*pi/2 randn()*pi/2])'; 0; 0; 0];
Y0 = [684436.915615222;3371704.43120338;-5840024.59108916;-7629.37020630040;387.128285881048;-670.593822301558;0.856939656428299;-0.452354332018004;-0.0327725378753924;-0.244859029456344;-1.70872106236373e-18;-3.34769893055079e-18;2.92487586236297e-17];
initstate = clPropagatorInitialValues();
initstate.proptime = 12*60*60;
initstate.step = 10;

initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 40;
initstate.m = 40;

initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 1; % drag
initstate.g_drag = 1;
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 1;
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
station.trackError = 0;
station.pointError = 0.5;
station.referencePower = 10000;
% station.bGroundBased = false;
% station.lla = [0 120 650];

%% some more target conf (needs knowledge of pulselength, wavelength)
sa = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa.bUseDeviation = true;
target.setSurfaceAttributes(sa);
target.hitmethod = eHitMethod.Beam;
    
%% Add Objects

propagator.addTarget(target);
propagator.addLaserStation(station);


%% Propagation Matlab's ODE
tic
profile on
odefcn = @RK45; 
% odefcn = @ode23;
eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
toc

%% Write Output
propagator.finishSim();
save2file = sprintf('%s.mat', propagator.outfilename);
if exist('eph', 'var')
    save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'timestamp');
else
    save(save2file, 'propagator', 'initstate', 'Y0', 'timestamp');
end


