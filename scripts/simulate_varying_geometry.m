% simulate_varying_geometry.m

clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 2.9].*pi./180;
q0 = eul2quat([0, 50, 0].*pi./180)';
% Y0 = [getInitialState(0,30,600); q0; w0];
Y0 = [6644264.68005270;-1955163.54312126;-842230.248328992;2292.78323879521;6199.46750663702;3667.02707298642; ...
    q0; w0];

initstate = clPropagatorInitialValues();
initstate.proptime = 1000;
initstate.mjd = 52389.5801736113;
initstate.sec = -1;
initstate.step = 1;
initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 40;
initstate.m = 40;

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



%% create Target
% target = createRandEllipsoidTarget(height, xyratio, yzratio, percentNoise, bSolid);
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0;
station.referencePower = 10000;
station.beamResolution = 200;
station.bOnlyInDark = true;
station.bGroundBased = true;
station.lla = [5+9/60+40.3/360, -52-38/60-49/360, 0]; %public swimming pool kourou


%% some more target conf (needs knowledge of pulselength, wavelength)
target.hitmethod = eHitMethod.Beam;

propagator.init(initstate);
propagator.station.bPulsed = true;    
propagator.station.trackError = 1.5;
propagator.station.referencePower = 20000;
propagator.station.PulseLength = 1E-9;
propagator.station.Wavelength = 1064E-9;

sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = false;
% sa(2) = sa(1);
target.setSurfaceAttributes(sa);

propagator.addLaserStation(station);
%% Add Objects

height = 0.1;
ratios = [];

noise = [0 0.1 1 2.5 5];

for i = 1:length(noise)
    target = createRandEllipsoidTarget(height, xyratio, yzratio, noise(i), 1);
    
    %% Start Diary
    propagator.initNewSim();
    propagator.addTarget(target);

    %% Propagation
    tic
    profile on
    odefcn = @RK45; 
    % odefcn = @ode23;
    eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
    toc

    %% Write Output
    propagator.finishSim();
    save2file = sprintf('%s.mat', propagator.outfilename);

    dim_describ = ["repfreq ", num2str(noise(i)); ...
        "groundbased ", num2str(station.bGroundBased); ...
        "pulsed laser ", num2str(station.bPulsed)];

    if exist('eph', 'var')
        save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'dim_describ');
    else
        save(save2file, 'propagator', 'initstate', 'Y0', 'dim_describ');
    end

end

