clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();

stella_file = 'sampledata/L7A06001L56.PRE';
[stella_eph, mjd0] = loadStellaEph(stella_file);
stella_eph_eci = ecef2eci_ephemerides(stella_eph, mjd0);

%% some initial values    
timestamp = now;
w0 = [0; 0; 0].*pi./180;
% q0 = angle2quat([0, 0, 0].*pi./180);
q0 = eul2quat([0, 0, 0]'.*pi./180);
% [Y0_dof3, mjd0, date] = loadEnvisatMahootiInitialState();
Y0_dof3 = stella_eph_eci(2:7,1);
Y0 = [Y0_dof3; q0; w0];
initstate = clPropagatorInitialValues();
initstate.proptime = stella_eph(1,end);
initstate.step = 30;
initstate.mjd = mjd0;
initstate.sec = -1;

initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 100;
initstate.m = 100;

initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 0; % drag
initstate.g_drag = 0;
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 0;
initstate.a_sun = 1; % gravity of sun
initstate.a_moon = 1; % gravity of moon
initstate.a_planets = 1; % gravity of other planets

initstate.a_solidEarthTides = 1; % solid earth tides
initstate.a_oceanTides = 1; % ocean tides
initstate.a_relativity = 1; % relativity effects



%% create Target
% target = createCubeTarget(0.2, false);
% target = loadDSPOSEEnvisat();
target = loadStellaSatellite();
% plate = makePlate(0.2, 0.05);
% target = createTarget(plate, 'Solid Plate');
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0.5;
station.referencePower = 10000;
station.bGroundBased = false;
station.lla = [0 120 650];

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
% sa(2) = sa(1);
target.setSurfaceAttributes(sa);
target.hitmethod = eHitMethod.Area;


% target2DSPOSE(target, Y0, propagator.AuxParam, '');
    
propagator.init(initstate);
propagator.addTarget(target);
% propagator.addLaserStation(station);

% target2DSPOSE(target, Y0, propagator.AuxParam, '');
for i = 1:3
    %% Start Diary
    propagator.initNewSim();
%     q0 = angle2quat([rand*360, rand*360, rand*360].*pi./180);
    q0 = eul2quat([rand*360, rand*360, rand*360]'.*pi./180);
    w0 = [rand*6-3; rand*6-3; rand*6-3].*pi./180;
    Y0 = [Y0_dof3; q0; w0];

    %% Propagation
    tic
    profile on
    odefcn = @RK45; 
    % odefcn = @ode23;
    eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
%     eph = propagator.propagateTarget_radau(Y0,initstate.step);
    toc

    %% compare
    % close all;
    compareEphemerides(stella_eph_eci, eph);

    %% Write Output
    propagator.finishSim();
    save2file = sprintf('%s.mat', propagator.outfilename);
%     saveAllFigures(propagator.outfilename);

    % cd dspose
    % !./dspose_exec
    % cd ..
    dspose_timestamp = '';

    % propagator = reconstructPerturbations(propagator, eph);

    if exist('eph', 'var')
        save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'dspose_timestamp', 'stella_eph_eci');
    else
        save(save2file, 'propagator', 'initstate', 'Y0', 'dspose_timestamp', 'stella_eph_eci');
    end
end


