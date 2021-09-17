clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 2.9].*pi./180;
q0 = eul2quat([0, 50, 0]'.*pi./180);
% Y0 = [getInitialState(0,30,600); q0; w0];
Y0 = [6644264.68005270;-1955163.54312126;-842230.248328992;2292.78323879521;6199.46750663702;3667.02707298642; ...
    q0; w0];

initstate = clPropagatorInitialValues();
initstate.proptime = 1000;
initstate.mjd = 52389.5801736113;
initstate.sec = -1;
initstate.step = 0.1;
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
% target = createCubeTarget(0.25, true);
% target = createConeTarget(0.3, 0.2, true);
% target = createCubeTarget(0.2,0.1,true);
% target = loadDSPOSEEnvisat();
% plate = makePlate(0.4, 0.05);
% target = createTarget(plate, 'Solid Plate');
target = createCubeSatTarget(3);
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0;
station.referencePower = 10000;
station.beamResolution = 200;
station.bOnlyInDark = true;
station.bGroundBased = true;
station.lla = [5+9/60+40.3/360, -52-38/60-49/360, 0]; %public swimming pool kourou
station.intstep = 1.0;

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = false;
% sa(2) = sa(1);
target.setSurfaceAttributes(sa);
target.hitmethod = eHitMethod.Beam;

% pointing_accuracy = 0:0.25:2;
% laser_power = 5000:5000:40000;
pointing_accuracy = 2;
laser_power = 5000;

propagator.init(initstate);
propagator.addLaserStation(station);
    
%% Add Objects
for st=1
    if st==1
        station.bPulsed = true;
    else
        station.bPulsed = false;
    end
    for i = 1:length(pointing_accuracy)
        propagator.station.trackError = pointing_accuracy(i);
        for j = 1:length(laser_power)
            propagator.station.referencePower = laser_power(j);
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

            dim_describ = ["power bias variation ", num2str(i); ...
                "power simulated ", num2str(propagator.station.referencePower); ...
                "groundbased ", num2str(station.bGroundBased); ...
                "pulsed laser ", num2str(station.bPulsed)];

            if exist('eph', 'var')
                save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'dim_describ');
            else
                save(save2file, 'propagator', 'initstate', 'Y0', 'dim_describ');
            end
        end

    end
end
