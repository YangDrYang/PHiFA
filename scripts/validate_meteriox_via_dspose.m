clear all

addpath('../TianGong1');
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 10].*pi./180;
q0 = angle2quat(0*pi/180, 50*pi/180, 0*pi/180,'ZYX')';
yy = 2016;
mon = 4;
dd = 24;
hr = 21;
mm = 55;
sec = 28;
epoch = Mjday(yy, mon, dd, hr, mm, sec);
orb_iss = initialisation_tle('./inputfiles/iss_tle.txt',epoch);
[sma,inc,ecc,argp,raan,theta] = state2coe(orb_iss(1:3),orb_iss(4:6));
[r,v] = coe2rv(sma*(1-ecc^2)/1000+100, ecc, inc, raan, argp, theta+1.0/180*pi, 0, 0, 0);
orb_target = [r;v]*1000;
Y0_teme = [orb_target; q0; w0];
mat = mat_teme2eci(yy, mon, dd, hr, mm, sec);%teme to j2000
Y0 = [mat*orb_target(1:3); mat*orb_target(4:6); q0; w0];
initstate = clPropagatorInitialValues();
initstate.proptime = 3600*12;
initstate.step = 1;

initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.n = 100;
initstate.m = 100;
initstate.a_drag = 0; % drag
initstate.g_drag = 0;
initstate.draggradient = 0;%only for area method
initstate.a_srad = 0; % solar radiation pressure
initstate.g_srad = 0;
initstate.a_erad = 0; % earth reflected and emmited radiation pressure: albedo and infrared
initstate.g_erad = 0;
initstate.g_mag = 0; % magnetic torque
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets
initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects
initstate.albedo_model = eAlbedoModel.CERES; % Albedo model

%% create Target
target = loadDSPOSEMeteorix();
% target.hitmethod = eHitMethod.Beam;
target.hitmethod = eHitMethod.Area;

propagator.init(initstate);

%% create Laserstation
station = createOHigginsLaserstation();

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
sa(3) = sa(1);
sa(4) = sa(1);
sa(5) = sa(1);
target.setSurfaceAttributes(sa);
% target.hitmethod = eHitMethod.Beam;


%% Propagation
for i = 1:6
    
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
    elseif i==4
        initstate.a_grav = 1; % harmonic terms of gravity
        initstate.g_grav = 1;
        initstate.a_drag = 0; % drag
        initstate.g_drag = 0;
        initstate.a_srad = 1; % solar radiation pressure
        initstate.g_srad = 1;
    elseif i==5
        initstate.a_grav = 0; % harmonic terms of gravity
        initstate.g_grav = 0;
        initstate.a_drag = 1; % drag
        initstate.g_drag = 1;
        initstate.a_srad = 1; % solar radiation pressure
        initstate.g_srad = 1;
    elseif i==6
        initstate.a_grav = 1; % harmonic terms of gravity
        initstate.g_grav = 1;
        initstate.a_drag = 1; % drag
        initstate.g_drag = 1;
        initstate.a_srad = 1; % solar radiation pressure
        initstate.g_srad = 1;
    end
    
    propagator.init(initstate);
    propagator.addTarget(target);
%     propagator.addLaserStation(station);  
    
    %% Start Diary
    propagator.initNewSim();

    %% Propagation
    tic
    profile on
    odefcn = @RK45; 
    eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
    t_phifa = toc;
    
    cd dspose
    target2DSPOSE(target, Y0_teme, propagator.AuxParam, '');  
    !./dspose_exec
    cd ..
    dspose_timestamp = getLatestDSPOSETimestamp();
    t_dspose = toc;
    
    %% Write Output
    save2file = sprintf('%s.mat', propagator.outfilename);
    if exist('eph', 'var')
        save(save2file, 'propagator', 'initstate', 'eph', 'Y0', 'dspose_timestamp');
    else
        save(save2file, 'propagator', 'initstate', 'Y0', 'dspose_timestamp');
    end    
    
    close all
    addpath('../Tiangong1');
    tmp = char(propagator.outfilename);
    [rotOrbAng_DSPOSE,rotOrbAng_PHiFA,diffRotOrbAng] = compareSimulationWithDSPOSE(epoch, dspose_timestamp, str2double(tmp(end-3:end)));
    
end


