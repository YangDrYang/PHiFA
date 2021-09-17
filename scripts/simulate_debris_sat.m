% simulate_meteorix.m

clear all
% clf
% close all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
% w0 = [0; 0.1*4.5/7; 0].*pi./180;
% q0 = eul2quat([50*-0.9/6; 0; 50*9.4/6].*pi./180);
w0 = [0; 0; 0].*pi./180;
q0 = eul2quat([0; 0; 0].*pi./180);
epoch = Mjday(2010,1,2,0,0,0);
[r,v] = coe2rv(6885,0.001/180*pi,80/180*pi,180/180*pi,120/180*pi,0,0,0,0);
orb_sat = [r;v]*1000;
[r,v] = coe2rv(7170,0.001,100/180*pi,0,180/180*pi,220/180*pi,0,0,0);
orb_target = [r;v]*1000;
% orb_target = initialisation_tle('./inputfiles/lemur2_tle.txt',epoch);
Y0 = [orb_target; q0; w0];
Diff = norm(orb_target - orb_sat);

initstate = clPropagatorInitialValues();
% initstate.proptime = 12000;
initstate.proptime = 150;
% initstate.year = 2010;
% initstate.mon = 1;
% initstate.day = 1;
% initstate.hour = 0;
% initstate.min = 0;
% initstate.sec = 0;
% initstate.mjd = Mjday(initstate.year,initstate.mon,initstate.day,...
%     initstate.hour,initstate.min,initstate.sec);
initstate.mjd = epoch;
initstate.sec = -1;
initstate.step = 1;
initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 100;
initstate.m = 100;

initstate.a_lase = 1; % laser engagements
initstate.g_lase = 1;
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



%% create Target
target = loadDSPOSEMeteorix();
% plotTarget(target);
%% create Laserstation
station = createISSLaserstation();
station.trackError = 3;
station.pointError = 1;
station.referencePower = 10000;
station.beamResolution = 200;
station.bOnlyInDark = true;
station.bGroundBased = false;
station.coor = orb_sat;
station.minHitProbability = 0.05;
% station.coor = initialisation_tle('./inputfiles/iss_tle.txt',epoch);

%% some more target conf (needs knowledge of pulselength, wavelength)
target.hitmethod = eHitMethod.Beam;

propagator.init(initstate);
propagator.station.bPulsed = true;    
propagator.station.trackError = 1.5;
propagator.station.referencePower = 20000;
propagator.station.PulseLength = 1E-9;
propagator.station.Wavelength = 1064E-9;

sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
sa(3) = sa(1);
sa(4) = sa(1);
sa(5) = sa(1);
target.setSurfaceAttributes(sa);
%% Add Objects
propagator.addLaserStation(station);

propagator.initNewSim();
propagator.addTarget(target);

%% Propagation
tic
%         profile on
odefcn = @RK45; 
% odefcn = @ode23;
eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
toc

%% Write Output
propagator.finishSim();
save2file = sprintf('%s.mat', propagator.outfilename);

if exist('eph', 'var')
    save(save2file, 'propagator', 'initstate', 'eph', 'Y0');
else
    save(save2file, 'propagator', 'initstate', 'Y0');
end

plot3Ephemerides_ecef(eph,station,initstate)
plotSimulationOutput(initstate,propagator.output_table,'debris');
