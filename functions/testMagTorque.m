% simulate_meteorix.m

clear all
clf
clc
close all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 0.3].*pi./180;
q0 = eul2quat([0, 50, 0]'.*pi./180);
epoch = Mjday(2017,8,22,1,0,0);
orb_iss = initialisation_tle('./inputfiles/iss_tle.txt',epoch);
[sma,inc,ecc,argp,raan,theta] = state2coe(orb_iss(1:3),orb_iss(4:6));
% % [r,v] = coe2rv(sma*(1-ecc^2)/1000+10, ecc, inc+1.0/180*pi, raan+1.0/180*pi, argp+1.0/180*pi, theta+1.0/180*pi, 0, 0, 0);
% [r,v] = coe2rv(sma*(1-ecc^2)/1000+100, ecc, inc, raan, argp, theta+1.0/180*pi, 0, 0, 0);
% orb_target = [r;v]*1000;
koe.sma = sma+100*1e3;koe.incl = inc;koe.ecc = ecc; koe.argp = argp; koe.raan = raan;koe.tran = theta+1.0/180*pi;
orb_target = kep2cart(koe);
mu = 3.986004418e14;
period = 2*pi*sqrt(koe.sma^3/mu);
% orb_target = initialisation_tle('./inputfiles/lemur2_tle.txt',epoch);
Y0 = [orb_target; q0; w0];
Diff = norm(orb_target - orb_iss);

initstate = clPropagatorInitialValues();
initstate.proptime = 7200;
% initstate.proptime = 100;
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
initstate.dof = 3;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 0;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 100;
initstate.m = 100;

initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 1; % drag
initstate.g_drag = 0;
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 0;
initstate.a_erad = 1; % earth reflected and emmited radiation pressure: albedo and infrared
initstate.g_erad = 0;
initstate.g_mag = 0; % magnetic torque
initstate.a_sun = 1; % gravity of sun
initstate.a_moon = 1; % gravity of moon
initstate.a_planets = 1; % gravity of other planets
initstate.a_solidEarthTides = 1; % solid earth tides
initstate.a_oceanTides = 1; % ocean tides
initstate.a_relativity = 1; % relativity effects
initstate.albedo_model = eAlbedoModel.Stephens; % Albedo model

initstate.b_logging = 1; % else takes ages
initstate.n_log_only_every_xth_step = 1;

%% create Target
target = loadDSPOSEMeteorix();
% target = loadDSPOSEMeteorixLostPanel();
% plotTarget(target);
%% create Laserstation
station = createISSLaserstation();
station.repetitionRate = 1;
station.trackError = 3;
station.pointError = 1;
% station.trackError = 0;
% station.pointError = 0;
station.referencePower = 10000;
station.beamResolution = 100;
station.beamType = eLaserBeam.Gaussian;
% station.beamType = eLaserBeam.TopHat;
station.bOnlyInDark = true;
station.bGroundBased = false;
if ~station.bGroundBased % if space-borne, to set for orbit integration step
    station.intstep = 1.0;
end
station.coor = orb_iss;
station.minHitProbability = 0.05;
% station.coor = initialisation_tle('./inputfiles/iss_tle.txt',epoch);

%% some more target conf (needs knowledge of pulselength, wavelength)
target.hitmethod = eHitMethod.Beam;
% target.hitmethod = eHitMethod.Area;

propagator.init(initstate);
propagator.station.bPulsed = true;    
propagator.station.trackError = 0;
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
% odefcn = @ABM8;
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
plotSimulationOutput(initstate,propagator.output_table,'Meteorix');
