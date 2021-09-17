% simulate_meteorix.m

clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
% w0 = [0; 0.1*4.5/7; 0].*pi./180;
% q0 = eul2quat([50*-0.9/6; 0; 50*9.4/6].*pi./180);
w0 = [0; 0; 0].*pi./180;
q0 = eul2quat([0; 0; 0].*pi./180);
% a = 6878136.3;
% ecc = 0.002;
% i = 1.7000029;
% argp = 1.5707963;
% omega = 4.4918767;
% nu = 0;
% p = a*(1-ecc^2);
% [r,v] = coe2rv(p./1000, ecc, i, omega, argp, nu, 0, 0, 0);
% r0 = r*1000;
% v0 = v*1000;
% Y0 = [r0; v0; q0; w0];
Y0 = [6644264.68005270;-1955163.54312126;-842230.248328992;2292.78323879521;6199.46750663702;3667.02707298642; ...
    q0; w0];

initstate = clPropagatorInitialValues();
% initstate.proptime = 12000;
initstate.proptime = 100;
% initstate.year = 2010;
% initstate.mon = 1;
% initstate.day = 1;
% initstate.hour = 0;
% initstate.min = 0;
% initstate.sec = 0;
% initstate.mjd = Mjday(initstate.year,initstate.mon,initstate.day,...
%     initstate.hour,initstate.min,initstate.sec);
initstate.mjd = 52389.5801736113;
initstate.sec = -1;
initstate.step = 1;
initstate.dof = 3;
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
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets

initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects



%% create Target
target = loadDSPOSEMeteorix();
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0;
station.referencePower = 10000;
station.beamResolution = 200;
station.bOnlyInDark = true;
station.bGroundBased = false;
% station.lla = [0 0 400]; %altitude of 400km
[sma,inc,ecc,argp,raan,theta] = state2coe(Y0(1:3),Y0(4:6));
[r,v] = coe2rv(sma*(1-ecc^2)/1000+50, ecc, inc, raan, argp, theta, 0, 0, 0);
station.coor = [r;v]*1000;

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
plotSimulationOutput(initstate,propagator.output_table,'Meteorix');
