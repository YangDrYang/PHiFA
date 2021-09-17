clear
clc

%% Create Propagator
% global propagator
propagator = clPropagator.instance();

%% some initial values
timestamp = now;
Y0 = [getInitialState(0,0,700); eul2quat([randn()*pi/2; randn()*pi/2; randn()*pi/2]); 0; 0; 0];
% Y0 = [ 40000000; 0; 0; 0; 0; 0; ...
%     1; 0; 0; 0; 0; 0; 0];
initstate = clPropagatorInitialValues();
initstate.proptime = 10;
initstate.step = 10;

initstate.a_grav = 0; % harmonic terms of gravity
initstate.g_grav = 0;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 40;
initstate.m = 40;

initstate.a_lase = 1; % laser engagements
initstate.g_lase = 1;
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

propagator.init(initstate);

%% create Target
target = createRandEllipsoidTarget(0.2, 3, 2, 30,false);
target.hitmethod = eHitMethod.Beam;
%% create Laserstation
station = createOHigginsLaserstation();
station.bPulsed = false;
station.trackError = 0;
station.pointError = 0;
station.referencePower = 10000;
station.bGroundBased = false;
station.lla = [0 0 600];

target.setSurfaceAttributes(...
    createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98));

%% Add Objects

propagator.addTarget(target);
propagator.addLaserStation(station);

%% Do Loop Over Attitude and Resolutions
profile on
% res = 10:5:100;
res = [1 5:5:105 110:10:300 325:25:450 500:50:600];
output_beam(1:length(res)) = clSimulationOutput();
output_net(1:length(res)) = clSimulationOutput();
% %% BEAM
for i = 1:length(res)
    tic;
    fprintf('###\tBeam resolution: %d\t###\n',res(i));
    propagator.station.beamResolution = res(i);
    propagator.sun.beamResolution = res(i);
    [~, output_beam(i)] = target.Accel(0,Y0);
    output_beam(i).aux1 = toc;
end
evalBeamConvergence(output_beam, res, propagator.AuxParam, 'Beam');
%% AREA
target.hitmethod = eHitMethod.Area;
tic;
fprintf('###\tArea method\t###\n');
[~, output_area] = target.Accel(0,Y0);
output_area.aux1 = toc;

%% NET
target.hitmethod = eHitMethod.Net;
for i = 1:length(res)
    tic;
    fprintf('###\tNet resolution: %d\t###\n',res(i));
    target.setResolution(res(i));
    [~, output_net(i)] = target.Accel(0,Y0);
    output_net(i).aux1 = toc;
end
evalBeamConvergence(output_net, res, propagator.AuxParam, 'Net');
%% EVAL
evalHitMethods(output_beam, output_net, output_area, res, propagator.AuxParam); 
prof_info = profile('info');
profile viewer

%% EVAL again
% evalBeamConvergence(output_beam, res, propagator.AuxParam);
% evalBeamConvergence(output_net, res, propagator.AuxParam);
% evalHitMethods(output_beam, output_net, output_area, res, propagator.AuxParam);