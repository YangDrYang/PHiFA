%% some initial values

Y0 = [ +7144843.808; +217687.110; -506463.296; +562.650611; -1616.516697; +7358.157263; ...
    1; 0; 0; 0; 0; 0; 0];

initstate = clPropagatorInitialValues();
initstate.proptime = 24*60*60;

%% Propagator
global propagator
propagator = clPropagator();

propagator.init(initstate);

%% create Objects

target = createCubeTarget();

station = createDefaultLaserstation();
station.init();
    
%% Add Objects

propagator.targetobjects = target;
propagator.addLaserStation(station);

%% Propagation Radau
% profile on
% eph = propagator.propagateTarget_radau(Y0,60);

%% Propagation Matlab's ODE
tic
profile on
options = odeset('RelTol',1e-11,'AbsTol',1e-6);
eph = propagator.propagateTarget_odex(@ode45,options,Y0,30);
toc