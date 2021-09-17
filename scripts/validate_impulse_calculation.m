%validate functionality

tarinput = clCompoundTargetInput;
inputstr = clCompoundSegmentsInput;

inputstr(1).bSTL = false;
inputstr(1).stlFile = "";

inputstr(1).surfAtr = clSA_const();
inputstr(1).surfAtr.init(25*10^(-6));
inputstr(1).surfAtr.name = 'testsurf1';
inputstr(1).offset = [0 0 0];
inputstr(1).name = 'test1';
inputstr(1).density = 2700;
inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(1).resolution = 100;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

%% initialice target
cmpt = clCompoundTarget();

%% test segments list

segobj(1) = makePlate(0.1, 0.005);
segobj(2) = makeEllipsoid(11, 0.05, 0.05, 0.05);

%% test segment
inputstr(1).segobj = makeCube(0.25);
% inputstr(1).segobj = makePlate(0.1, 0.005);
% inputstr(1).segobj = makeRandEllipsoid(7, 0.25, 0.25, 0.1, 0.01);
% inputstr(1).segobj = makeEllipsoid(11, 0.05, 0.05, 0.05);
% inputstr(1).segobj = makeCylinder(9,0.1,0.25);
% inputstr(1).segobj = makeCone(9,0.025,0.1);
% inputstr(1).segobj = segobj(1);

%% init
cmpt.init(inputstr, tarinput);
%% set hit methode
% cmpt.hitmethod = eHitMethod.Net;
% cmpt.hitmethod = eHitMethod.Area;
cmpt.hitmethod = eHitMethod.Beam;

%% init Laserstation
last = clLaserStation();
last.referencePower = 400000; % battiston2017systematic
last.repetitionRate = 10;
last.PulseLength = 5e-9;
last.xv = [0 0 1000 0 0 0];
last.beamType = eLaserBeam.TopHat;
last.bUsetophatIntensity = true;
last.tophatDiameter = 2;
last.beamResolution = 1000;
last.tophatIntensity = 25000*last.repetitionRate;
last.trackError = 0;
last.pointError = 0;

%% simulation
last.inducePulse(cmpt);
[laserimpulse, lasermomentum, laserforce, lasertorque, projectedArea] = ...
                            cmpt.getLaserAccel(last);

