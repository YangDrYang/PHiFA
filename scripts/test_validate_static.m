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
inputstr(1).resolution = 200;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

%% initialice target
cmpt = clCompoundTarget();

%% test segment
% inputstr(1).segobj = makeCube(0.1);
% inputstr(1).segobj = makePlate(0.1, 0.005);
% inputstr(1).segobj = makeEllipsoid(11, 0.05, 0.05, 0.05);
% inputstr(1).segobj = makeRandEllipsoid(9, 0.05, 0.05, 0.05,0.001);
inputstr(1).segobj = makeCylinder(9,0.025,0.1);
% inputstr(1).segobj = makeCone(9,0.1,0.25);

%% init
nsegments = length(inputstr);
cmpt.init(inputstr, tarinput);
%% set hit methode
% cmpt.hitmethod = eHitMethod.Net;
% cmpt.hitmethod = eHitMethod.Area;
cmpt.hitmethod = eHitMethod.Beam;


%% init Laserstation
last = clLaserStation();
last.referencePower = 400000; % \citep{battiston2017systematic}
last.repetitionRate = 10;
last.PulseLength = 5e-9;
% last.xv = [0 0 -1000 0 0 0];
last.beamType = eLaserBeam.Gaussian;
% last.bUsetophatIntensity = false;
% last.tophatDiameter = 2;
% last.tophatIntensity = 25000*last.repetitionRate;
last.trackError = 0;
last.pointError = 0;
time=0;

%% simulation

clear outputlog
outputlog = clOutputPulse();
deltat = 1/last.repetitionRate;
outputlog.init(1);
cmpt.xv = [0 0 0 0 0 0];
startxv = cmpt.xv;
cmpt.qw = [1 0 0 0 0 0 0];

wwdot = [quat2eul(quaternion(cmpt.qwdot(1:4))) QTForm(cmpt.qwdot(1:4).',cmpt.qwdot(5:7).').'];
outputlog.stepFinished(time,0,cmpt.xv,last.xv,wwdot,0)

snpstr = last.inducePulse(cmpt);
outputlog.ioSetNextPulse(snpstr);

[impulse, torque, gistru, hits] = cmpt.getLaserAccel(last);
outputlog.ioGetImpulse(impulse,torque,gistru);
cmpt.move(impulse, deltat);
cmpt.spin(torque, deltat);

wwdot = [quat2eul(quaternion(cmpt.qwdot(1:4))) QTForm(cmpt.qwdot(1:4).',cmpt.qwdot(5:7).').'];
outputlog.stepFinished(time,1,cmpt.xv,last.xv,wwdot,deltat);
%% plot
% plotSegmentWithHits(cmpt,hits);
plotSegment(cmpt);
