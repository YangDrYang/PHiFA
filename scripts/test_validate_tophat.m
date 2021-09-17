%validate functionality

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
% inputstr(1).segobj = makeCube(0.25);
% inputstr(1).segobj = makePlate(0.1, 0.005);
% inputstr(1).segobj = makeRandEllipsoid(7, 0.25, 0.25, 0.1, 0.01);
inputstr(1).segobj = makeEllipsoid(11, 0.05, 0.05, 0.05);
% inputstr(1).segobj = makeCylinder(9,0.1,0.25);
% inputstr(1).segobj = makeCone(9,0.025,0.1);

%% init
nsegments = length(inputstr);
cmpt.init(nsegments, inputstr);
%% set hit methode
cmpt.hitmethod = eHitMethod.Net;
% cmpt.hitmethod = eHitMethod.Area;

%% init Laserstation
last = clLaserStation();
last.referencePower = 400000; % battiston2017systematic
last.repetitionRate = 10;
last.PulseLength = 5e-9;
last.xv = [0 -1000 -1000 0 0 0];
last.beamType = eLaserBeam.TopHat;
last.bUsetophatIntensity = true;
last.tophatDiameter = 2;
last.tophatIntensity = 25000*last.repetitionRate;
last.trackError = 0;
last.pointError = 0;

%% simulation

clear outputlog
outputlog = clOutputPulse();
deltat = 1/last.repetitionRate;
time = 0;
nSteps = 500;
outputlog.init(nSteps);
cmpt.xv = [0 0 0 0 0 0];
startxv = cmpt.xv;
cmpt.qwdot = [1 0 0 0 0 0 0];

wwdot = [quat2eul(quaternion(cmpt.qwdot(1:4))) QTForm(cmpt.qwdot(1:4).',cmpt.qwdot(5:7).').'];
outputlog.stepFinished(time,0,cmpt.xv,last.xv,wwdot,0)

for i = 1:nSteps
    fprintf('##### Step %03i of %03i #####\n%s\n',i,nSteps,datestr(datetime('now')));

    time=time+deltat;
    
    snpstr = last.setNextPulse(cmpt);
    outputlog.ioSetNextPulse(snpstr);
    
    [impulse, torque, gistru, hits] = cmpt.getImpulse(last, last.raySAT);
    outputlog.ioGetImpulse(impulse,torque,gistru);
    cmpt.move(impulse, deltat);
    cmpt.spin(torque, deltat);
    
    wwdot = [quat2eul(quaternion(cmpt.qwdot(1:4))) QTForm(cmpt.qwdot(1:4).',cmpt.qwdot(5:7).').'];
    outputlog.stepFinished(time,i,cmpt.xv,last.xv,wwdot,deltat);
    
    if mod(i,nSteps/4) == 0
        plotSegmentWithHits(cmpt,hits);
    end
    
end

outputlog.plotMovement();
