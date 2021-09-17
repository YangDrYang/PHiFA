%test script

clear;

%% init Target
inputstr = clCompoundSegmentsInput;

% from stl
% inputstr(1).stlFile = string('stlfiles\tetrahedron.stl');
% inputstr.stlFile = string('stlfiles\femur.stl');
% inputstr(1).bSTL = true;
%from function
testseg = makeCube(0.1);
inputstr(1).segobj = testseg;
inputstr(1).bSTL = false;

inputstr(1).surfAtr = clSA_wilken2015();
a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
t = [15.04 1.19 2.41 23.38];
inputstr(1).surfAtr.init(a, t);
inputstr(1).surfAtr.name = 'testsurf1';
inputstr(1).offset = [0 0 0];
inputstr(1).name = 'test1';
inputstr(1).density = 1000;
inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(1).resolution = 100;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

cmpt = clCompoundTarget();
nSegments = length(inputstr);
cmpt.init(nSegments, inputstr);

%%
cmpt.hitmethod = eHitMethod.Net;
% cmpt.hitmethod = eHitMethod.Area;

%% init Laserstation
last = clLaserStation();
last.referencePower = 400000; % battiston2017systematic
last.repetitionRate = 5;
last.PulseLength = 5e-9;
last.Wavelength = last.Wavelength;
last.xv = [0 6371*sin(pi/4) 6371*cos(pi/4) 0 0 0];


%% Simulation
clear tmov tspi hmom imp lmov lmom fom poe mcm maoi mint marea mimp mamom
deltat = 1/last.repetitionRate;
nSteps = 100;
tmov(1:nSteps+1,1:6) = 0;
tspi(1:nSteps+1,1:6) = 0;
hmom(1:nSteps+1,1:3) = 0;
imp(1:nSteps+1,1:6) = 0;
lmov(1:nSteps+1,1:6) = 0;
lmom(1:nSteps+1,1:3) = 0;
fom(1:nSteps) = 0;
poe(1:nSteps) = 0;
mcm(1:nSteps) = 0;
maoi(1:nSteps) = 0;
mint(1:nSteps) = 0;
marea(1:nSteps) = 0;
mimp(1:nSteps) = 0;
mamom(1:nSteps) = 0;
cmpt.xv = [0 0 400000 0 0 0];
startxv = cmpt.xv;
cmpt.qwdot = [1 0 0 0 0 0 0];
imp(1:nSteps,1:6) = 0;

% dispstat('','init');
for i = 1:nSteps
%     dispstat(sprintf('Simulating %d%%',round(i/nSteps*100)));
    fprintf('##### Step %03i of %03i #####\n%s\n',i,nSteps,datestr(datetime('now')));

    tmov(i,:) = cmpt.xv - startxv;
    tspi(i,1:3) = quat2eul(quaternion(cmpt.qwdot(1:4)));
    tspi(i,4:6) = QTForm(cmpt.qwdot(1:4).',cmpt.qwdot(5:7).');
    hmom(i,1:3) = QTForm(cmpt.qwdot(1:4).',cmpt.moi*cmpt.qwdot(5:7).');
    lmom(i,1:3) = cmpt.xv(4:6)*cmpt.mass;
    
    snpstr = last.setNextPulse(cmpt);
    poe(i) = snpstr.poe;
    
    [impulse, torque, gistru] = cmpt.getImpulse(last, last.raySAT);
    imp(i+1,:) = [impulse torque];
    cmpt.move(impulse, deltat);
    cmpt.spin(torque, deltat);
    fom(i) = gistru.fom;
    mcm(i) = gistru.mcm;
    maoi(i) = gistru.maoi;
    mint(i) = gistru.mint;
    marea(i) = gistru.marea;
    mimp(i) = gistru.mimp;
    mamom(i) = gistru.mamom;
end

tmov(nSteps+1,:) = cmpt.xv - startxv;
tspi(nSteps+1,1:3) = quat2eul(quaternion(cmpt.qwdot(1:4)));
tspi(nSteps+1,4:6) = cmpt.qwdot(5:7);
hmom(nSteps+1,1:3) = QTForm(cmpt.qwdot(1:4).',cmpt.moi*cmpt.qwdot(5:7).').';
lmom(nSteps+1,1:3) = cmpt.xv(4:6)*cmpt.mass;

%% Plot
close all;
time = 0:deltat:nSteps*deltat;

subplot(2,2,1);
plot(time, tmov(:,1),time, tmov(:,2),time, tmov(:,3));
title('Position');
legend({'x','y','z'});
xlabel('time [s]');
ylabel('distance [m]');
grid;

subplot(2,2,2);
plot(time, tmov(:,4),time, tmov(:,5),time, tmov(:,6));
title('Velocity');
legend({'x','y','z'});
xlabel('time [s]');
ylabel('velocity [m/s]');
grid;

% subplot(2,3,3);
% plot(time, lmom(:,1),time, lmom(:,2),time, lmom(:,3));
% title('Linear Momentum');
% legend({'x','y','z'});
% xlabel('time [s]');
% ylabel('Momentum [Nm]');

subplot(2,2,3);
yyaxis left;
plot(time, tspi(:,3),time, tspi(:,2),time, tspi(:,1));
title('Rotation');
xlabel('time [s]');
ylabel('angle [rad]');
yyaxis right;
plot(time,[0 fom], '+');
ylabel('Radial/Axial Impulse [-]');
legend({'x','y','z','FOM'});
grid;

subplot(2,2,4);
plot(time, tspi(:,4),time, tspi(:,5),time, tspi(:,6));
title('Anglular Rate');
legend({'x','y','z'});
xlabel('time [s]');
ylabel('anglular rate [rad/s]');
grid;
% subplot(2,3,6);
% plot(time, hmom(:,1),time, hmom(:,2),time, hmom(:,3));
% title('Angular Momentum');
% legend({'x','y','z'});
% xlabel('time [s]');
% ylabel('Momentum [Nm]');

% saveas(gcf, 'figures\test_getimpulse.png');
savefig('figures\test_getimpulseplusbeamprop_xv.fig');

time = time(1:end-1);
figure;
yyaxis left;
plot(time,fom);
title('Figure of merit');
xlabel('time [s]');
ylabel('Radial/Axial Impulse [-]');
yyaxis right;
plot(time,poe)
ylabel('Pointing error [murad]');
legend('Figure of Merit','Pointing error');
savefig('figures\test_getimpulseplusbeamprop_fompoe.fig');
