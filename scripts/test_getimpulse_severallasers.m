%test script

clear;

%% init Target
inputstr = clCompoundSegmentsInput;

inputstr(1).stlFile = string('stlfiles\tetrahedron.stl');
% inputstr.stlFile = string('stlfiles\femur.stl');

inputstr(1).surfAtr = clSurfaceAttributes();
a = [-0.01314 0.01753 0.01904 449.81 632.06 913.51];
t = [15.04 1.19 2.41 23.38];
inputstr(1).surfAtr.constantCouplingCoefficient = 20;
inputstr(1).surfAtr.init(a, t);
inputstr(1).surfAtr.name = 'testsurf1';
inputstr(1).offset = [0 0 0];
inputstr(1).name = 'test1';
inputstr(1).density = 2500;
inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
inputstr(1).resolution = 500;
inputstr(1).solid = true;
inputstr(1).twosided = false;
inputstr(1).thickness = 0.001;

cmpt = clCompoundTarget();
nSegments = length(inputstr);
cmpt.init(nSegments, inputstr);

cmpt.hitmethod = eHitMethod.Area;

%% init Laserstation
last(1:6) = clLaserStation();
nlast(1:6) = 0;
nlasers = 6;
for i = 1:nlasers
    last(i) = clLaserStation();
    last(i).referenceDiameter = 1;
    last(i).repetitionRate = 10;
    last(i).PulseLength = 5e-9;
    last(i).Wavelength = 2 * last(i).Wavelength;
end
last(1).xv = [100000000000 0 0 0 0 0];
last(2).xv = [-100000000000 0 0 0 0 0];
last(3).xv = [0 100000000000 0 0 0 0];
last(4).xv = [0 -100000000000 0 0 0 0];
last(5).xv = [0 0 100000000000 0 0 0];
last(6).xv = [0 0 -100000000000 0 0 0];

%% Simulation
deltat = 1/last(1).repetitionRate;
nSteps = 10000;
tmov(1:nSteps+1,1:6) = 0;
imp(1:nSteps+1,1:6) = 0;
% lmov1(1:nSteps+1,1:6) = 0;
% lmov2(1:nSteps+1,1:6) = 0;
ray = clRay();
cmpt.xv = [0 0 0 0 0 0];

dist(1:nSteps,1:6) = 0;

dispstat('','init');
for i = 1:nSteps
    dispstat(sprintf('Simulating %d%%',round(i/nSteps*100)));
    tmov(i,:) = cmpt.xv;
%     lmov1(i,:) = last1.xv;
%     lmov2(i,:) = last2.xv;
    for j = 1:6
        dist(i,j) = norm(cmpt.xv(1:3) - last(j).xv(1:3));
    end
    [~, ind] = min(dist(i,:));
%     ind = [1 2 3 4 5 6];
    impulse = [0 0 0];
    torque = [0 0 0];
    for k = 1:length(ind)
        ray.origin = last(ind(k)).xv(1:3);
        ray.direction = cmpt.xv(1:3) - last(ind(k)).xv(1:3);
        ray.direction = ray.direction/norm(ray.direction);
        [tmpimp, tmptor] = cmpt.getImpulse(last(ind(k)), ray);
        impulse = impulse + tmpimp;
        torque = torque + tmptor;
        nlast(ind(k)) = nlast(ind(k)) + 1;
    end
%     if dist1 < dist2
%         nlast1 = nlast1 + 1;
%         ray.origin = last1.xv(1:3);
%         ray.direction = cmpt.xv(1:3) - last1.xv(1:3);
%         ray.direction = ray.direction/norm(ray.direction);
%         [impulse, torque] = cmpt.getImpulse(last1, ray);
%     elseif dist2 < dist1
%         nlast2 = nlast2 + 1;
%         ray.origin = last2.xv(1:3);
%         ray.direction = cmpt.xv(1:3) - last2.xv(1:3);
%         ray.direction = ray.direction/norm(ray.direction);
%         [impulse, torque] = cmpt.getImpulse(last2, ray);
%     else
%         ray.origin = last1.xv(1:3);
%         ray.direction = cmpt.xv(1:3) - last1.xv(1:3);
%         ray.direction = ray.direction/norm(ray.direction);
%         [impulse1, torque1] = cmpt.getImpulse(last1, ray);
%         ray.origin = last2.xv(1:3);
%         ray.direction = cmpt.xv(1:3) - last2.xv(1:3);
%         ray.direction = ray.direction/norm(ray.direction);
%         [impulse2, torque2] = cmpt.getImpulse(last2, ray);
%         impulse = impulse1 + impulse2;
%         torque = torque1 + torque2;
%     end
    imp(i,:) = [impulse torque];
    cmpt.move(impulse, deltat);
%     last1.xv = [ (last1.xv(1:3) + deltat*last1.xv(4:6)) (last1.xv(4:6)) ];
%     last2.xv = [ (last2.xv(1:3) + deltat*last2.xv(4:6)) (last2.xv(4:6)) ];
%     last1.xv(2:3) = cmpt.xv(2:3);
%     last2.xv(2:3) = cmpt.xv(2:3);
end

tmov(nSteps+1,:) = cmpt.xv;

%% Plot
time = 0:deltat:nSteps*deltat;
% figure;
% plot(time, tmov(:,1),time, lmov1(:,1),time, lmov2(:,1));
% plot(time, tmov(:,1),time, tmov(:,2),time, tmov(:,3));
% plot(time, tmov(:,1));

% figure;
% hs(1) = subplot(2,1,1);
% hs(2) = subplot(2,1,2);
% scatter3(hs(1),lmov1(:,1),lmov1(:,2),lmov1(:,3),2,time);
% scatter3(hs(2),tmov(:,1),tmov(:,2),tmov(:,3),2,time);

figure;
scatter3(tmov(:,1),tmov(:,2),tmov(:,3),2,time);
xlabel('X');
ylabel('Y');
zlabel('Z');

% figure;
% % plot(time, imp(:,1),time, imp(:,2),time, imp(:,3));
% plot(time, imp(:,3));


