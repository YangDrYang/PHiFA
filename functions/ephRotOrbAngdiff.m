function [rotOrbAng_DSPOSE,rotOrbAng_PHiFA,diffRotOrbAng] = ephRotOrbAngdiff(varargin)

nEpoch = length(varargin{1}(1,:));
rotOrbAng_DSPOSE = zeros(nEpoch,6);
rotOrbAng_PHiFA = zeros(nEpoch,6);
diffRotOrbAng = zeros(nEpoch,6);
iner = varargin{3};
for i = 1:nEpoch
    [~,~,inc,Omega,~,~]=kepel(varargin{1}(2:4,i),varargin{1}(5:7,i));
    rotm = angle2dcm(varargin{1}(8,i)/180*pi,varargin{1}(9,i)/180*pi,varargin{1}(10,i)/180*pi,'ZYX');
    [ang1,ang2,ang3] = dcm2angle(rotm,'ZYZ');
    EulerAng = [ang1,ang2,ang3]';
    EulerAngVel = rotm*varargin{1}(11:13,i)./180*pi;
    rotOrbAng_DSPOSE(i,:) = Euler2RotOrb(EulerAng,EulerAngVel,Omega,inc,iner);
%     EulerAng = SpinCalc('EA321toEA313',varargin{1}(8:10,i)',1e-6,0)';%from 311 Euler angles to 313
%     EulerAngVel = SpinCalc('EA321toEA313',varargin{1}(11:13,i)',1e-6,0)';%from 311 Euler angular velocity to 313
%     rotOrbAng_DSPOSE(i,:) = Euler2RotOrb(EulerAng/180*pi,EulerAngVel/180*pi,Omega,inc,iner);
    [~,~,inc,Omega,~,~]=kepel(varargin{2}(2:4,i),varargin{2}(5:7,i));
    rotm = angle2dcm(varargin{2}(8,i)/180*pi,varargin{2}(9,i)/180*pi,varargin{2}(10,i)/180*pi,'ZYX');
    [ang1,ang2,ang3] = dcm2angle(rotm,'ZYZ');
    EulerAng = [ang1,ang2,ang3]';    
    EulerAngVel = rotm*varargin{2}(11:13,i)./180*pi;
    rotOrbAng_PHiFA(i,:) = Euler2RotOrb(EulerAng,EulerAngVel,Omega,inc,iner);
%     EulerAng = SpinCalc('EA321toEA313',varargin{2}(8:10,i)',1e-6,0)';%from 311 Euler angles to 313
%     EulerAngVel = SpinCalc('EA321toEA313',varargin{2}(11:13,i)',1e-6,0)';%from 311 Euler angular velocity to 313
%     rotOrbAng_PHiFA(i,:) = Euler2RotOrb(EulerAng/180*pi,EulerAngVel/180*pi,Omega,inc,iner);
    diffRotOrbAng(i,1) = rotOrbAng_DSPOSE(i,1) - rotOrbAng_PHiFA(i,1);
    diffRotOrbAng(i,2:6) = angdiff(rotOrbAng_DSPOSE(i,2:6),rotOrbAng_PHiFA(i,2:6))*180/pi;%radian -> degree
%     diffRotOrbAng(i,2:6) = mod(rotOrbAng_DSPOSE(i,2:6) - rotOrbAng_PHiFA(i,2:6), 2*pi)*180/pi;%radian -> degree
end

figure('Name', 'Rotational Parameters')
subplot(3,2,1)
plot(varargin{1}(1,:),rotOrbAng_PHiFA(:,1))
% xlabel('Time [sec]','Fontsize',14);
ylabel({'Magnitude of Angular' ' Momentum [kg m^2 s^{-1}]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,2)
plot(varargin{1}(1,:),rotOrbAng_PHiFA(:,2).*(180/pi))
ylabel({'\psi_H [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,3)
plot(varargin{1}(1,:),rotOrbAng_PHiFA(:,3).*(180/pi))
xlabel('Time [sec]','Fontsize',14);
ylabel({'\theta_H [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,4)
plot(varargin{1}(1,:),rotOrbAng_PHiFA(:,4).*(180/pi))
ylabel({'\phi_H [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,5)
plot(varargin{1}(1,:),rotOrbAng_PHiFA(:,5).*(180/pi))
ylabel({'\theta^{\prime} [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,6)
plot(varargin{1}(1,:),rotOrbAng_PHiFA(:,6).*(180/pi))
xlabel('Time [sec]','Fontsize',14);
ylabel({'\phi^{\prime} [deg]'},'Fontsize',14)
xlim([0 nEpoch])


figure('Name', 'Difference of Rotational Parameters')
subplot(3,2,1)
plot(varargin{1}(1,:),diffRotOrbAng(:,1))
% xlabel('Time [sec]','Fontsize',14);
% ylabel({'Magnitude of Angular' ' Momentum [kg m^2 s^{-1}]'},'Fontsize',14)
ylabel({'H [kg m^2 s^{-1}]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,2)
plot(varargin{1}(1,:),diffRotOrbAng(:,2))
ylabel({'\psi_H [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,3)
plot(varargin{1}(1,:),diffRotOrbAng(:,3))
xlabel('Time [sec]','Fontsize',14);
ylabel({'\theta_H [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,4)
plot(varargin{1}(1,:),diffRotOrbAng(:,4))
ylabel({'\phi_H [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,5)
plot(varargin{1}(1,:),diffRotOrbAng(:,5))
ylabel({'\theta^{\prime} [deg]'},'Fontsize',14)
xlim([0 nEpoch])
subplot(3,2,6)
plot(varargin{1}(1,:),diffRotOrbAng(:,6))
xlabel('Time [sec]','Fontsize',14);
ylabel({'\phi^{\prime} [deg]'},'Fontsize',14)
xlim([0 nEpoch])