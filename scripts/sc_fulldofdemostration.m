%% load
load logfiles/sim0001.mat

time = eph(1,:)./60;
state = eph(2:4,:)./1000;
velo = eph(5:7,:)./1000;

timeperiod = (eph(1,end) - eph(1,1))/60;

ephsize = size(eph);
if ephsize(1)>7
    rows = 2;
    angles = (quat2eul(quaternion(norm_quat(eph(8:11,:))'))').*180/pi;
    arates = eph(12:14,:).*180/pi;   
else
    rows = 1;
end
%% print
figure;
subplot(1,2,1);
plot(time, angles(3,:),time, angles(2,:),time, angles(1,:));
xlim([0, 20]);
set(gca,'FontSize',16)
legend('x','y','z');
xlabel('Minutes [m]');
ylabel('Angle [deg]');
subplot(1,2,2);
plot(time, arates(1,:),time, arates(2,:),time, arates(3,:));
xlim([0, 20]);
set(gca,'FontSize',16)
legend('x','y','z');
xlabel('Minutes [m]');
ylabel('Angular velocity [deg/s]');