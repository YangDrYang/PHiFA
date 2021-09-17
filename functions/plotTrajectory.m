function plotTrajectory(eph,propagator)

%plot position components and relative distance of the laser station and target


laserEph = eph(1:4,:);
targetEph = propagator.station.ephemerides;
targetEph = targetEph(:,31:end-30);
relDis = NaN(1,length(laserEph(1,:)));
for i = 1:length(laserEph(1,:))
    relDis(i) = norm(laserEph(2:4,i) - targetEph(2:4,i));
end
   

figure
subplot(4,1,1)
plot(laserEph(1,:),laserEph(2,:));
hold on
plot(targetEph(1,:),targetEph(2,:));


subplot(4,1,2)
plot(laserEph(1,:),laserEph(3,:));
hold on
plot(targetEph(1,:),targetEph(3,:));

subplot(4,1,3)
plot(laserEph(1,:),laserEph(4,:));
hold on
plot(targetEph(1,:),targetEph(4,:));

subplot(4,1,4)
plot(laserEph(1,:),relDis(1,:));
