sd = loadSimData(1,4);

for i=1:length(sd)
    plot3Ephemerides(sd{i}.propagator.station.ephemerides, sd{i}.eph);
end

%% position
figure('Name', 'Position');
subplot(3,1,1);
for i = 1:length(sd)
    h = plot(sd{i}.eph(1,:)./3600, sd{i}.eph(2,:));
    set(h,'DisplayName',[num2str(sd{i}.propagator.station.referencePower*...
        (1+sd{i}.propagator.station.powerBias/100)/1000) ' kW']);
    hold on;
end
legend();
xlabel('Time [h]');
ylabel('dx [m]');
grid;
hold off;

subplot(3,1,2);
for i = 1:length(sd)
    h = plot(sd{i}.eph(1,:)./3600, sd{i}.eph(3,:));
    set(h,'DisplayName',[num2str(sd{i}.propagator.station.referencePower*...
        (1+sd{i}.propagator.station.powerBias/100)/1000) ' kW']);
    hold on;
end
legend();
xlabel('Time [h]');
ylabel('dy [m]');
grid;
hold off;

subplot(3,1,3);
for i = 1:length(sd)
    h = plot(sd{i}.eph(1,:)./3600, sd{i}.eph(4,:));
    set(h,'DisplayName',[num2str(sd{i}.propagator.station.referencePower*...
        (1+sd{i}.propagator.station.powerBias/100)/1000) ' kW']);
    hold on;
end
legend();
xlabel('Time [h]');
ylabel('dz [m]');
grid;
hold off;

%% velocity
figure('Name', 'Velocity');
subplot(3,1,1);
for i = 1:length(sd)
    h = plot(sd{i}.eph(1,:)./3600, sd{i}.eph(5,:));
    set(h,'DisplayName',[num2str(sd{i}.propagator.station.referencePower*...
        (1+sd{i}.propagator.station.powerBias/100)/1000) ' kW']);
    hold on;
end
legend();
xlabel('Time [h]');
ylabel('dx [m/s]');
grid;
hold off;

subplot(3,1,2);
for i = 1:length(sd)
    h = plot(sd{i}.eph(1,:)./3600, sd{i}.eph(6,:));
    set(h,'DisplayName',[num2str(sd{i}.propagator.station.referencePower*...
        (1+sd{i}.propagator.station.powerBias/100)/1000) ' kW']);
    hold on;
end
legend();
xlabel('Time [h]');
ylabel('dy [m/s]');
grid;
hold off;

subplot(3,1,3);
for i = 1:length(sd)
    h = plot(sd{i}.eph(1,:)./3600, sd{i}.eph(7,:));
    set(h,'DisplayName',[num2str(sd{i}.propagator.station.referencePower*...
        (1+sd{i}.propagator.station.powerBias/100)/1000) ' kW']);
    hold on;
end
legend();
xlabel('Time [h]');
ylabel('dz [m/s]');
grid;
hold off;