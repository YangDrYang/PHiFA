sd = loadSimData(180,3);

[Y0, mjd0, date] = loadEnvisatMahootiInitialState();
true_eph_ecef = loadEnvisatMahootiTrueEphemerides();
true_eph_eci = ecef2eci_ephemerides(true_eph_ecef, mjd0);

dd = zeros(7, length(true_eph_eci), length(sd));

for i = 1:length(sd)
    dd(1,:,i) = sd{i}.eph(1,:);
    dd(2:7,:,i) = true_eph_eci(2:7,:) - sd{i}.eph(2:7,:);
end

true_eph_eci(1,:) = true_eph_eci(1,:)./60;

figure('Name', 'Envisat Propagation');
subplot(3,1,1);
for i = 1:size(dd,3)
    plot(true_eph_eci(1,:), dd(2,:,i));
    hold on;
end

xlabel('Time [min]');
ylabel('dx [m]');
grid;
hold off;

subplot(3,1,2);
for i = 1:size(dd,3)
    plot(true_eph_eci(1,:), dd(3,:,i));
    hold on;
end
xlabel('Time [min]');
ylabel('dy [m]');
grid;
hold off;

subplot(3,1,3);
for i = 1:size(dd,3)
    plot(true_eph_eci(1,:), dd(4,:,i));
    hold on;
end
xlabel('Time [min]');
ylabel('dz [m]');
grid;
hold off;