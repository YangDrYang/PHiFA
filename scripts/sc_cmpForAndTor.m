% sc_cmpForAndTor.m
% compare forces and torques plus impulse and ang momentum change
%% load data
% sim_data = loadSimDataFolder('poinacc_and_power_plus_logging');
sim_data = loadSimDataFolder('');
[pert_data, prop_data] = reformatSimulationData(sim_data{3}.initstate, ...
            sim_data{3}.propagator.output_table);

GM = 398600.4418e9; 
mass = sim_data{3}.propagator.rso.mass;
%% make plotable data
n = length(pert_data);
time = zeros(n,1);
impuls = zeros(n,1);
angmom = zeros(n,1);
force = zeros(n,1);
torque = zeros(n,1);
grav = zeros(n,1);
gravt = zeros(n,1);
srp = zeros(n,1);
srpt = zeros(n,1);
drag = zeros(n,1);
dragt = zeros(n,1);

for i = 1:n
    time(i) = pert_data(1,i);
    impuls(i) = log10(norm(pert_data(2:4,i)));
    angmom(i) = log10(norm(pert_data(5:7,i)));
    force(i) = log10(norm(pert_data(8:10,i)));
    torque(i) = log10(norm(pert_data(11:13,i)));
    grav(i) = log10(norm(pert_data(14:16,i)+...
        GM*prop_data(2:4,i)/norm(prop_data(2:4,i))^3.*mass));
    gravt(i) = log10(norm(pert_data(17:19,i)));
    drag(i) = log10(norm(pert_data(20:22,i)));
    dragt(i) = log10(norm(pert_data(23:25,i)));
    srp(i) = log10(norm(pert_data(26:28,i)));
    srpt(i) = log10(norm(pert_data(29:31,i)));

end
time = time-200;
%% plot the sh*t

figure;
subplot(2,1,1);
plot(time(:), impuls(:));
ylabel('Impulse [Ns]');
legend('Pulsed laser','Location', 'northwest');
grid;
xlim([0 400]);
subplot(2,1,2);
plot(time(:), force(:));
hold on;
% plot(time(:), grav(:));
plot(time(:), drag(:));
plot(time(:), srp(:));
hold off;
ylabel('Force [N]');
xlabel('Seconds [s]');
grid;
xlim([0 400]);
legend('CW Laser', 'Drag', 'SRP', 'Location', 'northwest');

plot2tikz('cmpForces', 0.6);

%% plus torques
figure;
subplot(2,1,1);
plot(time(:), angmom(:));
ylabel('Angular Momentum [Nms]');
legend('Pulsed laser','Location', 'northwest');
grid;
xlim([0 400]);
subplot(2,1,2);
plot(time(:), torque(:));
hold on;
% plot(time(:), grav(:));
plot(time(:), dragt(:));
plot(time(:), srpt(:));
hold off;
ylabel('Torque [Nm]');
xlabel('Seconds [s]');
grid;
xlim([0 400]);
legend('CW Laser', 'Drag', 'SRP', 'Location', 'northwest');

plot2tikz('cmpTorques', 0.6);