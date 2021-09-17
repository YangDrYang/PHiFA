close all

sim_out_original = propagator.output_table;
sim_out = sort_sim_out_radau(sim_out_original);
tep = zeros(1,length(sim_out));
grav_for = zeros(3,length(sim_out));
grav_tor = zeros(3,length(sim_out));
grav_for_norm = zeros(1,length(sim_out));
grav_tor_norm = zeros(1,length(sim_out));
drag_for = zeros(3,length(sim_out));
drag_tor = zeros(3,length(sim_out));
drag_tor_norm = zeros(1,length(sim_out));
drag_for_norm = zeros(1,length(sim_out));
srad_for = zeros(3,length(sim_out));
srad_tor = zeros(3,length(sim_out));
srad_for_norm = zeros(1,length(sim_out));
srad_tor_norm = zeros(1,length(sim_out));
alb_for = zeros(3,length(sim_out));
alb_tor = zeros(3,length(sim_out));
alb_for_norm = zeros(1,length(sim_out));
alb_tor_norm = zeros(1,length(sim_out));
ifd_for = zeros(3,length(sim_out));
ifd_tor = zeros(3,length(sim_out));
ifd_for_norm = zeros(1,length(sim_out));
ifd_tor_norm = zeros(1,length(sim_out));
for i = 1:length(sim_out)
    tep(i) = sim_out(i).t;
    grav_for(:,i) = sim_out(i).grav.for;
    grav_tor(:,i) = sim_out(i).grav.tor;
    grav_for_norm(i) = norm(grav_for(:,i));
    grav_tor_norm(i) = norm(grav_tor(:,i));
    drag_for(:,i) = sim_out(i).drag.for;
    drag_tor(:,i) = sim_out(i).drag.tor;
    drag_for_norm(i) = norm(drag_for(:,i));
    drag_tor_norm(i) = norm(drag_tor(:,i));
    if drag_tor_norm(i) > 1
        drag_for(:,i) = zeros(3,1);
        drag_tor(:,i) = zeros(3,1);
        drag_tor_norm(i) = 0;
    end
    srad_for(:,i) = sim_out(i).srad.for;
    srad_tor(:,i) = sim_out(i).srad.tor;
    srad_for_norm(i) = norm(srad_for(:,i));
    srad_tor_norm(i) = norm(srad_tor(:,i));
    alb_for(:,i) = sim_out(i).albedo.for;
    alb_tor(:,i) = sim_out(i).albedo.tor;
    alb_for_norm(i) = norm(alb_for(:,i));
    alb_tor_norm(i) = norm(alb_tor(:,i));    
    ifd_for(:,i) = sim_out(i).infrared.for;
    ifd_tor(:,i) = sim_out(i).infrared.tor;
    ifd_for_norm(i) = norm(ifd_for(:,i));
    ifd_tor_norm(i) = norm(ifd_tor(:,i));    
end

figure(1)
subplot(5,1,1)
plot(tep,grav_tor_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Gravity-gradient', 'Torque Magnitude [Nm]'},'Fontsize',14);
subplot(5,1,2)
plot(tep,drag_tor_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Torque', 'Magnitude [Nm]'},'Fontsize',14);
subplot(5,1,3)
plot(tep,srad_tor_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Torque', 'Magnitude [Nm]'},'Fontsize',14);
subplot(5,1,4)
plot(tep,alb_tor_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Torque', 'Magnitude [Nm]'},'Fontsize',14);
subplot(5,1,5)
plot(tep,ifd_tor_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Torque', 'Magnitude [Nm]'},'Fontsize',14);

figure(2)
subplot(5,1,1)
plot(tep,grav_for_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Grav Force', 'Magnitude [N]'},'Fontsize',14);
subplot(5,1,2)
plot(tep,drag_for_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Force', 'Magnitude [N]'},'Fontsize',14);
subplot(5,1,3)
plot(tep,srad_for_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Force', 'Magnitude [N]'},'Fontsize',14);
subplot(5,1,4)
plot(tep,alb_for_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Force', 'Magnitude [N]'},'Fontsize',14);
subplot(5,1,5)
plot(tep,ifd_for_norm);
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Force', 'Magnitude [N]'},'Fontsize',14);

figure(3)
subplot(5,1,1)
plot(tep,grav_for(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Grav Force', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,grav_for(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Grav Force', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,grav_for(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Grav Force', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,2)
plot(tep,drag_for(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Force', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,drag_for(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Force', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,drag_for(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Force', 'Torque Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);


subplot(5,1,3)
plot(tep,srad_for(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Force', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,srad_for(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Force', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,srad_for(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Force', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,4)
plot(tep,alb_for(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Force', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,alb_for(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Force', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,alb_for(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Force', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,5)
plot(tep,ifd_for(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Force', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,ifd_for(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Force', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,ifd_for(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Force', 'Z Component [Nm]'},'Fontsize',14);
legend('x','y','z','Fontsize',14);

figure(4)
subplot(5,1,1)
plot(tep,grav_tor(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Gravity-gradient Torque', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,grav_tor(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Gravity-gradient Torque', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,grav_tor(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Gravity-gradient Torque', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,2)
plot(tep,drag_tor(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Torque', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,drag_tor(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Torque', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,drag_tor(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Drag Torque', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,3)
plot(tep,srad_tor(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Torque', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,srad_tor(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Torque', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,srad_tor(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Torque', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,4)
plot(tep,alb_tor(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Torque', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,alb_tor(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Torque', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,alb_tor(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'Albedo Torque', 'Z Component [Nm]'},'Fontsize',14);
% legend('x','y','z','Fontsize',14);

subplot(5,1,5)
plot(tep,ifd_tor(1,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Torque', 'X Component [Nm]'},'Fontsize',14);
hold on
plot(tep,ifd_tor(2,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Torque', 'Y Component [Nm]'},'Fontsize',14);
plot(tep,ifd_tor(3,:))
xlim([0 tep(end)])
xlabel('Time [sec]','Fontsize',14);
ylabel({'IR Torque', 'Z Component [Nm]'},'Fontsize',14);
legend('x','y','z','Fontsize',14);