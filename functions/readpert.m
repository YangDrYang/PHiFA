function torfor_dspose = readpert(file)
filedata = load(file);
pert_dspose = filedata.pert_dspose;
for_aero = pert_dspose(:,2:4);
tor_aero = pert_dspose(:,5:7);
for_grav = pert_dspose(:,8:10);
tor_grav = pert_dspose(:,11:13);
tor_eddy = pert_dspose(:,14:16);
for_sun = pert_dspose(:,17:18);
for_moon = pert_dspose(:,19:21);
for_srp = pert_dspose(:,22:24);
for_alb = pert_dspose(:,25:27);
for_ir = pert_dspose(:,28:30);
tor_srp = pert_dspose(:,31:33);
tor_alb = pert_dspose(:,34:36);
tor_ir = pert_dspose(:,37:39);

figure
subplot(5,1,1)
plot(pert_dspose(:,1),for_grav);
ylabel('grav force')
subplot(5,1,2)
plot(pert_dspose(:,1),for_aero);
ylabel('drag force')
subplot(5,1,3)
plot(pert_dspose(:,1),for_srp);
ylabel('srp force')
subplot(5,1,4)
plot(pert_dspose(:,1),for_alb);
ylabel('albedo force')
subplot(5,1,5)
plot(pert_dspose(:,1),for_ir);
ylabel('infrared force')

figure
subplot(6,1,1)
plot(pert_dspose(:,1),tor_grav);
ylabel('grav torque')
subplot(6,1,2)
plot(pert_dspose(:,1),tor_aero);
ylabel('drag torque')
subplot(6,1,3)
plot(pert_dspose(:,1),tor_srp);
ylabel('srp torque')
subplot(6,1,4)
plot(pert_dspose(:,1),tor_alb);
ylabel('albedo torque')
subplot(6,1,5)
plot(pert_dspose(:,1),tor_ir);
ylabel('infrared torque')
subplot(6,1,6)
plot(pert_dspose(:,1),tor_eddy);
ylabel('eddy current torque')




