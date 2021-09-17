function [prop, pert] = plotSimulationOutput_PHiFA(is, sim_out, figname)

% close all

time(1,length(sim_out)-1) = 0;
laser_imp(3,length(sim_out)-1) = 0;
laser_for(3,length(sim_out)-1) = 0;
laser_projarea(1,length(sim_out)-1) = 0;
grav(3,length(sim_out)-1) = 0;
sun_for(3,length(sim_out)-1) = 0;
moon_for(3,length(sim_out)-1) = 0;
srad(3,length(sim_out)-1) = 0;
srad_projarea(1,length(sim_out)-1) = 0;
alb(3,length(sim_out)-1) = 0;
alb_projarea(1,length(sim_out)-1) = 0;
ir(3,length(sim_out)-1) = 0;
ir_projarea(1,length(sim_out)-1) = 0;
drag(3,length(sim_out)-1) = 0;
drag_projarea(1,length(sim_out)-1) = 0;
laser_mom(3,length(sim_out)-1) = 0;
laser_tor(3,length(sim_out)-1) = 0;
gravt(3,length(sim_out)-1) = 0;
sradt(3,length(sim_out)-1) = 0;
albt(3,length(sim_out)-1) = 0;
irt(3,length(sim_out)-1) = 0;
dragt(3,length(sim_out)-1) = 0;
magt(3,length(sim_out)-1) = 0;
energy(5,length(sim_out)-1) = 0;
eph(14,length(sim_out)-1) = 0;

pert = zeros(31, length(sim_out)-1);

for i = 2:length(sim_out)
    if i~=1 && sim_out(i).t==0
        time = time(:,1:i-2);
        laser_imp = laser_imp(:,1:i-2);
        laser_for = laser_for(:,1:i-2);
        laser_projarea = laser_projarea(:,1:i-2);
        grav = grav(:,1:i-2);
        sun_for = sun_for(:,1:i-2);
        moon_for = moon_for(:,1:i-2);
        srad = srad(:,1:i-2);
        alb = alb(:,1:i-2);
        ir = ir(:,1:i-2);
        drag = drag(:,1:i-2);
        laser_mom = laser_mom(:,1:i-2);
        laser_tor = laser_tor(:,1:i-2);
        gravt = gravt(:,1:i-2);
        sradt = sradt(:,1:i-2);
        albt = albt(:,1:i-2);
        irt = irt(:,1:i-2);
        dragt = dragt(:,1:i-2);
        magt = magt(:,1:i-2);
        energy = energy(:,1:i-2);
        eph = eph(:,1:i-2);
        pert = pert(:,1:i-2);
        break;
    end
    time(i-1) = sim_out(i).t(end);
    eph(1,i-1) = time(i-1);
    eph(2:14,i-1) = sim_out(i).Y;
    if is.a_lase
        laser_imp(1:3,i-1) = sim_out(i).laser.imp;
        laser_imp(4,i-1) = norm(sim_out(i).laser.imp);
        laser_for(1:3,i-1) = sim_out(i).laser.for;
        laser_for(4,i-1) = norm(sim_out(i).laser.for);
        laser_projarea(1,i-1) = sim_out(i).laser.projarea;
    else
        laser_imp(1:4,i-1) = [0;0;0;0];
        laser_for(1:4,i-1) = [0;0;0;0];
        laser_projarea(1,i-1) = 0;
    end
    if is.a_drag
        drag(1:3,i-1) = sim_out(i).drag.for;
        drag(4,i-1) = norm(sim_out(i).drag.for);
        drag_projarea(1,i-1) = sim_out(i).drag.projarea;
    else
        drag(1:4,i-1) = [0;0;0;0];
        drag_projarea(1,i-1) = 0;
    end
    if is.a_srad
        srad(1:3,i-1) = sim_out(i).srad.for;
        srad(4,i-1) = norm(sim_out(i).srad.for);
        srad_projarea(1,i-1) = sim_out(i).srad.projarea;
    else
        srad(1:4,i-1) = [0;0;0;0];
        srad_projarea(1,i-1) = 0;
    end
    if is.a_erad
        alb(1:3,i-1) = sim_out(i).albedo.for;
        alb(4,i-1) = norm(sim_out(i).albedo.for);
        alb_projarea(1,i-1) = sim_out(i).albedo.projarea;
        ir(1:3,i-1) = sim_out(i).infrared.for;
        ir(4,i-1) = norm(sim_out(i).infrared.for);
        ir_projarea(1,i-1) = sim_out(i).infrared.projarea;        
    else
        ir(1:4,i-1) = [0;0;0;0];
        ir_projarea(1,i-1) = 0;
    end    
    if is.a_grav
        grav(1:3,i-1) = sim_out(i).grav.for;
        grav(4,i-1) = norm(sim_out(i).grav.for);
    else
        grav(1:4,i-1) = [0;0;0;0];
    end
    if is.a_sun
        sun_for(1:3,i-1) = sim_out(i).sun.for;
        sun_for(4,i-1) = norm(sim_out(i).sun.for);
    else
        sun_for(1:4,i-1) = [0;0;0;0];
    end    
    if is.a_moon
        moon_for(1:3,i-1) = sim_out(i).moon.for;
        moon_for(4,i-1) = norm(sim_out(i).moon.for);
    else
        moon_for(1:4,i-1) = [0;0;0;0];
    end    
    
    if is.g_lase
        laser_mom(1:3,i-1) = sim_out(i).laser.mom;
        laser_mom(4,i-1) = norm(sim_out(i).laser.mom);
        laser_tor(1:3,i-1) = sim_out(i).laser.tor;
        laser_tor(4,i-1) = norm(sim_out(i).laser.tor);
    else
        laser_mom(1:4,i-1) = [0;0;0;0];
        laser_tor(1:4,i-1) = [0;0;0;0];
    end
    if is.g_drag
        dragt(1:3,i-1) = sim_out(i).drag.tor;
        dragt(4,i-1) = norm(sim_out(i).drag.tor);
    else
        dragt(1:4,i-1) = [0;0;0;0];
    end
    if is.g_srad
        sradt(1:3,i-1) = sim_out(i).srad.tor;
        sradt(4,i-1) = norm(sim_out(i).srad.tor);
    else
        sradt(1:4,i-1) = [0;0;0;0];
    end
    if is.g_erad
        albt(1:3,i-1) = sim_out(i).albedo.tor;
        albt(4,i-1) = norm(sim_out(i).albedo.tor);
        irt(1:3,i-1) = sim_out(i).infrared.tor;
        irt(4,i-1) = norm(sim_out(i).infrared.tor);        
    else
        albt(1:4,i-1) = [0;0;0;0];        
        irt(1:4,i-1) = [0;0;0;0];
    end
    if is.g_grav
        gravt(1:3,i-1) = sim_out(i).grav.tor;
        gravt(4,i-1) = norm(sim_out(i).grav.tor);
    else
        gravt(1:4,i-1) = [0;0;0;0];
    end
    if is.g_mag
        magt(1:3,i-1) = sim_out(i).mag.tor;
        magt(4,i-1) = norm(sim_out(i).mag.tor);
    else
        magt(1:4,i-1) = [0;0;0;0];
    end    
    energy(1:5,i-1) = sim_out(i).energy;
end
pert(1, :) = time;
% time = time./60;

nrows = 2;
pert(nrows:nrows+2, :) = laser_imp(1:3,:);
pert(nrows+3:nrows+5, :) = laser_mom(1:3,:);
nrows = nrows+6;
pert(nrows:nrows+2, :) = laser_for(1:3,:);
pert(nrows+3:nrows+5, :) = laser_tor(1:3,:);
nrows = nrows+6;
pert(nrows:nrows+2, :) = grav(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = gravt(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = drag(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = dragt(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = srad(1:3,:);
nrows = nrows+3;
pert(nrows:nrows+2, :) = sradt(1:3,:);

% prop = zeros(size(eph,1)+1,length(eph));
prop = eph;

% plotEphemerides(eph);
plotEphemerides_new(eph);
npnt = eph(1,end);

% if is.a_srad && is.a_drag && is.a_lase && is.a_erad
if is.a_srad && is.a_drag && is.a_erad
    figure('Name',figname);
    
    plot(time, drag_projarea(1:length(time)));
    hold on
    plot(time, srad_projarea(1:length(time)));
    plot(time, alb_projarea(1:length(time)));
    plot(time, ir_projarea(1:length(time)));
    xlabel('Time [sec]','Fontsize',14);
    ylabel('Projected Area [m^2]','Fontsize',14);
%     title('Projected Areas for Different Perturbations','Fontsize',14);
%     legend('laser','drag','srp','albedo','infrared','Fontsize',14);
    legend('drag','srp','albedo','infrared');
%     ylim([-0.05 0.20]);
    xlim([0 npnt])
%     xlim([0 20])
    grid;
    hold off;
end

% %% albedo and infrared
% if is.a_erad
%     figure('Name', figname);
%     subplot(2,1,1)
%     plot(time,alb(1,:));
%     hold on 
%     plot(time,alb(2,:));
%     plot(time,alb(3,:));
%     plot(time,alb(4,:));
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel('Albedo Radiation Pressure [N]','Fontsize',14);
%     
%     subplot(2,1,2)
%     plot(time,ir(1,:));
%     hold on 
%     plot(time,ir(2,:));
%     plot(time,ir(3,:));
%     plot(time,ir(4,:));
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel('Infrared Radiation Pressure [N]','Fontsize',14);
% end
% 
% if is.g_erad
%     figure('Name', figname);
%     subplot(2,1,1)
%     plot(time,albt(1,:));
%     hold on 
%     plot(time,albt(2,:));
%     plot(time,albt(3,:));
%     plot(time,albt(4,:));
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel('Albedo Radiation Torque [Nm]','Fontsize',14);
%     
%     subplot(2,1,2)
%     plot(time,irt(1,:));
%     hold on 
%     plot(time,irt(2,:));
%     plot(time,irt(3,:));
%     plot(time,irt(4,:));
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel('Infrared Radiation Torque [Nm]','Fontsize',14);
% end

%% grav srp albedo infrared drag in subplot force
% if is.a_grav && is.a_srad && is.a_drag && is.a_erad && is.a_lase
% if is.a_grav && is.a_srad && is.a_drag && is.a_erad
%     figure('Name', figname);
%     h(1) = subplot(7,1,1);
%     plot(time, grav(1,:));
%     hold on;
%     plot(time, grav(2,:));
%     plot(time, grav(3,:));
% %     plot(time, grav(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Gravity';'Force [N]'},'Fontsize',14);
% %     ylabel('Gravitational Force [N]','Fontsize',14);
% %     axis([0 100 0 35])
% %     legend('x','y','z');
% %     title('Gravitational Force','Fontsize',14);
%     grid;
%     hold off;
% 
%     h(2) = subplot(7,1,2);
%     plot(time, sun_for(1,:));
%     hold on;
%     plot(time, sun_for(2,:));
%     plot(time, sun_for(3,:));
% %     plot(time, srad(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Sun'; 'Attraction [N]'},'Fontsize',14);
% %     legend('x','y','z','Fontsize',14);
%     grid;
%     hold off;   
%     
%     h(3) = subplot(7,1,3);
%     plot(time, moon_for(1,:));
%     hold on;
%     plot(time, moon_for(2,:));
%     plot(time, moon_for(3,:));
% %     plot(time, srad(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Moon'; 'Attraction [N]'},'Fontsize',14);
% %     legend('x','y','z','Fontsize',14);
%     grid;
%     hold off; 
%     
%     h(4) = subplot(7,1,4);
%     plot(time, srad(1,:));
%     hold on;
%     plot(time, srad(2,:));
%     plot(time, srad(3,:));
% %     plot(time, srad(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'SRP'; 'Force[N]'},'Fontsize',14);
% %     legend('x','y','z','Fontsize',14);
% %     title('SRP Force','Fontsize',14);
%     grid;
%     hold off;   
%     
%     h(5) = subplot(7,1,5);
%     plot(time, alb(1,:));
%     hold on;
%     plot(time, alb(2,:));
%     plot(time, alb(3,:));
% %     plot(time, alb(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Albedo'; 'Force[N]'},'Fontsize',14);
% 
% %     ylabel('Albedo Radiation Pressure [N]','Fontsize',14);
% %     legend('x','y','z');
% %     title('Drag Force','Fontsize',14);
%     grid;
%     hold off;
%     
%     h(6) = subplot(7,1,6);
%     plot(time, ir(1,:));
%     hold on;
%     plot(time, ir(2,:));
%     plot(time, ir(3,:));
% %     plot(time, ir(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Infrared'; 'Force[N]'},'Fontsize',14);    
% %     ylabel('Infrared Radiation Pressure [N]','Fontsize',14);
% %     legend('x','y','z');
% %     title('Drag Force','Fontsize',14);
%     grid;
%     hold off;
%     
%     h(7) = subplot(7,1,7);
%     plot(time, drag(1,:));
%     hold on;
%     plot(time, drag(2,:));
%     plot(time, drag(3,:));
% %     plot(time, drag(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Drag'; 'Force[N]'},'Fontsize',14);
% %     legend('x','y','z');
% %     title('Drag Force','Fontsize',14);
%     grid;
%     hold off;
% end

% if is.a_grav && is.a_srad && is.a_drag && is.a_erad
    figure('Name', figname);
    h(1) = subplot(3,1,1);
    plot(time, grav(1,:));
    hold on;
    plot(time, grav(2,:));
    plot(time, grav(3,:));
%     plot(time, grav(4,:));
    xlim([0 npnt])
%     xlim([0 20])
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Gravity';'Force [N]'},'Fontsize',14);
%     ylabel('Gravitational Force [N]','Fontsize',14);
%     axis([0 100 0 35])
%     legend('x','y','z');
%     title('Gravitational Force','Fontsize',14);
    grid;
    hold off;
    
    h(4) = subplot(3,1,2);
    plot(time, srad(1,:));
    hold on;
    plot(time, srad(2,:));
    plot(time, srad(3,:));
%     plot(time, srad(4,:));
    xlim([0 npnt])
%     xlim([0 20])
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'SRP'; 'Force[N]'},'Fontsize',14);
%     legend('x','y','z','Fontsize',14);
%     title('SRP Force','Fontsize',14);
    grid;
    hold off;   
    
    h(7) = subplot(3,1,3);
    plot(time, drag(1,:));
    hold on;
    plot(time, drag(2,:));
    plot(time, drag(3,:));
%     plot(time, drag(4,:));
    xlim([0 npnt])
%     xlim([0 20])
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Drag'; 'Force[N]'},'Fontsize',14);
%     legend('x','y','z');
%     title('Drag Force','Fontsize',14);
    grid;
    hold off;
% end
% %% grav srp drag in subplot torque
% % if is.g_grav && is.g_srad && is.g_drag && is.g_erad && is.g_lase
% if is.g_grav && is.g_srad && is.g_drag && is.g_erad
%     figure('Name', figname);
%     h(1) = subplot(6,1,1);
%     plot(time, gravt(1,:));
%     hold on;
%     plot(time, gravt(2,:));
%     plot(time, gravt(3,:));
% %     plot(time, gravt(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel({'Gravity';'Gradient [Nm]'},'Fontsize',14);
% %     legend('x','y','z');
% %     title('Gravitational Torque','Fontsize',14);
%     grid;
%     hold off;
% 
%     h(2) = subplot(6,1,2);
%     plot(time, dragt(1,:));
%     hold on;
%     plot(time, dragt(2,:));
%     plot(time, dragt(3,:));
% %     plot(time, dragt(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
% %     ylabel('Drag [Nm]','Fontsize',14);
% %     legend('x','y','z');
%     ylabel({'Drag';'Torque [Nm]'},'Fontsize',14);
%     grid;
%     hold off;
%     
%     h(3) = subplot(6,1,3);
%     plot(time, sradt(1,:));
%     hold on;
%     plot(time, sradt(2,:));
%     plot(time, sradt(3,:));
% %     plot(time, sradt(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
% %     ylabel('SRP [Nm]','Fontsize',14);
% %     legend('x','y','z');
%     ylabel({'SRP';'Torque [Nm]'},'Fontsize',14);
%     grid;
%     hold off;
%     
%     h(4) = subplot(6,1,4);
%     plot(time, albt(1,:));
%     hold on;
%     plot(time, albt(2,:));
%     plot(time, albt(3,:));
% %     plot(time, albt(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
% %     ylabel('SRP [Nm]','Fontsize',14);
% %     legend('x','y','z');
%     ylabel({'Albedo';'Torque [Nm]'},'Fontsize',14);
%     grid;
%     hold off; 
%     
%     h(5) = subplot(6,1,5);
%     plot(time, irt(1,:));
%     hold on;
%     plot(time, irt(2,:));
%     plot(time, irt(3,:));
% %     plot(time, irt(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
% %     ylabel('Drag [Nm]','Fontsize',14);
% %     legend('x','y','z');
%     ylabel({'Infrared';'Torque [Nm]'},'Fontsize',14);
%     grid;
%     hold off;
%     
%     h(6) = subplot(6,1,6);
%     plot(time, magt(1,:));
%     hold on;
%     plot(time, magt(2,:));
%     plot(time, magt(3,:));
% %     plot(time, magt(4,:));
%     xlim([0 npnt])
% %     xlim([0 20])
%     xlabel('Time [sec]','Fontsize',14);
% %     ylabel('Drag [Nm]','Fontsize',14);
% %     legend('x','y','z');
%     ylabel({'Magnetic';'Torque [Nm]'},'Fontsize',14);
%     grid;
%     hold off;
%     
%     
% %     pos = get(h,'Position');
% %     new = mean(cellfun(@(v)v(1),pos(1:2)));
% %     set(h(7),'Position',[new,pos{end}(2:end)])
% end

% if is.g_grav && is.g_srad && is.g_drag && is.g_erad
    figure('Name', figname);
    h(1) = subplot(3,1,1);
    plot(time, gravt(1,:));
    hold on;
    plot(time, gravt(2,:));
    plot(time, gravt(3,:));
%     plot(time, gravt(4,:));
    xlim([0 npnt])
%     xlim([0 20])
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Gravity';'Gradient [Nm]'},'Fontsize',14);
%     legend('x','y','z');
%     title('Gravitational Torque','Fontsize',14);
    grid;
    hold off;

    h(2) = subplot(3,1,2);
    plot(time, sradt(1,:));
    hold on;
    plot(time, sradt(2,:));
    plot(time, sradt(3,:));
%     plot(time, sradt(4,:));
    xlim([0 npnt])
%     xlim([0 20])
    xlabel('Time [sec]','Fontsize',14);
%     ylabel('SRP [Nm]','Fontsize',14);
%     legend('x','y','z');
    ylabel({'SRP';'Torque [Nm]'},'Fontsize',14);
    grid;
    hold off;
    
    h(3) = subplot(3,1,3);
    plot(time, dragt(1,:));
    hold on;
    plot(time, dragt(2,:));
    plot(time, dragt(3,:));
%     plot(time, dragt(4,:));
    xlim([0 npnt])
%     xlim([0 20])
    xlabel('Time [sec]','Fontsize',14);
%     ylabel('Drag [Nm]','Fontsize',14);
%     legend('x','y','z');
    ylabel({'Drag';'Torque [Nm]'},'Fontsize',14);
    grid;
    hold off;
    

% end

% %% magnetic torque
% if is.g_mag
%     figure
%     plot(time, magt(1,:));
%     hold on;
%     plot(time, magt(2,:));
%     plot(time, magt(3,:));
%     plot(time, magt(4,:));
%     xlabel('Time [sec]','Fontsize',14);
%     ylabel('Magnetic Torque [Nm]','Fontsize',14);
%     legend('x','y','z','norm','Fontsize',14);
%     grid;
%     hold off;
% end

%work-energy balance plot
% figure
% subplot(2,2,[1,2])
% plot(time, energy(5,:));
% xlabel('Time [sec]','Fontsize',14);
% subplot(2,2,3)
% plot(time, energy(1,:));
% hold on
% plot(time, energy(2,:));
% xlabel('Time [sec]','Fontsize',14);
% subplot(2,2,4)
% plot(time, energy(3,:));
% hold on
% plot(time, energy(4,:));
% xlabel('Time [sec]','Fontsize',14);

figure('Name', 'Work-energy Balance')
subplot(3,2,[1,2])
semilogx(time, energy(5,:));
xlabel('Time [sec]','Fontsize',14);
ylabel({'Work-energy' 'Balance [J]'},'Fontsize',14)
subplot(3,2,3)
semilogx(time, energy(1,:));
xlabel('Time [sec]','Fontsize',14);
ylabel({'Total Mechanical' 'Energy [J]'},'Fontsize',14)
subplot(3,2,4)
semilogx(time, energy(2,:));
xlabel('Time [sec]','Fontsize',14);
ylabel({'Translational' 'Energy [J]'},'Fontsize',14)
subplot(3,2,5)
plot(time, energy(3,:));
xlabel('Time [sec]','Fontsize',14);
ylabel({'Rotational' 'Energy [J]'},'Fontsize',14)
subplot(3,2,6)
plot(time, energy(4,:));
xlabel('Time [sec]','Fontsize',14);
ylabel({'Rotational Work-energy' 'Balance [J]'},'Fontsize',14)
end