function compareSimulations_(varargin)

if length(varargin) < 2
    fprintf('Nothing to compare, add simulations.\n');
end

prop_data = cell(1, nargin);
pert_data = cell(1, nargin);
mass = -1;
masstomult(1:nargin) = false;
skipped = 0;
for i = 1:nargin
    if isnumeric(varargin{i})
        tmpname = sprintf('logfiles//sim%04d.mat',varargin{i});
        if ~isfile(tmpname)
            fprintf('Skipping simulation. File not found\n');
            skipped = skipped + 1;
            if nargin-skipped<2
                fprintf('Nothing to compare to left, add simulation.\n');
                break;
            else
                continue;
            end
        end

        clear('propagator', 'initstate', 'eph', 'Y0', 'timestamp', 'profiler_info');
        tmpname = sprintf('logfiles//sim%04d.mat',varargin{i});
        load(tmpname, 'propagator', 'initstate');
        sim_data = struct('propagator', propagator,...
            'initstate', initstate);
        
        [pert_data{i}, prop_data{i}] = reformatSimulationData_radau(sim_data.initstate, ...
            sim_data.propagator.output_table);
        
        if mass<0
            mass = propagator.rso.mass;
        end
    else
        if mass>0
            [pert_data{i}, prop_data{i}] = getDSPOSESimulationData(varargin{i}, mass);
        else
            [pert_data{i}, prop_data{i}] = getDSPOSESimulationData(varargin{i}, 1);
            masstomult(i) = true;
        end
    end

end

if mass>0
    for i = 1:length(masstomult)
        if masstomult(i)
            pert_data{i}(2:4,:) = pert_data{i}(2:4,:).*mass;
            pert_data{i}(8:10,:) = pert_data{i}(8:10,:).*mass;
            pert_data{i}(14:16,:) = pert_data{i}(14:16,:).*mass;
        end
    end
else
    fprintf('Missing mass. Can not construct forces.\n');
    return;
end
pert_data = pert_data(~cellfun('isempty',pert_data));
prop_data = prop_data(~cellfun('isempty',prop_data));

measure_prop = prop_quat2ang(prop_data{1});
measure_pert = pert_data{1};
comp_prop = cell(1, length(prop_data)-1);
comp_pert = cell(1, length(prop_data)-1);

if size(measure_pert,1)<31
    measure_pert_new = zeros(31,size(measure_pert,2));
    measure_pert_new(1,:) = measure_pert(1,:);
    measure_pert_new(14:end,:) = measure_pert(2:end,:);
    measure_pert = measure_pert_new;
end

for i = 2:length(prop_data)
    pert = pert_data{i};
    prop = prop_data{i};
    
    if size(measure_pert,1)>size(pert,1)
        pert_new = zeros(size(measure_pert,1),size(pert,2));
        pert_new(1,:) = pert(1,:);
        pert_new(14:end,:) = pert(2:end,:);
        pert = pert_new;
    elseif size(measure_pert,1)<size(pert,1)
        pert_new = zeros(size(measure_pert,1),size(pert,2));
        pert_new(1,:) = pert(1,:);
        pert_new(2:end,:) = pert(14:end,:);
        pert = pert_new;
    end
    
    comp_prop{i-1} = zeros(size(measure_prop));
    prop = prop_quat2ang(prop);
    if ~isequal(size(measure_prop),size(prop))
        comp_prop{i-1}(1,:) = measure_prop(1,:);
        for j = 2:size(prop,1)
            comp_prop{i-1}(j,:) = interp1(prop(1,:),prop(j,:),measure_prop(1,:));
        end
    else
        comp_prop{i-1} = prop;
    end
    
    comp_pert{i-1} = zeros(size(measure_pert));
    if ~isequal(size(measure_pert),size(pert))
        comp_pert{i-1}(1,:) = measure_pert(1,:);
        for j = 2:size(pert,1)
            comp_pert{i-1}(j,:) = interp1(pert(1,:),pert(j,:),measure_pert(1,:));
        end
    else
        comp_pert{i-1} = pert;
    end
    
%     compareMeasurements(measure_prop(1,:), measure_pert(20:22,:),comp_pert{i-1}(20:22,:), ...
%         'Drag Force [N]', 'Time [sec]');
%     compareMeasurements(measure_prop(1,:), measure_pert(23:25,:),comp_pert{i-1}(23:25,:), ...
%         'Drag Torque [Nm]', 'Time [sec]');
    
    comp_prop{i-1}(2:7,:) = comp_prop{i-1}(2:7,:) - measure_prop(2:7,:);
%     comp_prop{i-1}(8:13,:) = angdiff(comp_prop{i-1}(8:13,:), measure_prop(8:13,:));
    comp_prop{i-1}(8:10,:) = angdiff(comp_prop{i-1}(8:10,:), measure_prop(8:10,:))/pi*180;
    comp_prop{i-1}(11:13,:) = (comp_prop{i-1}(11:13,:) - measure_prop(11:13,:))/pi*180;
    comp_pert{i-1} = comp_pert{i-1} - measure_pert;
    
    comp_prop{i-1}(isnan(comp_prop{i-1}))=0;
    comp_pert{i-1}(isnan(comp_pert{i-1}))=0;
end

com_pos_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_vel_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_ang_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_angvel_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_grafor_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_grator_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_dragfor_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_dragtor_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_srpfor_norm = zeros(length(comp_prop),size(comp_prop{1},2));
com_srptor_norm = zeros(length(comp_prop),size(comp_prop{1},2));
for j = 1:length(comp_prop)
    for i = 1:size(comp_prop{1},2)
        com_pos_norm(j,i) = norm(comp_prop{j}(2:4,i));
        com_vel_norm(j,i) = norm(comp_prop{j}(5:7,i));
%         com_ang_norm(j,i) = wrapToPi(norm(comp_prop{j}(8:10,i)));
        com_ang_norm(j,i) = norm(comp_prop{j}(8:10,i));
        com_angvel_norm(j,i) = norm(comp_prop{j}(11:13,i));
        com_grafor_norm(j,i) = norm(comp_pert{j}(14:16,i));
        com_grator_norm(j,i) = norm(comp_pert{j}(17:19,i));
        com_dragfor_norm(j,i) = norm(comp_pert{j}(20:22,i));
        com_dragtor_norm(j,i) = norm(comp_pert{j}(23:25,i));
        com_srpfor_norm(j,i) = norm(comp_pert{j}(26:28,i));
        com_srptor_norm(j,i) = norm(comp_pert{j}(29:31,i));        
    end
end

npnt = measure_pert(1,end);
% figure('Name', 'Pos & Vel Comparison');
% subplot(2,1,1)
% for i = 1:length(comp_prop)
%     h(1) = plot(measure_prop(1,:), comp_prop{i}(2,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_prop{i}(3,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_prop{i}(4,:),'.');
%     h(4) = plot(measure_prop(1,:), com_pos_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% ylabel('Position [m]','Fontsize',14);
% % title('Translational Motion Comparison','Fontsize',14);
% grid;
% hold off;
% subplot(2,1,2)
% for i = 1:length(comp_prop)
%     h(1) = plot(measure_prop(1,:), comp_prop{i}(5,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_prop{i}(6,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_prop{i}(7,:),'.');
%     h(4) = plot(measure_prop(1,:), com_vel_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% ylabel('Velocity [m/s]','Fontsize',14);
% legend('x','y','z','norm','Fontsize',14,'Location','northwest');
% grid;
% hold off;
% 
% figure('Name', 'Ang & Ang. Vel. Comparison');
% subplot(2,1,1)
% for i = 1:length(comp_prop)
%     h(1) = plot(measure_prop(1,:), comp_prop{i}(8,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_prop{i}(9,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_prop{i}(10,:),'.');
%     h(4) = plot(measure_prop(1,:), com_ang_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% ylabel('Angles [deg]','Fontsize',14);
% % title('Rotational Motion Comparison','Fontsize',14);
% grid;
% hold off;
% 
% subplot(2,1,2)
% for i = 1:length(comp_prop)
%     h(1) = plot(measure_prop(1,:), comp_prop{i}(11,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_prop{i}(12,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_prop{i}(13,:),'.');
%     h(4) = plot(measure_prop(1,:), com_angvel_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% % legend('x','y','z','Fontsize',14,'Location','northwest');
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% ylabel('Angular velocity [deg/s]','Fontsize',14);
% legend('yaw','pitch','roll','norm','Fontsize',14,'Location','northwest');
% 
% grid;
% hold off;
% 
% figure('Name', 'Gravity Comparison');
% % yyaxis left
% subplot(2,1,1)
% for i = 1:length(comp_pert)
%     h(1) = plot(measure_prop(1,:), comp_pert{i}(14,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_pert{i}(15,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_pert{i}(16,:),'.');
%     h(4) = plot(measure_prop(1,:), com_grafor_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['xf' num2str(i)];['yf' num2str(i)];['zf' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% ylabel('Gravitational Force [N]','Fontsize',14);
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% % title('Gravitational Force/Torque Comparison','Fontsize',14);
% grid;
% hold off;
% 
% subplot(2,1,2)
% for i = 1:length(comp_pert)
%     h(1) = plot(measure_prop(1,:), comp_pert{i}(17,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_pert{i}(18,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_pert{i}(19,:),'.');
%     h(4) = plot(measure_prop(1,:), com_grator_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['xt' num2str(i)];['yt' num2str(i)];['zt' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% ylabel('Gravitational Torque [Nm]','Fontsize',14);
% legend('x','y','z','norm','Fontsize',14);
% grid;
% hold off;
% 
% figure('Name', 'Drag Comparison');
% % yyaxis left
% subplot(2,1,1)
% for i = 1:length(comp_pert)
%     h(1) = plot(measure_prop(1,:), comp_pert{i}(20,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_pert{i}(21,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_pert{i}(22,:),'.');
%     h(4) = plot(measure_prop(1,:), com_dragfor_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['xf' num2str(i)];['yf' num2str(i)];['zf' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% ylabel('Drag Force [N]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% xlabel('Time [sec]','Fontsize',14);
% % title('Drag Force/Torque Comparison','Fontsize',14);
% grid;
% hold off;
% 
% subplot(2,1,2)
% for i = 1:length(comp_pert)
%     h(1) = plot(measure_prop(1,:), comp_pert{i}(23,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_pert{i}(24,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_pert{i}(25,:),'.');
%     h(4) = plot(measure_prop(1,:), com_dragtor_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['xt' num2str(i)];['yt' num2str(i)];['zt' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% % legend show;
% % legend('x','y','z','Fontsize',14);
% % legend('boxoff')
% xlabel('Time [sec]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% ylabel('Drag Torque [Nm]','Fontsize',14);
% legend('x','y','z','norm','Fontsize',14);
% grid;
% hold off;
% 
% figure('Name', 'SRP Comparison');
% % yyaxis left
% subplot(2,1,1)
% for i = 1:length(comp_pert)
%     h(1) = plot(measure_prop(1,:), comp_pert{i}(26,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_pert{i}(27,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_pert{i}(28,:),'.');
%     h(4) = plot(measure_prop(1,:), com_srpfor_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['xf' num2str(i)];['yf' num2str(i)];['zf' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% ylabel('SRP Force [N]','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% xlabel('Time [sec]','Fontsize',14);
% % title('SRP Force/Torque Comparison','Fontsize',14);
% grid;
% hold off;
% 
% subplot(2,1,2)
% for i = 1:length(comp_pert)
%     h(1) = plot(measure_prop(1,:), comp_pert{i}(29,:),'.');
%     hold on;
%     h(2) = plot(measure_prop(1,:), comp_pert{i}(30,:),'.');
%     h(3) = plot(measure_prop(1,:), comp_pert{i}(31,:),'.');
%     h(4) = plot(measure_prop(1,:), com_srptor_norm(i,:),'.');
% %     set(h,{'DisplayName'},{['xt' num2str(i)];['yt' num2str(i)];['zt' num2str(i)]})
% %     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
% end
% legend('x','y','z','norm','Fontsize',14);
% xlim([0 npnt])
% % xlim([0 20])
% xlabel('Time [sec]','Fontsize',14);
% ylabel('SRP Torque [Nm]','Fontsize',14);
% grid;
% hold off;

figure('Name', 'State Comparison');
subplot(2,2,1)
for i = 1:length(comp_prop)
    h(1) = plot(measure_prop(1,:), comp_prop{i}(2,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_prop{i}(3,:),'.');
    h(3) = plot(measure_prop(1,:), comp_prop{i}(4,:),'.');
    h(4) = plot(measure_prop(1,:), com_pos_norm(i,:),'.');
%     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
ylabel('Position [m]','Fontsize',14);
% title('Translational Motion Comparison','Fontsize',14);
grid;
hold off;
subplot(2,2,2)
for i = 1:length(comp_prop)
    h(1) = plot(measure_prop(1,:), comp_prop{i}(5,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_prop{i}(6,:),'.');
    h(3) = plot(measure_prop(1,:), comp_prop{i}(7,:),'.');
    h(4) = plot(measure_prop(1,:), com_vel_norm(i,:),'.');
%     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
ylabel('Velocity [m/s]','Fontsize',14);
legend('x','y','z','norm','Fontsize',14,'Location','northwest');
grid;
hold off;

subplot(2,2,3)
for i = 1:length(comp_prop)
    h(1) = plot(measure_prop(1,:), comp_prop{i}(8,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_prop{i}(9,:),'.');
    h(3) = plot(measure_prop(1,:), comp_prop{i}(10,:),'.');
    h(4) = plot(measure_prop(1,:), com_ang_norm(i,:),'.');
%     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
ylabel('Angles [deg]','Fontsize',14);
% title('Rotational Motion Comparison','Fontsize',14);
grid;
hold off;

subplot(2,2,4)
for i = 1:length(comp_prop)
    h(1) = plot(measure_prop(1,:), comp_prop{i}(11,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_prop{i}(12,:),'.');
    h(3) = plot(measure_prop(1,:), comp_prop{i}(13,:),'.');
    h(4) = plot(measure_prop(1,:), com_angvel_norm(i,:),'.');
%     set(h,{'DisplayName'},{['x' num2str(i)];['y' num2str(i)];['z' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
% legend('x','y','z','Fontsize',14,'Location','northwest');
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
ylabel({'Angular', 'velocity [deg/s]'},'Fontsize',14);
legend('yaw','pitch','roll','norm','Fontsize',14,'Location','northwest');

grid;
hold off;

figure('Name', 'Forces Comparison');
% yyaxis left
subplot(3,2,1)
for i = 1:length(comp_pert)
    h(1) = plot(measure_prop(1,:), comp_pert{i}(14,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_pert{i}(15,:),'.');
    h(3) = plot(measure_prop(1,:), comp_pert{i}(16,:),'.');
    h(4) = plot(measure_prop(1,:), com_grafor_norm(i,:),'.');
%     set(h,{'DisplayName'},{['xf' num2str(i)];['yf' num2str(i)];['zf' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
ylabel({'Gravitational', 'Force [N]'},'Fontsize',14);
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
% title('Gravitational Force/Torque Comparison','Fontsize',14);
grid;
hold off;

subplot(3,2,2)
for i = 1:length(comp_pert)
    h(1) = plot(measure_prop(1,:), comp_pert{i}(17,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_pert{i}(18,:),'.');
    h(3) = plot(measure_prop(1,:), comp_pert{i}(19,:),'.');
    h(4) = plot(measure_prop(1,:), com_grator_norm(i,:),'.');
%     set(h,{'DisplayName'},{['xt' num2str(i)];['yt' num2str(i)];['zt' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
ylabel({'Gravitational', 'Torque [Nm]'},'Fontsize',14);
legend('x','y','z','norm','Fontsize',14);
grid;
hold off;

% yyaxis left
subplot(3,2,3)
for i = 1:length(comp_pert)
    h(1) = plot(measure_prop(1,:), comp_pert{i}(20,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_pert{i}(21,:),'.');
    h(3) = plot(measure_prop(1,:), comp_pert{i}(22,:),'.');
    h(4) = plot(measure_prop(1,:), com_dragfor_norm(i,:),'.');
%     set(h,{'DisplayName'},{['xf' num2str(i)];['yf' num2str(i)];['zf' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
ylabel({'Drag Force','[N]'},'Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
xlabel('Time [sec]','Fontsize',14);
% title('Drag Force/Torque Comparison','Fontsize',14);
grid;
hold off;

subplot(3,2,4)
for i = 1:length(comp_pert)
    h(1) = plot(measure_prop(1,:), comp_pert{i}(23,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_pert{i}(24,:),'.');
    h(3) = plot(measure_prop(1,:), comp_pert{i}(25,:),'.');
    h(4) = plot(measure_prop(1,:), com_dragtor_norm(i,:),'.');
%     set(h,{'DisplayName'},{['xt' num2str(i)];['yt' num2str(i)];['zt' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
% legend show;
% legend('x','y','z','Fontsize',14);
% legend('boxoff')
xlabel('Time [sec]','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
ylabel({'Drag Torque', '[Nm]'},'Fontsize',14);
legend('x','y','z','norm','Fontsize',14);
grid;
hold off;

% yyaxis left
subplot(3,2,5)
for i = 1:length(comp_pert)
    h(1) = plot(measure_prop(1,:), comp_pert{i}(26,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_pert{i}(27,:),'.');
    h(3) = plot(measure_prop(1,:), comp_pert{i}(28,:),'.');
    h(4) = plot(measure_prop(1,:), com_srpfor_norm(i,:),'.');
%     set(h,{'DisplayName'},{['xf' num2str(i)];['yf' num2str(i)];['zf' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
ylabel({'SRP Force', '[N]'},'Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
xlabel('Time [sec]','Fontsize',14);
% title('SRP Force/Torque Comparison','Fontsize',14);
grid;
hold off;

subplot(3,2,6)
for i = 1:length(comp_pert)
    h(1) = plot(measure_prop(1,:), comp_pert{i}(29,:),'.');
    hold on;
    h(2) = plot(measure_prop(1,:), comp_pert{i}(30,:),'.');
    h(3) = plot(measure_prop(1,:), comp_pert{i}(31,:),'.');
    h(4) = plot(measure_prop(1,:), com_srptor_norm(i,:),'.');
%     set(h,{'DisplayName'},{['xt' num2str(i)];['yt' num2str(i)];['zt' num2str(i)]})
%     set(h,{'DisplayName'},{['x' ];['y' ];['z' ]})
end
legend('x','y','z','norm','Fontsize',14);
xlim([0 npnt])
% xlim([0 20])
xlabel('Time [sec]','Fontsize',14);
ylabel({'SRP Torque', '[Nm]'},'Fontsize',14);
grid;
hold off;
end