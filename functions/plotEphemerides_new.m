function plotEphemerides_new(varargin)

for i = 1:length(varargin)
    time = varargin{i}(1,:);
    state = varargin{i}(2:4,:)./1000;
    velo = varargin{i}(5:7,:)./1000;

    timeperiod = varargin{i}(1,end) - varargin{i}(1,1) + 1;

    ephsize = size(varargin{i},1);
    
    angles = zeros(3, size(varargin{i},2));
    arates = zeros(3, size(varargin{i},2));
    if ephsize>7
        for j = 1:size(varargin{i},2)
            angles(:,j) = eulerd(quaternion(varargin{i}(8:11,j)'),'ZYX','frame');%radian -> degree
    %         angles(:,j) = (Q2Eul(norm_quat(varargin{i}(8:11,j)))).*180/pi;
            arates(:,j) = varargin{i}(12:14,j).*180/pi;   
        end
    end


    figure(1);
    subplot(3,2,1);
    plot(time, state(1,:));
    grid;
    xlim([0, timeperiod]);
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Position X', 'Component [km]'},'Fontsize',14);
    subplot(3,2,3);
    plot(time, state(2,:));
    grid;
    xlim([0, timeperiod]);
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Position Y', 'Component [km]'},'Fontsize',14);
    subplot(3,2,5);
    plot(time, state(3,:));
    grid;
    xlim([0, timeperiod]);
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Position Z', 'Component [km]'},'Fontsize',14); 
    
    subplot(3,2,2);
    plot(time, velo(1,:));
    grid;
    xlim([0, timeperiod]);
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Velocity X', 'Component [km/s]'},'Fontsize',14);
    subplot(3,2,4);
    plot(time, velo(2,:));
    grid;
    xlim([0, timeperiod]);
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Velocity Y', 'Component [km/s]'},'Fontsize',14);
    subplot(3,2,6);
    plot(time, velo(3,:));
    grid;
    xlim([0, timeperiod]);
    xlabel('Time [sec]','Fontsize',14);
    ylabel({'Velocity Z', 'Component [km/s]'},'Fontsize',14);    
%     axis([0 100 -6.5 6])
    if ephsize(1)>7
        figure(2)
        subplot(3,2,1);
        plot(time, angles(1,:));
        grid;
        xlim([0, timeperiod]);
        xlabel('Time [sec]','Fontsize',14);
        ylabel('Yaw Angle [deg]','Fontsize',14);
        subplot(3,2,3);
        plot(time, angles(2,:));
        grid;
        xlim([0, timeperiod]);
        xlabel('Time [sec]','Fontsize',14);
        ylabel('Pitch Angle [deg]','Fontsize',14);
        subplot(3,2,5);
        plot(time, angles(3,:));
        grid;
        xlim([0, timeperiod]);
        xlabel('Time [sec]','Fontsize',14);
        ylabel('Roll Angle [deg]','Fontsize',14);
        
        subplot(3,2,2);
        plot(time, arates(1,:));
%         plot(time(1:10:end), arates(1,1:10:end),time(1:10:end), arates(2,1:10:end),time(1:10:end), arates(3,1:10:end));
        grid;
        xlim([0, timeperiod]);
        xlabel('Time [sec]','Fontsize',14);
        ylabel({'Yaw Angular', 'Velocity [deg/s]'},'Fontsize',14);
        subplot(3,2,4);
        plot(time, arates(2,:));
%         plot(time(1:10:end), arates(1,1:10:end),time(1:10:end), arates(2,1:10:end),time(1:10:end), arates(3,1:10:end));
        grid;
        xlim([0, timeperiod]);
        xlabel('Time [sec]','Fontsize',14);
        ylabel({'Pitch Angular', 'Velocity [deg/s]'},'Fontsize',14);
        subplot(3,2,6);
        plot(time, arates(3,:));
%         plot(time(1:10:end), arates(1,1:10:end),time(1:10:end), arates(2,1:10:end),time(1:10:end), arates(3,1:10:end));
        grid;
        xlim([0, timeperiod]);
        xlabel('Time [sec]','Fontsize',14);
        ylabel({'Roll Angular', 'Velocity [deg/s]'},'Fontsize',14);        
    end
%     suptitle(sprintf('Ephemerides of %2.2f-second duration',timeperiod))
end

end