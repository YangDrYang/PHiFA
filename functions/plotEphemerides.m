function plotEphemerides(varargin)

for i = 1:length(varargin)
    time = varargin{i}(1,:);
    state = varargin{i}(2:4,:)./1000;
    velo = varargin{i}(5:7,:)./1000;

    timeperiod = varargin{i}(1,end) - varargin{i}(1,1) + 1;

    ephsize = size(varargin{i},1);
    if ephsize(1)>7
        rows = 2;
        angles = zeros(3, size(varargin{i},2));
        arates = zeros(3, size(varargin{i},2));
        for j = 1:size(varargin{i},2) 
            angles(:,j) = eulerd(quaternion(varargin{i}(8:11,j)'),'ZYX','frame')./180*pi;%degree -> radian
%               angles(:,j) = (Q2Eul(norm_quat(varargin{i}(8:11,j)))).*180/pi;
            arates(:,j) = varargin{i}(12:14,j);   
            
        end
    else
        rows = 1;
    end

    figure;
    subplot(rows,2,1);
    plot(time, state(1,:),time, state(2,:),time, state(3,:));
    grid;
    xlabel('Time [sec]','Fontsize',14);
    ylabel('Position [km]','Fontsize',14);
%     legend('x','y','z');
    subplot(rows,2,2);
    plot(time, velo(1,:),time, velo(2,:),time, velo(3,:));
    grid;
    legend('x','y','z','Fontsize',14);
    xlabel('Time [sec]','Fontsize',14);
    ylabel('Velocity [km/s]','Fontsize',14);
%     axis([0 100 -6.5 6])
    if ephsize(1)>7
        subplot(rows,2,3);
        plot(time, angles(1,:),time, angles(2,:),time, angles(3,:));
        grid;
%         legend('x','y','z');
        xlabel('Time [sec]','Fontsize',14);
        ylabel('Angle [rad]','Fontsize',14);
        subplot(rows,2,4);
%         plot(time, arates(1,:),time, arates(2,:),time, arates(3,:));
        plot(time(1:10:end), arates(1,1:10:end),time(1:10:end), arates(2,1:10:end),time(1:10:end), arates(3,1:10:end));
        grid;
        legend('yaw','pitch','roll','Fontsize',14);
        xlabel('Time [sec]','Fontsize',14);
        ylabel('Angular velocity [rad/s]','Fontsize',14);
    end
%     suptitle(sprintf('Ephemerides of %2.2f-second duration',timeperiod))
end

end
