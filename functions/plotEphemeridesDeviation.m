function fig_handles = plotEphemeridesDeviation(varargin)

if nargin<2
    fprintf('At least 2 ephemerides needed.\n');
    return;
end
ephsize = size(varargin{1});
good2go = cell(length(varargin)-1);
fig_handles = cell(length(varargin)-1,1);
for i = 1:length(good2go)
    good2go{i} = varargin{i+1};
    good2go{i}(2:7,:) = good2go{i}(2:7,:) - varargin{1}(2:7,:);
    if ephsize(1)>7
        good2go{i}(12:14,:) = good2go{i}(12:14,:) - varargin{1}(12:14,:);
    end

    if good2go{i}(1,1) == 0
        time = good2go{i}(1,:)./60;
    else
        time = (good2go{i}(1,:)-good2go{i}(1,1))*24*60;
    end
    state = good2go{i}(2:4,:);
    velo = good2go{i}(5:7,:);

    timeperiod = (good2go{i}(1,end) - good2go{i}(1,1))/60;

    ephsize = size(good2go{i});
    if ephsize(1)>7
        rows = 2;
        angles = (quat2eul(quaternion(normalizeQuad(good2go{i}(8:11,:))'))').*180/pi-...
            (quat2eul(quaternion(normalizeQuad(varargin{1}(8:11,:))'))').*180/pi;
        arates = good2go{i}(12:14,:).*180/pi;   
    else
        rows = 1;
    end

    fig_handles{i} = figure;
    title(sprintf('Ephemerides of %2.2f duration',timeperiod))
    subplot(rows,2,1);
    plot(time, state(1,:),time, state(2,:),time, state(3,:));
    xlabel('Minutes [m]');
    ylabel('Position [m]');
    legend('x','y','z');
    grid;
    subplot(rows,2,2);
    plot(time, velo(1,:),time, velo(2,:),time, velo(3,:));
    legend('x','y','z');
    xlabel('Minutes [m]');
    ylabel('Velocity [m/s]');
    grid;
    if ephsize(1)>7
        subplot(rows,2,3);
        plot(time, angles(1,:),time, angles(2,:),time, angles(3,:));
        legend('x','y','z');
        xlabel('Minutes [m]');
        ylabel('Angle [deg]');
        grid;
        subplot(rows,2,4);
        plot(time, arates(1,:),time, arates(2,:),time, arates(3,:));
        legend('x','y','z');
        xlabel('Minutes [m]');
        ylabel('Rotational velocity [deg/s]');
        grid;
    end
end

end