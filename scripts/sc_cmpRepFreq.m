% sc_cmpRepFreq.m
%% load data
sim_data = loadSimDataFolder('repfreq');

%% load blind probe
tmpname = sprintf('logfiles/poinacc_and_power/blind.mat');
load(tmpname, 'propagator', 'initstate', 'eph', 'Y0');
blind_data = cell(1,1);
if exist('eph','var')
    blind_data{1} = struct('propagator', propagator,...
        'initstate', initstate,...
        'eph', eph,...
        'Y0', Y0);
end
clear('propagator', 'initstate', 'eph', 'Y0');
koe = rv2coe_ephemerides(blind_data{1}.eph);
dp_blind = koe(3,:);
dp_blind_dp = koe(3,end) - koe(3,1);
time = blind_data{1}.eph(1,:)';

%% plot

dp = zeros(length(sim_data)-1, length(dp_blind));
dp_rel = zeros(length(sim_data)-1, length(dp_blind));
repfreq = zeros(length(sim_data),1);

for i = 1:length(sim_data)
    koe = rv2coe_ephemerides(sim_data{i}.eph);
    repfreq(i) = sim_data{i}.propagator.station.repetitionRate;
    if length(dp)~=length(koe)
        koe3_interp = interp1(sim_data{i}.eph(1,:), koe(3,:), time');
        dp(i,:) = koe3_interp;
        dp_rel(i,:) = koe3_interp - dp_blind;
    else
        dp(i,:) = koe(3,:);
        dp_rel(i,:) = koe(3,:) - dp_blind;
    end
end
%% plot figure
th_fac = 15;
figure;
% tri = delaunay(time, rotvel);
% trisurf(tri, time, rotvel, dp_rel);
% view(2);
% shading interp;
% c = colorbar;
% c.Label.String = '\Deltaa [m]';

time_plot = time - 300;
for i = 1:length(sim_data)
    plot(time_plot(1:th_fac:length(time)), dp_rel(i,1:th_fac:length(time)), getPlotMarkers(i), ...
        'DisplayName',[num2str(repfreq(i)) ' pulses/s']);
    hold on;
end

xlabel('Time after simulationstart [s]');
ylabel('\Deltaa [m]');
legend('Location', 'southwest');
xlim([0 300]);
grid
% zlabel('\Deltaa [m]');

plot2tikz('repfreq',0.6);

