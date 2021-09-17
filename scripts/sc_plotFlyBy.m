% sc_plotFlyBy.m

function sc_plotFlyBy(varargin)


%% load vars
filename = sprintf('logfiles//sim%04d.mat',varargin{1});
load(filename);
[koe,perigee] = rv2coe_ephemerides(eph);
[emitted_energy, transfered_energy, pulses, summary] = compareEnergyTransfer(filename);
[elevation, distance] = calculateElevationFromEphemerides(propagator.station, eph);
incAng = zeros(length(propagator.output_table),1);
for i = 1:length(propagator.output_table)
    incAng(i) = propagator.output_table(i).incAng;
end
[pert, prop, vars] = reformatSimulationData(initstate, propagator.output_table);
pertnorm = zeros(1,length(pert));
pertnorm(1,:) = pert(1,:);
for i = 1:length(pert)
    pertnorm(2,i) = norm( pert(2:4,i) );
end

inter = initstate.n_log_only_every_xth_step;
%% image
figure 
yyaxis left
plot(eph(1,1:inter:end-1),transfered_energy(2,:)*1e4);
ylabel('Transfered energy [1e-4J]','Fontsize',14);
yyaxis right
vars(2,vars(2,:)==0) = NaN;
plot(eph(1,1:inter:end-1),vars(2,:),'.');%laser engagement indicator
ylabel('Laser Pulse Epoch','FontSize',14)
set(gca,'ytick',[0 2])
legend('transfered energy','pulse epoch','Fontsize',14);
xlabel('Time [sec]','Fontsize',14);
grid on

figure
yyaxis left
% plot(perigee(1,:),perigee(2,:)/1000);%perigee 
% ylabel('Perigee Height [km]','Fontsize',14);
% plot(elevation(1,:),elevation(2,:)*180/pi,'.');%elevation angle 
% ylabel('Elevaton Angle [deg]','Fontsize',14);
plot(distance(1,:),(pi-incAng(1:end-1)).*180/pi,'.');%incident angle 
ylabel('Interaction Angle [deg]','Fontsize',14);
yyaxis right
plot(distance(1,:),distance(2,:)./1000,'.');%distance 
% legend('cubesat perigee','relative distance','Fontsize',14);
legend('interaction angle','relative distance','Fontsize',14);
xlim([0 7200])
ylabel('Relative Distance [km]','Fontsize',14)
xlabel('Time [sec]','Fontsize',14);
grid on 



% figure
% subplot(4,1,1);
% % yyaxis left
% plot(elevation(2,:).*180./pi);
% ylabel(sprintf('Elevation [deg]'));
% grid;
% hold on;
% % line([276 276],[0 90], 'Color', 'green', ...
% %         'HandleVisibility','off');
% % line([488 488],[0 90], 'Color', 'red', ...
% %         'HandleVisibility','off');
% % ylim([0 60]);
% xlim([0 length(pert)]);
% subplot(4,1,2);
% % yyaxis right
% plot(distance(2,:)./1000,'-');
% grid;
% hold on;
% % line([276 276],[0 4000], 'Color', 'green', ...
% %         'HandleVisibility','off');
% % line([488 488],[0 4000], 'Color', 'red', ...
% %         'HandleVisibility','off');
% xlabel('Time [sec]');
% ylabel(sprintf('Distance [km]'));
% % legend('Elevation', 'Distance', 'Start Engagement', 'End Engagement');
% % legend('Location','northeast')
% xlim([0 length(pert)]);
% 
% subplot(4,1,3);
% % yyaxis left
% plot(pertnorm(2,:));
% hold on;
% % line([276 276],[0 0.04], 'Color', 'green', ...
% %         'HandleVisibility','off');
% % line([488 488],[0 0.04], 'Color', 'red', ...
% %         'HandleVisibility','off');
% ylabel(sprintf('Impulse [Ns]'));
% grid;
% xlim([0 length(pert)]);
% % legend('Location','northwest')
% % yyaxis right
% subplot(4,1,4);
% plot(pulses(2,:),'-');
% hold on;
% % line([276 276],[0 400], 'Color', 'green', ...
% %         'HandleVisibility','off');
% % line([488 488],[0 400], 'Color', 'red', ...
% %         'HandleVisibility','off');
% % grid;
% xlabel('Time [sec]');
% ylabel(sprintf('Pulsecount [-]'));
% % legend('Impulse', 'Pulses');
% % legend('Location','northwest')
% xlim([0 length(pert)]);

plot2tikz('examflyby', 0.8);