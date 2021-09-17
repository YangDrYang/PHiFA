%% load vars
load logfiles\sim0207.mat
koe = rv2coe_ephemerides(eph);
[emitted_energy, transfered_energy, pulses, summary] = compareEnergyTransfer(207);
elevation = calculateElevationFromEphemerides(propagator.station, eph);
r_p = (1-koe(2,:)).*koe(3,:);

%% cutting

transfered_energy = transfered_energy(:,57850:58100);
elevation = elevation(:,57850:58100);
%% image
figure
yyaxis left
plot(elevation(2,:).*180./pi);
xlim([0 250]);
set(gca,'FontSize',16)
yyaxis right
semilogy(transfered_energy(2,:),'-');
set(gca,'FontSize',16)
hold on
semilogy(transfered_energy(4,:),'m-');
set(gca,'FontSize',16)
grid
xline(24);
xline(239);
ylabel('Energy change due to radiation [J/s]');
xlabel('Time [sec]');
yyaxis left
ylabel('Pulsenumber of fly-by [-]');
legend('Elevation', 'Pulsed', 'CW');
legend('Location','northwest')
yyaxis left
ylabel('Elevation [deg]');