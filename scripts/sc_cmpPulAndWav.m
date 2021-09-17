% sc_cmpPulAndWav.m
%% load data
sim_data = loadSimDataFolder('pulen_and_walen');

%% divide by pulsed or not
pulsed = zeros(size(sim_data));
for i = 1:length(sim_data)
    if sim_data{i}.propagator.station.bPulsed
        pulsed(i) = 1;
    end
end
npulsed = nnz(pulsed);
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
delpblind = koe(3,:);
delpblinddp = koe(3,end) - koe(3,1);

%% plot pulsed

poip = zeros(npulsed,1);
powp = zeros(npulsed,1);
delpp = zeros(npulsed,length(sim_data{1}.eph));
delppdp = zeros(npulsed,1);

nsim = 1;
for i = 1:length(sim_data)
    if sim_data{i}.propagator.station.bPulsed
        koe = rv2coe_ephemerides(sim_data{i}.eph);
        powp(nsim,:) = sim_data{i}.propagator.station.Wavelength;
        poip(nsim,:) = sim_data{i}.propagator.station.PulseLength;
        delpp(nsim,:) = koe(3,:);
        delppdp(nsim,1) = koe(3,end) - koe(3,1) - delpblinddp;
        nsim = nsim + 1;
    end
end
powp = powp*1E9;
poip = poip*1E9;
%% plot figure
figure;
tri = delaunay(poip, powp);
trisurf(tri, poip, powp, delppdp);
shading interp;
c = colorbar;
c.Label.String = '\Deltaa [m]';
set(gca,'xscale','log');
xlabel('Pulse Length [ns]');
ylabel('Wavelength [nm]');
% zlabel('\Deltaa [m]');
view(2);

plot2tikz('pulwalp',0.6);

