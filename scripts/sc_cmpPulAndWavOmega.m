% sc_cmpPulAndWavOmega.m
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
domegablind = norm(blind_data{1}.eph(12:14,end))-norm(blind_data{1}.eph(12:14,1));

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
        delppdp(nsim,1) = (norm(sim_data{i}.eph(12:14,end)) - ...
            norm(sim_data{i}.eph(12:14,1)) - domegablind)*180/pi;
        nsim = nsim + 1;
    end
end
powp = powp*1E9;
poip = poip*1E9;
%% plot figure
figure;
tri = delaunay(poip, powp);
trisurf(tri, poip, powp, delppdp);
view(2);
shading interp;
c = colorbar;
c.Label.String = '\Delta\omega [°/s]';
set(gca,'xscale','log');
xlabel('Pulse Length [ns]');
ylabel('Wavelength [nm]');
zlabel('\Delta\omega[°/s]');

plot2tikz('pulwalp',0.6);

