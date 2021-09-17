% simulate_meteorix.m

clear all
clf
clc
close all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();

[Y0_dof3, mjd_utc_start, date] = loadEnvisatMahootiInitialState();
q0 = angle2quat(rand*360*pi/180, rand*360*pi/180, rand*360*pi/180,'ZYX')';
w0 = [rand*6-3; rand*6-3; rand*6-3].*pi./180;
Y0 = [Y0_dof3; q0; w0];

initstate = clPropagatorInitialValues();
initstate.mjd = mjd_utc_start;
initstate.sec = -1;
initstate.step = 1;
initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.n = 20;
initstate.m = 20;
initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 1; % drag
initstate.g_drag = 1;
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 1;
initstate.a_erad = 0; % earth reflected and emmited radiation pressure: albedo and infrared
initstate.g_erad = 0;
initstate.g_mag = 0; % magnetic torque
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets
initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects
initstate.albedo_model = eAlbedoModel.CERES; % Albedo model
propagator.init(initstate);

%% create Target
target = loadDSPOSEEnvisat();
%% some more target conf (needs knowledge of pulselength, wavelength)
% target.hitmethod = eHitMethod.Beam;
target.hitmethod = eHitMethod.Area;
% plotTarget(target);

propagator.init(initstate);
station = createOHigginsLaserstation();
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
target.setSurfaceAttributes(sa);
%% Add Objects
propagator.initNewSim();
propagator.addTarget(target);


envisat_eph = loadEnvisatMahootiTrueEphemerides();
envisat_eph_eci = ecef2eci_ephemerides(envisat_eph(:,1:1441), initstate.mjd);
 
%% Propagation
tic
orbDiff = [];
t_exe = [];
profile off
for i = 2:3
    if i == 1
        propagator.used_integrator_fcn = "RADAU II";
        for j = 1:0.5:6
            
            initstate.proptime = j*60*60;
            tepochs = 0:initstate.step:initstate.proptime;
            tic
            eph = propagator.propagateTarget_radau(Y0,tepochs);
            toc
            t_radau = toc;
            eph_length = size(eph,2);
            ind = (eph_length - 1)/(60/initstate.step) + 1;
            orbDiff = cat(1,orbDiff,[norm(envisat_eph_eci(2:4,ind) - eph(2:4,end));envisat_eph_eci(2:7,ind) - eph(2:7,end)]');
            t_exe = cat(1,t_exe,t_radau);
            sprintf('%.1f seconds for Raudu II integration\n.',t_radau);            
        end
        
    elseif i == 2
        
        odefcn = @RK45; 
        
        for j = 1:0.5:6
            initstate.proptime = j*60*60;
            propagator.init(initstate);
            tic
            eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
            toc
            t_rk45 = toc;
            eph_length = size(eph,2);
            ind = (eph_length - 1)/(60/initstate.step) + 1;            
            orbDiff = cat(1,orbDiff,[norm(envisat_eph_eci(2:4,ind) - eph(2:4,end));envisat_eph_eci(2:7,ind) - eph(2:7,end)]');
            t_exe = cat(1,t_exe,t_rk45);
            sprintf('%.1f seconds for RK45 integration\n.',t_rk45);
        end
        
    else
        tic       
        odefcn = @ABM8;
        
        for j = 1:0.5:6
            initstate.proptime = j*60*60;
            propagator.init(initstate);
            tic
            eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
            toc
            t_abm8 = toc;
            eph_length = size(eph,2);
            ind = (eph_length - 1)/(60/initstate.step) + 1;            
            orbDiff = cat(1,orbDiff,[norm(envisat_eph_eci(2:4,ind) - eph(2:4,end));envisat_eph_eci(2:7,ind) - eph(2:7,end)]');
            t_exe = cat(1,t_exe,t_abm8);
            sprintf('%.1f seconds for ABM8 integration\n.',t_abm8);
        end

    end
end

% figure
% subplot(1,2,1)
% yyaxis left
% plot(1:0.5:6,t_exe(1:11),'o');
% axis([0.5,6.5,100,1800])
% ylabel('Time Cost (s)','FontSize',14)
% yyaxis right
% plot(1:0.5:6,t_exe(12:22),'>');
% axis([0.5,6.5,100,1800])
% xlabel('Integration Time (h)','FontSize',14)
% set(gca,'XTick',0.5:0.5:6.5);
% ylabel('Time Cost (s)','FontSize',14)
% legend('RK45','ABM8','FontSize',14)
% 
% subplot(1,2,2)
% yyaxis left
% plot(1:0.5:6,orbDiff(1:11,1),'o');
% axis([0.5,6.5,10,60])
% ylabel('3D Orbit Difference (m)','FontSize',14)
% yyaxis right
% plot(1:0.5:6,orbDiff(12:22,1),'>');
% axis([0.5,6.5,10,60])
% set(gca,'XTick',0.5:0.5:6.5);
% xlabel('Integration Time (h)','FontSize',14)
% ylabel('3D Orbit Difference (m)','FontSize',14)
% legend('RK45','ABM8','FontSize',14)

figure
subplot(1,2,1)
plot(1:0.5:6,t_exe(1:11),'o');
axis([0.5,6.5,100,1800])
ylabel('Time Cost (s)','FontSize',14)
hold on
plot(1:0.5:6,t_exe(12:22),'>');
xlabel('Integration Time (h)','FontSize',14)
set(gca,'XTick',0.5:0.5:6.5);
ylabel('Time Cost (s)','FontSize',14)
legend('RK45','ABM8','FontSize',14)

subplot(1,2,2)
plot(1:0.5:6,orbDiff(1:11,1),'o');
axis([0.5,6.5,10,60])
ylabel('3D Orbit Difference (m)','FontSize',14)
hold on
plot(1:0.5:6,orbDiff(12:22,1),'>');
axis([0.5,6.5,10,60])
set(gca,'XTick',0.5:0.5:6.5);
xlabel('Integration Time (h)','FontSize',14)
ylabel('3D Orbit Difference (m)','FontSize',14)
legend('RK45','ABM8','FontSize',14)