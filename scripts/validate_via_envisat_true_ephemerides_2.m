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
% w0 = [rand*20-10; rand*20-10; rand*20-10].*pi./180;
w0 = [rand*2-1; rand*2-1; rand*2-1].*pi./180;
Y0 = [Y0_dof3; q0; w0];

%% some initial values    
initstate = clPropagatorInitialValues();

initstate.proptime = 7200;
initstate.step = 1;
initstate.mjd = mjd_utc_start;
initstate.sec = -1;

initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.n = 100;
initstate.m = 100;
initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 1; % drag
initstate.g_drag = 1;
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 1;
initstate.a_erad = 1; % earth reflected and emmited radiation pressure: albedo and infrared
initstate.g_erad = 1;
initstate.g_mag = 0; % magnetic torque
initstate.a_sun = 1; % gravity of sun
initstate.a_moon = 1; % gravity of moon
initstate.a_planets = 1; % gravity of other planets
initstate.a_solidEarthTides = 1; % solid earth tides
initstate.a_oceanTides = 1; % ocean tides
initstate.a_relativity = 1; % relativity effects
initstate.albedo_model = eAlbedoModel.CERES; % Albedo model

% initstate.b_logging = 1; % else takes ages
% initstate.n_log_only_every_xth_step = 1;

%% create Target
target = loadDSPOSEEnvisat();
%% some more target conf (needs knowledge of pulselength, wavelength)
% target.hitmethod = eHitMethod.Beam;
target.hitmethod = eHitMethod.Area;

propagator.init(initstate);
station = createOHigginsLaserstation();
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
target.setSurfaceAttributes(sa);
%% Add Objects

propagator.addTarget(target);

%% Propagation
tic
eph = cell(3,1);
%         profile on
for i = 1:3
    if i == 1
        tic
        
        timestamp = now;
        propagator.initNewSim();
        propagator.addTarget(target);
        propagator.used_integrator_fcn = "RADAU II";
        tepochs = 0:initstate.step:initstate.proptime;
        eph_radau = propagator.propagateTarget_radau(Y0,tepochs);
        
        %% Write Output
        save2file = sprintf('%s.mat', propagator.outfilename);
        save(save2file, 'propagator', 'initstate', 'eph_radau', 'Y0', 'timestamp');
        
        eph{i} = eph_radau;
        t_radau = toc;
    elseif i == 2
        tic
        
        propagator.initNewSim();
        propagator.addTarget(target);        
        odefcn = @RK45; 
        eph_rk45 = propagator.propagateTarget(odefcn,Y0,initstate.step);
        
        %% Write Output
        save2file = sprintf('%s.mat', propagator.outfilename);
        save(save2file, 'propagator', 'initstate', 'eph_rk45', 'Y0', 'timestamp');

        eph{i} = eph_rk45;
        t_rk45 = toc;
    else
        tic 
        
        propagator.initNewSim();
        propagator.addTarget(target);        
        
        odefcn = @ABM8;
        eph_abm8 = propagator.propagateTarget(odefcn,Y0,initstate.step);

        %% Write Output
        save2file = sprintf('%s.mat', propagator.outfilename);
        save(save2file, 'propagator', 'initstate', 'eph_abm8', 'Y0', 'timestamp');
        
        eph{i} = eph_abm8;
        t_abm8 = toc;
    end
end
sprintf('%.1f seconds for Raudu II integration\n.',t_radau);
sprintf('%.1f seconds for RK45 integration\n.',t_rk45);
sprintf('%.1f seconds for ABM8 integration\n.',t_abm8);


%% compare
envisat_eph = loadEnvisatMahootiTrueEphemerides();
envisat_eph_eci = ecef2eci_ephemerides(envisat_eph(:,1:121), mjd_utc_start);
orbDiff1 = (envisat_eph_eci(2:7,:) - eph{1,1}(2:7,1:60:end))';
orbDiff2 = (envisat_eph_eci(2:7,:) - eph{2,1}(2:7,1:60:end))';
orbDiff3 = (envisat_eph_eci(2:7,:) - eph{3,1}(2:7,1:60:end))';
% compareEphemerides(envisat_eph_eci, eph{1,1}(:,1:60:end));
% [rotOrbAng_Ref,rotOrbAng_PHiFA,diffRotOrbAng] = ephRotOrbAngdiff(envisat_eph_eci,eph{1},diag(propagator.rso.moi));

norm(orbDiff1(end,1:3))
norm(orbDiff2(end,1:3))
norm(orbDiff3(end,1:3))

angles = zeros(3,length(eph{1}(1,:)));
for j = 1:size(angles,2)
    [angZ, angY, angX] = quat2angle(norm_quat(eph{1}(8:11,j))', 'ZYX' );
    angles(1:3,j) = [angZ, angY, angX].*180/pi;
end

npt = envisat_eph_eci(1,end);

figure
subplot(2,2,1)
plot(eph{1}(1,:),eph{1,1}(2:4,:)')
xlim([0,npt]);
grid on
xlabel('Time of Propagation (Sec)','Fontsize',14);
ylabel({'Position', 'Components (m)'},'Fontsize',14)
legend('x','y','z','Fontsize',14)


subplot(2,2,2)
plot(eph{1}(1,:),eph{1,1}(5:7,:)')
xlim([0,npt]);
grid on
xlabel('Time of Propagation (Sec)','Fontsize',14);
ylabel({'Velocity', 'Components (m/s)'},'Fontsize',14)
legend('v_x','v_y','v_z','Fontsize',14)

subplot(2,2,3)
plot(eph{1}(1,:),angles')
xlim([0,npt]);
grid on
xlabel('Time of Propagation (Sec)','Fontsize',14);
ylabel({'Angular', 'Components (rad)'},'Fontsize',14)
legend('yaw(\phi)','pitch(\theta)','roll(\psi)','Fontsize',14)

subplot(2,2,4)
plot(eph{1}(1,:),eph{1,1}(12:14,:)')
xlim([0,npt]);
grid on
xlabel('Time of Propagation (Sec)','Fontsize',14);
ylabel({'Augular Velocity', 'Components (rad/s)'},'Fontsize',14)
legend('\omega_{\phi}','\omega_{\theta}','\omega_{\psi}','Fontsize',14)
title('6-DOF Propagation','Fontsize',14)


figure
subplot(2,1,1)
plot(envisat_eph_eci(1,:),orbDiff1(:,1:3),'DisplayName','orbDiff1(:,1:3)')
xlim([0,npt]);
grid on
ylabel({'Position', 'Components (m)'},'Fontsize',14)
legen('x','y','z','Fontsize',14)
title('Difference wrt. True Ephemerides','Fontsize',14)
subplot(2,1,2)
plot(envisat_eph_eci(1,:),orbDiff1(:,4:6),'DisplayName','orbDiff1(:,4:6)')
xlim([0,npt]);
grid on
xlabel('Time of Propagation (Sec)','Fontsize',14);
ylabel({'Velocity', 'Components (m/s)'},'Fontsize',14)
legend('v_x','v_y','v_z','Fontsize',14)

