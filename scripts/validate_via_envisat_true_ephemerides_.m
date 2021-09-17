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

%% some initial values    
initstate = clPropagatorInitialValues();

initstate.proptime = 1588*60;
initstate.step = 10;
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
initstate.g_mag = 1; % magnetic torque
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

propagator.initNewSim();
propagator.addTarget(target);

%% Propagation
tic
eph = cell(3,1);
%         profile on
for i = 1:3
    if i == 1
        tic
        propagator.used_integrator_fcn = "RADAU II";
        tepochs = 0:initstate.step:initstate.proptime;
        eph{i} = propagator.propagateTarget_radau(Y0,tepochs);
        t_radau = toc;
    elseif i == 2
        tic
        odefcn = @RK45; 
        eph{i} = propagator.propagateTarget(odefcn,Y0,initstate.step);
        t_rk45 = toc;
    else
        odefcn = @ABM8;
        eph{i} = propagator.propagateTarget(odefcn,Y0,initstate.step);
        t_abm8 = toc;
    end
    
end
sprintf('%.1f seconds for Raudu II integration\n.',t_radau);
sprintf('%.1f seconds for RK45 integration\n.',t_rk45);
sprintf('%.1f seconds for ABM8 integration\n.',t_abm8);

tt = eph{1}(1,:);
dd1 = eph{1} - eph{2};
dd2 = eph{3} - eph{2};
dd3 = eph{3} - eph{1};

ephsize = size(eph{1},1);
angles1 = zeros(3, size(eph{1},2));
arates1 = zeros(3, size(eph{1},2));
for j = 1:size(eph{1},2)
    angles1(:,j) = eulerd(quaternion(eph{1}(8:11,j)'),'ZYX','frame');%radian -> degree
%     angles1(:,j) = (Q2Eul(norm_quat(eph{1}(8:11,j)))).*180/pi;
    arates1(:,j) = eph{1}(12:14,j).*180/pi;   
end
angles2 = zeros(3, size(eph{2},2));
arates2 = zeros(3, size(eph{2},2));
for j = 1:size(eph{2},2)
    angles2(:,j) = eulerd(quaternion(eph{2}(8:11,j)'),'ZYX','frame');%radian -> degree
%     angles2(:,j) = (Q2Eul(norm_quat(eph{2}(8:11,j)))).*180/pi;
    arates2(:,j) = eph{2}(12:14,j).*180/pi;   
end
angles3 = zeros(3, size(eph{3},2));
arates3 = zeros(3, size(eph{3},2));
for j = 1:size(eph{3},2)
    angles3(:,j) = eulerd(quaternion(eph{3}(8:11,j)'),'ZYX','frame');%radian -> degree    
%     angles3(:,j) = (Q2Eul(norm_quat(eph{3}(8:11,j)))).*180/pi;
    arates3(:,j) = eph{3}(12:14,j).*180/pi;   
end
    

% plotEphemerides(eph{1});
% plotEphemerides(eph{2});
% plotEphemerides(eph{3});

ddang1 = angles1-angles2;
figure(1)
subplot(2,2,1)
plot(tt,dd1(2,:));
hold on
plot(tt,dd1(3,:));
plot(tt,dd1(4,:));
subplot(2,2,2)
plot(tt,dd1(5,:));
hold on
plot(tt,dd1(6,:));
plot(tt,dd1(7,:));
subplot(2,2,3)
plot(tt,ddang1(1,:));
hold on
plot(tt,ddang1(2,:));
plot(tt,ddang1(3,:));
subplot(2,2,4)
plot(tt,dd1(12,:));
plot(tt,dd1(13,:));
plot(tt,dd1(14,:));

ddang2 = angles3-angles2;
figure(2)
subplot(2,2,1)
plot(tt,dd2(2,:));
hold on
plot(tt,dd2(3,:));
plot(tt,dd2(4,:));
subplot(2,2,2)
plot(tt,dd2(5,:));
hold on
plot(tt,dd2(6,:));
plot(tt,dd2(7,:));
subplot(2,2,3)
plot(tt,ddang2(1,:));
hold on
plot(tt,ddang2(2,:));
plot(tt,ddang2(3,:));
subplot(2,2,4)
plot(tt,dd2(12,:));
plot(tt,dd2(13,:));
plot(tt,dd2(14,:));

ddang3 = angles3-angles1;
figure(3)
subplot(2,2,1)
plot(tt,dd3(2,:));
hold on
plot(tt,dd3(3,:));
plot(tt,dd3(4,:));
subplot(2,2,2)
plot(tt,dd3(5,:));
hold on
plot(tt,dd3(6,:));
plot(tt,dd3(7,:));
subplot(2,2,3)
plot(tt,ddang3(1,:));
hold on
plot(tt,ddang3(2,:));
plot(tt,ddang3(3,:));
subplot(2,2,4)
plot(tt,dd3(12,:));
plot(tt,dd3(13,:));
plot(tt,dd3(14,:));
toc
