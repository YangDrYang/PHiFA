% simulate_meteorix.m

clear all
clf
close all

addpath(genpath('../TianGong1'));

%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 10].*pi./180;
q0 = eul2quat([0, 50, 0]'.*pi./180);
% q0 = SpinCalc('EA321toQ',[0, 50, 0],1e-6,0)';
% epoch = Mjday(2017,8,22,1,0,0);
% orb_iss = initialisation_tle('./inputfiles/iss_tle.txt',epoch);
% [sma,inc,ecc,argp,raan,theta] = state2coe(orb_iss(1:3),orb_iss(4:6));
% % [r,v] = coe2rv(sma*(1-ecc^2)/1000+10, ecc, inc+1.0/180*pi, raan+1.0/180*pi, argp+1.0/180*pi, theta+1.0/180*pi, 0, 0, 0);
% [r,v] = coe2rv(sma*(1-ecc^2)/1000+100, ecc, inc, raan, argp, theta+1.0/180*pi, 0, 0, 0);
% orb_target = [r;v]*1000;
% orb_target = initialisation_tle('./inputfiles/lemur2_tle.txt',epoch);
epoch = 57621.0858224100;
[year, mon, day, hr, mm, sec] = invjday(epoch + 2400000.5);
mat = mat_teme2eci(year, mon, day, hr, mm, sec);%teme to j2000
orb_target(1:3,1) = mat*[-1831.88671586 -6469.93463154 -615.91210989]'*1e3;
orb_target(4:6,1) = mat*[4.761100513 -0.799650983 -6.003467612]'*1e3;
% orb_target = [-1831.88671586 -6469.93463154 -615.91210989 4.761100513 -0.799650983 -6.003467612]'*1e3;
 
Y0 = [orb_target; q0; w0];

initstate = clPropagatorInitialValues();
initstate.proptime = 7200;
initstate.mjd = epoch;
initstate.sec = -1;
initstate.step = 0.1;
initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.n = 100;
initstate.m = 100;

initstate.a_drag = 1; % drag
initstate.g_drag = 1;
initstate.draggradient = 0;%only for area method
initstate.a_srad = 1; % solar radiation pressure
initstate.g_srad = 1;
initstate.a_erad = 1; % earth reflected and emmited radiation pressure: albedo and infrared
initstate.g_erad = 1;
initstate.g_mag = 0; % magnetic torque
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets
initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects
initstate.albedo_model = eAlbedoModel.CERES; % Albedo model



%% create Target
target = loadDSPOSEMeteorix();
% plotTarget(target);
%% some more target conf (needs knowledge of pulselength, wavelength)
% target.hitmethod = eHitMethod.Beam;
target.hitmethod = eHitMethod.Area;

%% create Laserstation
station = createISSLaserstation();

propagator.init(initstate);

sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
sa(3) = sa(1);
sa(4) = sa(1);
sa(5) = sa(1);
target.setSurfaceAttributes(sa);

propagator.initNewSim();
propagator.addTarget(target);

%% Propagation
tic
%         profile on
% odefcn = @RK45; 
% odefcn = @ode23;
% odefcn = @ABM8; 
% propagator.used_integrator_fcn = "RK45";
% eph = propagator.propagateTarget(odefcn,Y0,initstate.step);
propagator.used_integrator_fcn = "RADAU II";
tepochs = 0:initstate.step:initstate.proptime;
eph = propagator.propagateTarget_radau(Y0,tepochs);
toc

%% Write Output
propagator.finishSim();
save2file = sprintf('%s.mat', propagator.outfilename);

if exist('eph', 'var')
    save(save2file, 'propagator', 'initstate', 'eph', 'Y0');
else
    save(save2file, 'propagator', 'initstate', 'Y0');
end
% plotSimulationOutput(initstate,propagator.output_table,'Meteorix');
plotSimulationOutput_PHiFA_radau(initstate,propagator.output_table,'Meteorix');

% load('./logfiles/prop_dspose.mat')
% load('./logfiles/prop_dspose_grav.mat')
load('./logfiles/prop_dspose_grav_drag.mat')
prop_dspose = prop_dspose(:,[1,5:7,2:4,11:14,8:10]);
eph_j2000 = DSPOSE_TEME2J2000(epoch,prop_dspose);
[eph_j2000,eph] = quaternion2angles321(eph_j2000',eph);
diffEph(:,[1:7,11:13]) = eph_j2000([1:7,11:13],:)' - eph([1:7,11:13],:)';
% diffEph(:,8:10) = angdiff(eph_j2000(8:10,:)', eph(8:10,:)');
% diffAng = ephAngdiff(eph_j2000,eph);
diffEph(:,8:10) = ephAngdiff(eph_j2000,eph);

[rotOrbAng_DSPOSE,rotOrbAng_PHiFA,diffRotOrbAng] = ephRotOrbAngdiff(eph_j2000,eph,diag(propagator.rso.moi));

Itensor = diag(propagator.rso.moi);
deg2rad = pi/180;
Hx = Itensor(1) .* eph(11,:).*deg2rad;
Hy = Itensor(2) .* eph(12,:).*deg2rad;
Hz = Itensor(3) .* eph(13,:).*deg2rad;
H_PHiFA = (Hx.^2 + Hy.^2 + Hz.^2).^0.5;
Hx = Itensor(1) .* eph_j2000(11,:).*deg2rad;
Hy = Itensor(2) .* eph_j2000(12,:).*deg2rad;
Hz = Itensor(3) .* eph_j2000(13,:).*deg2rad;
H_DSPOSE = (Hx.^2 + Hy.^2 + Hz.^2).^0.5;
Hx = Itensor(1) .* prop_dspose(:,12);
Hy = Itensor(2) .* prop_dspose(:,13);
Hz = Itensor(3) .* prop_dspose(:,14);
H_DSPOSE2 = (Hx.^2 + Hy.^2 + Hz.^2).^0.5;

figure
subplot(2,2,1)
plot(eph(1,:),diffEph(:,2))
hold on 
plot(eph(1,:),diffEph(:,3))
plot(eph(1,:),diffEph(:,4))
xlabel('Time (Sec)','Fontsize',14)
ylabel('Position Difference (m)','Fontsize',14)
subplot(2,2,2)
plot(eph(1,:),diffEph(:,5))
hold on
plot(eph(1,:),diffEph(:,6))
plot(eph(1,:),diffEph(:,7))
xlabel('Time (Sec)','Fontsize',14)
ylabel('Velocity Difference (m/s)','Fontsize',14)
subplot(2,2,3)
% plot(eph(1,:),diffEph(:,8))
% hold on
% plot(eph(1,:),diffEph(:,9))
% plot(eph(1,:),diffEph(:,10))
plot(eph(1,:),diffRotOrbAng(:,1))
hold on
plot(eph(1,:),diffRotOrbAng(:,2))
plot(eph(1,:),diffRotOrbAng(:,3))
xlabel('Time (Sec)','Fontsize',14)
ylabel('Angules Difference (deg)','Fontsize',14)
subplot(2,2,4)
% plot(eph(1,:),diffEph(:,11))
% hold on
% plot(eph(1,:),diffEph(:,12))
% plot(eph(1,:),diffEph(:,13))
plot(eph(1,:),diffRotOrbAng(:,4))
hold on
plot(eph(1,:),diffRotOrbAng(:,5))
plot(eph(1,:),diffRotOrbAng(:,6))
xlabel('Time (Sec)','Fontsize',14)
ylabel('Angular Rates Difference (deg/s)','Fontsize',14)
% plot3Ephemerides_ecef(eph,station,initstate)
plotSimulationOutput(initstate,propagator.output_table,'Meteorix');
