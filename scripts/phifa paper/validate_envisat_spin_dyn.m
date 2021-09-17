clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
% w0 = [0; 0; 2.9].*pi./180;
% q0 = eul2quat([0; 0; 0].*pi./180);
Az = 269.22*pi/180;
El = -28.14*pi/180;
dcm = angle2dcm(179.22*pi/180, 0, 61.86*pi/180);%%from the local frame to the body frame
w0 = dcm*[cos(El)*cos(Az); cos(El)*sin(Az); sin(El)]*2*pi/138;
q0 = angle2quat(179.22*pi/180, 0, 61.86*pi/180)';
epoch = 56664.8704;
orb_envisat = initialisation_tle('./inputfiles/envisat_tle.txt',epoch);
Y0 = [orb_envisat; q0; w0];
initstate = clPropagatorInitialValues();
initstate.proptime = 10*24*60*60;
initstate.mjd = 56664.8704;
initstate.year = 2014;
initstate.mon = 1;
initstate.day = 6;
initstate.hour = 20;
initstate.min = 53;
initstate.sec = 0.3760;
initstate.step = 1;
initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 40;
initstate.m = 40;

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
initstate.a_relativity = 1; % relativity effects
initstate.albedo_model = eAlbedoModel.CERES; % Albedo model
propagator.init(initstate);


%% create Target

target = loadDSPOSEEnvisat();
target.hitmethod = eHitMethod.Area;
% plotTarget(target);

%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0;
station.referencePower = 10000;
station.bGroundBased = true;
station.lla = [0, 0, 0];
station.beamResolution = 25;
station.bOnlyInDark = true;
station.intstep = initstate.step;


% %% some more target conf (needs knowledge of pulselength, wavelength)
% sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
% sa(1).bUseDeviation = false;
% sa(2) = sa(1);
% target.setSurfaceAttributes(sa);
% target.hitmethod = eHitMethod.Beam;
    
%% Add Objects

            
propagator.init(initstate);
propagator.addTarget(target);

propagator.initNewSim();

%% Propagation
tic
profile on

%             propagator.used_integrator_fcn = "RADAU II";
%             tepochs = 0:initstate.step:initstate.proptime;
%             eph_radau = propagator.propagateTarget_radau(Y0,tepochs);            
odefcn = @RK45; 
% odefcn = @ode23;
eph_rk45 = propagator.propagateTarget(odefcn,Y0,initstate.step);
toc
            
Az = [];
El = [];
w = [];
sp = [];
for i = 1:size(eph_rk45,2)
%     eul = quat2eul(eph_radau(8:11,i)');
    w = cat(1,w,norm(eph_rk45(12:14,i)));
    dcm = quat2dcm(eph_rk45(8:11,i));
    wl = dcm'*eph_rk45(12:14,i);%%angular velocity in the local frame from the body frame
    Az_tmp = atan2(wl(2),wl(1));
    El_tmp = asin(wl(3)/norm(wl));
    sp_tmp = 2*pi/norm(wl);%spinning period
    Az = cat(1,Az,Az_tmp*180/pi);
    El = cat(1,El,El_tmp*180/pi);
    sp = cat(1,sp,sp_tmp);
end

yy = 2013;
mm = 1;
dd = 1;
hh = 0;
mn = 0;
ss = 0;
mjd_ref = Mjday(yy,mm,dd,hh,mn,ss);
sp_ref = 0.0367320*(initstate.mjd + eph_rk45(1,:)'/86400 - mjd_ref) + 124.883;%%reference values of spinning period
sp_diff = sp - sp_ref;

figure
subplot(2,2,1)
plot(Az(1:60:round(end/2)),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel('Azimuth Angle (degrees)','FontSize',14)
subplot(2,2,2)
plot(El(1:60:round(end/2)),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel('Elevation Angle (degrees)','FontSize',14)
subplot(2,2,3);
plot(sp(1:60:round(end/2)),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel('Spinning Period (seconds)','FontSize',14)
subplot(2,2,4);
plot(sp_diff(1:60:round(end/2)),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel({'Difference of', 'Spinning Periods (seconds)'},'FontSize',14)

figure
subplot(2,2,1)
plot(Az(1:60:end),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel('Azimuth Angle (degrees)','FontSize',14)
subplot(2,2,2)
plot(El(1:60:end),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel('Elevation Angle (degrees)','FontSize',14)
subplot(2,2,3);
plot(sp(1:60:end),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel('Spinning Period (seconds)','FontSize',14)
subplot(2,2,4);
plot(sp_diff(1:60:end),'.')
xlabel('Time (minutes)','FontSize',14)
ylabel({'Difference of', 'Spinning Periods (seconds)'},'FontSize',14)


