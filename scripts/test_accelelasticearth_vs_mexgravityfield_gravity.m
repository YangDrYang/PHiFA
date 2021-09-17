clear all
%% Create Propagator
% global propagator
propagator = clPropagator.instance();


%% some initial values    
timestamp = now;
w0 = [0; 0; 0].*pi./180;
q0 = angle2quat([0, 0, 0].*pi./180);
% [Y0_dof3, mjd_utc_start, date] = loadEnvisatMahootiInitialState();
Y0 = [getInitialState(0,0,400); q0; w0];
initstate = clPropagatorInitialValues();
initstate.proptime = 1588*60/2;
initstate.step = 1;

initstate.dof = 6;
initstate.a_grav = 1; % harmonic terms of gravity
initstate.g_grav = 1;
initstate.DSPOSE = 1; % using DSPOSE functions       
initstate.n = 80;
initstate.m = 80;

initstate.a_lase = 0; % laser engagements
initstate.g_lase = 0;
initstate.a_drag = 0; % drag
initstate.g_drag = 0;
initstate.a_srad = 0; % solar radiation pressure
initstate.g_srad = 0;
initstate.a_sun = 0; % gravity of sun
initstate.a_moon = 0; % gravity of moon
initstate.a_planets = 0; % gravity of other planets

initstate.a_solidEarthTides = 0; % solid earth tides
initstate.a_oceanTides = 0; % ocean tides
initstate.a_relativity = 0; % relativity effects



%% create Target
% target = createCubeTarget(0.2, false);
target = loadDSPOSEEnvisat();
% plate = makePlate(0.2, 0.05);
% target = createTarget(plate, 'Solid Plate');
%% create Laserstation
station = createOHigginsLaserstation();
station.trackError = 0;
station.pointError = 0.5;
station.referencePower = 10000;
% station.bGroundBased = false;
% station.lla = [0 120 650];

%% some more target conf (needs knowledge of pulselength, wavelength)
sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = true;
sa(2) = sa(1);
target.setSurfaceAttributes(sa);
target.hitmethod = eHitMethod.Beam;

propagator.init(initstate);
propagator.addTarget(target);
propagator.addLaserStation(station);

s_C = size(propagator.Cnm);
s_S = size(propagator.Snm);
C = zeros(101);
S = zeros(101);
C(1:s_C(1),1:s_C(2)) = propagator.Cnm;
S(1:s_S(1),1:s_S(2)) = propagator.Snm;
[C,S]=unnorm_grav_coef(C,S);
%% grav comparison
t = 0:2:360;
mexfd = zeros(3,length(t));
fd = zeros(3,length(t));
dspomex = zeros(3,length(t));
dspo = zeros(3,length(t));
maho = zeros(3,length(t));
gravmex = zeros(3,length(t));
grav = zeros(3,length(t));
gd = zeros(3,length(t));
for i = 1:length(t)
    MJD_UTC = propagator.AuxParam.Mjd_UTC+0/86400;
    target.currentEpoche = MJD_UTC;
    propagator.startmjd = MJD_UTC;
    propagator.t = 0;
    ang = t(i)*pi/180;
    R = [cos(ang), -sin(ang), 0;sin(ang), cos(ang), 0; 0, 0, 1];
    p = R*(Y0(1:3)+[0;0;100]);
    target.xv = [p;0;0;0];
    Y0(7:10) = norm_quat(Y0(7:10)); 
    target.qw = Y0(7:13);
    
    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(propagator.eopdata,MJD_UTC,'l');
    [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    MJD_TT = MJD_UTC+TT_UTC/86400;
    MJD_UT1 = MJD_UTC+UT1_UTC/86400;
    P = PrecMatrix(propagator.const.MJD_J2000,MJD_TT);
    N = NutMatrix(MJD_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(MJD_UT1) * T;

    p_ecef=E*p;
    C_i2b = diag([1,1,1]);
    Inertia = target.moi;
    d = norm(p_ecef);                     % distance
    latgc = asin(p_ecef(3)/d);
    lon = atan2(p_ecef(2),p_ecef(1));
    h = (d - 6378.1366e3)/1000;
    LLA = [latgc,latgc,lon,h];
    
    l_max_a = propagator.AuxParam.n;
    l_max_g = l_max_a;
    n_grav_a = propagator.AuxParam.a_grav;
    n_grav_g = propagator.AuxParam.g_grav;
    
    MJD_TDB = Mjday_TDB(MJD_TT);
    [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
        r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE430(MJD_TDB);

    [dspose_grav, dspose_torque] = ...
                mexGravityField(p_ecef, p, LLA, C_i2b, Inertia, C, S, ...
                l_max_a, l_max_g, n_grav_a, n_grav_g);
            
    mahooti_grav = AccelHarmonic_ElasticEarth(MJD_UTC, [0;0;0], [0;0;0], p, ...
        E, UT1_UTC, TT_UTC, x_pole, y_pole);
    
    [m_grav2, m_dspose_torque] = ...
                AccelGravityGradientTorque(p_ecef, p, LLA, C_i2b, Inertia, C, S, ...
                l_max_a, l_max_g, n_grav_a, n_grav_g);

    maho(:,i) = mahooti_grav;
    dspomex(:,i) = -propagator.const.GM_Earth*p/norm(p)^3 + dspose_grav';
    dspo(:,i) = m_grav2;
    mexfd(:,i) = maho(:,i) - dspomex(:,i);
    fd(:,i) = maho(:,i) - dspo(:,i);
    
    gravmex(:,i) = dspose_torque;
    grav(:,i) = m_dspose_torque;
    gd(:,i) = dspose_torque - m_dspose_torque;
end

figure;
plot(t,maho);
legend('x','y','z');
title('maho');
grid;
figure;
plot(t,dspo);
legend('x','y','z');
title('gravity');
grid;
figure
plot(t,dspomex);
legend('x','y','z');
title('mexdspo');
grid;
figure
plot(t,mexfd);
legend('x','y','z');
title('mexfd');
grid;
figure
plot(t,fd);
legend('x','y','z');
title('fd');
grid;
figure
plot(t,gravmex);
legend('x','y','z');
title('gravmex');
grid;
figure
plot(t,grav);
legend('x','y','z');
title('grav');
grid;
figure
plot(t,gd);
legend('x','y','z');
title('gd');
grid;

