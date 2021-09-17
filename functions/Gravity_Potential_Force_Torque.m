function [gravity_potential, a_gravity_inertial, g_gravity_body] = Gravity_Potential_Force_Torque(Mjd_UTC,r_Sun,r_Moon,...
    p,E,UT1_UTC,TT_UTC,x_pole,y_pole,C_i2b,Inertia,mass)
% (Mjd_UTC,r_Sun,r_Moon,r,E,UT1_UTC,TT_UTC,x_pole,y_pole)

prpgtr = clPropagator.instance();

r_ref = 6378.1366e3;   % Earth's radius [m]; ITG-Grace03
gm    = 398600.4415e9; % [m^3/s^2]; ITG-Grace03

C = prpgtr.Cnm;
S = prpgtr.Snm;
if (prpgtr.AuxParam.a_solidEarthTides) || (prpgtr.AuxParam.a_oceanTides)
    r_Moon = E*r_Moon;
    [lM, phiM, rM] = CalcPolarAngles(r_Moon);
    r_Sun = E*r_Sun;
    [lS, phiS, rS] = CalcPolarAngles(r_Sun);

    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;

    T  = (Mjd_TT-prpgtr.const.MJD_J2000)/36525;
    T2 = T*T;
    T3 = T2*T;
    rev = 360*3600;  % arcsec/revolution
end

if (prpgtr.AuxParam.a_solidEarthTides)
    % Effect of Solid Earth Tides (elastic Earth)
    % For dC21 and dS21
    % The coefficients we choose are in-phase(ip) amplitudes and out-of-phase amplitudes of the
    % corrections for frequency dependence, and multipliers of the Delaunay variables
    % Refer to Table 6.5a in IERS2010
    coeff0 = [...
	%  l   l'  F   D   Om  Amp(R) Amp(I)
       2,  0,  2,  0,  2,  -0.1,    0;    
       0,  0,  2,  2,  2,  -0.1,    0;    
       1,  0,  2,  0,  1,  -0.1,    0;    
       1,  0,  2,  0,  2,  -0.7,    0.1;  
      -1,  0,  2,  2,  2,  -0.1,    0;    
       0,  0,  2,  0,  1,  -1.3,    0.1;  
       0,  0,  2,  0,  2,  -6.8,    0.6;  
       0,  0,  0,  2,  0,   0.1,    0;    
       1,  0,  2, -2,  2,   0.1,    0;    
      -1,  0,  2,  0,  1,   0.1,    0;    
      -1,  0,  2,  0,  2,   0.4,    0;    
       1,  0,  0,  0,  0,   1.3,   -0.1;  
       1,  0,  0,  0,  1,   0.3,    0;    
      -1,  0,  0,  2,  0,   0.3,    0;    
      -1,  0,  0,  2,  1,   0.1,    0;    
       0,  1,  2, -2,  2,  -1.9,    0.1;  
       0,  0,  2, -2,  1,   0.5,    0;    
       0,  0,  2, -2,  2,  -43.4,   2.9;
       0, -1,  2, -2,  2,   0.6,    0;    
       0,  1,  0,  0,  0,   1.6,   -0.1;  
      -2,  0,  2,  0,  1,   0.1,    0;    
       0,  0,  0,  0, -2,   0.1,    0;    
       0,  0,  0,  0, -1,  -8.8,    0.5;  
       0,  0,  0,  0,  0,   470.9, -30.2; 
       0,  0,  0,  0,  1,   68.1,  -4.6;  
       0,  0,  0,  0,  2,  -1.6,    0.1;  
      -1,  0,  0,  1,  0,   0.1,    0;    
       0, -1,  0,  0, -1,  -0.1,    0;    
       0, -1,  0,  0,  0,  -20.6,  -0.3;  
       0,  1, -2,  2, -2,   0.3,    0;    
       0, -1,  0,  0,  1,  -0.3,    0;    
      -2,  0,  0,  2,  0,  -0.2,    0;    
      -2,  0,  0,  2,  1,  -0.1,    0;    
       0,  0, -2,  2, -2,  -5.0,    0.3;  
       0,  0, -2,  2, -1,   0.2,    0;    
       0, -1, -2,  2, -2,  -0.2,    0;    
       1,  0,  0, -2,  0,  -0.5,    0;    
       1,  0,  0, -2,  1,  -0.1,    0;    
      -1,  0,  0,  0, -1,   0.1,    0;    
      -1,  0,  0,  0,  0,  -2.1,    0.1;  
      -1,  0,  0,  0,  1,  -0.4,    0;    
       0,  0,  0, -2,  0,  -0.2,    0;    
      -2,  0,  0,  0,  0,  -0.1,    0;    
       0,  0, -2,  0, -2,  -0.6,    0;    
       0,  0, -2,  0, -1,  -0.4,    0;    
       0,  0, -2,  0,  0,  -0.1,    0;    
      -1,  0, -2,  0, -2,  -0.1,    0;    
      -1,  0, -2,  0, -1,  -0.1,    0;    
    ];
    % For dC20
    % The nominal value k20 for the zonal tides is taken as 0.30190
	% Refer to Table 6.5b in IERS2010
    coeff1 = [...
	% l   l'  F   D   Om  Amp(R)  Amp(I)
      0,  0,  0,  0,  1,  16.6,   -6.7;
      0,  0,  0,  0,  2,  -0.1,    0.1;
      0, -1,  0,  0,  0,  -1.2,    0.8;
      0,  0, -2,  2, -2,  -5.5,    4.3;
      0,  0, -2,  2, -1,   0.1,   -0.1;
      0, -1, -2,  2, -2,  -0.3,    0.2;
      1,  0,  0, -2,  0,  -0.3,    0.7;
     -1,  0,  0,  0, -1,   0.1,   -0.2;
     -1,  0,  0,  0,  0,  -1.2,    3.7;
     -1,  0,  0,  0,  1,   0.1,   -0.2;
      1,  0, -2,  0, -2,   0.1,   -0.2;
      0,  0,  0, -2,  0,   0.0,    0.6;
     -2,  0,  0,  0,  0,   0.0,    0.3;
      0,  0, -2,  0, -2,   0.6,    6.3;
      0,  0, -2,  0, -1,   0.2,    2.6;
      0,  0, -2,  0,  0,   0.0,    0.2;
      1,  0, -2, -2, -2,   0.1,    0.2;
     -1,  0, -2,  0, -2,   0.4,    1.1;
     -1,  0, -2,  0, -1,   0.2,    0.5;
      0,  0, -2, -2, -2,   0.1,    0.2;
     -2,  0, -2,  0, -2,   0.1,    0.1;
    ];
    % For dC22 and dS22
    % Refer to Table 6.5c in IERS2010
    coeff2 = [...
    % l  l' F  D  Om   Amp
      1, 0, 2, 0, 2,  -0.3;
      0, 0, 2, 0, 2,  -1.2;
    ];
    
    % Mean arguments of luni-solar motion
    %
    %   l   mean anomaly of the Moon
    %   l'  mean anomaly of the Sun
    %   F   mean argument of latitude
    %   D   mean longitude elongation of the Moon from the Sun
    %   Om  mean longitude of the ascending node of the Moon
    l  = mod (  485866.733 + (1325.0*rev +  715922.633)*T...
                              + 31.310*T2 + 0.064*T3, rev );
    lp = mod ( 1287099.804 + (  99.0*rev + 1292581.224)*T...
                              -  0.577*T2 - 0.012*T3, rev );
    F  = mod (  335778.877 + (1342.0*rev +  295263.137)*T...
                              - 13.257*T2 + 0.011*T3, rev );
    D  = mod ( 1072261.307 + (1236.0*rev + 1105601.328)*T...
                              -  6.891*T2 + 0.019*T3, rev );
    Om = mod (  450160.280 - (   5.0*rev +  482890.539)*T...
                              +  7.455*T2 + 0.008*T3, rev );
    
    % STEP1 CORRECTIONS
    [lgM, dlgM] = Legendre(4,4,phiM);
    [lgS, dlgS] = Legendre(4,4,phiS);
    dCnm20 = (0.29525/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = (0.29470/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = (0.29470/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*(sin(lM))...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*(sin(lS)) );
    dCnm22 = (0.29801/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm22 = (0.29801/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*(sin(2*lM))...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*(sin(2*lS)) );
    dCnm30 = (0.093/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = (0.093/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = (0.093/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = (0.093/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = (0.093/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = (0.094/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = (0.094/7)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = (-0.00087/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,1) );
    dCnm41 = (-0.00079/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,2)*cos(lS) );
    dSnm41 = (-0.00079/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,2)*sin(lS) );
    dCnm42 = (-0.00057/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,3)*cos(2*lS) );
    dSnm42 = (-0.00057/5)*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(5,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(5,3)*sin(2*lS) );
    
    % STEP2 CORRECTIONS
    dC20 = 0;
    for i=1:21
        theta_f = -coeff1(i,1:5)*[l lp F D Om]';
        dC20 = dC20 + 1e-12*(coeff1(i,6)*cos(theta_f)-coeff1(i,7)*sin(theta_f));
    end
    dCnm20 = dCnm20 + dC20;
    
    theta_g = gmst(Mjd_UT1);
    dC21 = 0;
    dS21 = 0;
    for i=1:48
        theta_f = (theta_g+pi)-coeff0(i,1:5)*[l lp F D Om]';
        dC21 = dC21 + 1e-12*coeff0(i,6)*sin(theta_f);
        dS21 = dS21 + 1e-12*coeff0(i,6)*cos(theta_f);
    end
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    dC22 = 0;
    dS22 = 0;
    for i=1:2
        theta_f = 2*(theta_g+pi)-coeff2(i,1:5)*[l lp F D Om]';
		dC22 = dC22 + 1e-12*coeff2(i,6)*sin(theta_f);
        dS22 = dS22 + 1e-12*coeff2(i,6)*cos(theta_f);
    end
    dCnm22 = dCnm22 + dC22;
    dSnm22 = dSnm22 + dS22;
    
    % Treatment of the Permanent Tide (elastic Earth)
    dC20 = 4.4228e-8*(-0.31460)*0.29525;
    dCnm20 = dCnm20 - dC20;
    
    % Effect of Solid Earth Pole Tide (elastic Earth)
    dC21 = -1.290e-9*(x_pole);
    dS21 = 1.290e-9*(y_pole);
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;    
end

if (prpgtr.AuxParam.a_oceanTides)
    % Ocean Tides
    [lgM, dlgM] = Legendre(6,6,phiM);
    [lgS, dlgS] = Legendre(6,6,phiS);
    
    dCnm20 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = -0.3075/5*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
    dCnm22 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm22 = -0.3075/5*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );
    dCnm30 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = -0.195/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = -0.195/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = -0.195/7*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,1) );
    dCnm41 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,2)*cos(lS) );
    dSnm41 = -0.132/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,2)*sin(lS) );
    dCnm42 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,3)*cos(2*lS) );
    dSnm42 = -0.132/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,3)*sin(2*lS) );
    dCnm43 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,4)*cos(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,4)*cos(3*lS) );
    dSnm43 = -0.132/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,4)*sin(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,4)*sin(3*lS) );
    dCnm44 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,5)*cos(4*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,5)*cos(4*lS) );
    dSnm44 = -0.132/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,5)*sin(4*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,5)*sin(4*lS) );
    dCnm50 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,1) );
    dCnm51 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,2)*cos(lS) );
    dSnm51 = -0.1032/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,2)*sin(lS) );
    dCnm52 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,3)*cos(2*lS) );
    dSnm52 = -0.1032/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,3)*sin(2*lS) );
    dCnm53 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,4)*cos(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,4)*cos(3*lS) );
    dSnm53 = -0.1032/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,4)*sin(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,4)*sin(3*lS) );
    dCnm54 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,5)*cos(4*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,5)*cos(4*lS) );
    dSnm54 = -0.1032/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,5)*sin(4*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,5)*sin(4*lS) );
    dCnm55 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,6)*cos(5*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,6)*cos(5*lS) );
    dSnm55 = -0.1032/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,6)*sin(5*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,6)*sin(5*lS) );
    dCnm60 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,1)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,1) );
    dCnm61 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,2)*cos(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,2)*cos(lS) );
    dSnm61 = -0.0892/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,2)*sin(lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,2)*sin(lS) );
    dCnm62 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,3)*cos(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,3)*cos(2*lS) );
    dSnm62 = -0.0892/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,3)*sin(2*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,3)*sin(2*lS) );
    dCnm63 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,4)*cos(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,4)*cos(3*lS) );
    dSnm63 = -0.0892/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,4)*sin(3*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,4)*sin(3*lS) );
    dCnm64 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,5)*cos(4*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,5)*cos(4*lS) );
    dSnm64 = -0.0892/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,5)*sin(4*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,5)*sin(4*lS) );
    dCnm65 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,6)*cos(5*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,6)*cos(5*lS) );
    dSnm65 = -0.0892/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,6)*sin(5*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,6)*sin(5*lS) );
    dCnm66 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,7)*cos(6*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,7)*cos(6*lS) );
    dSnm66 = -0.0892/9*( (prpgtr.const.GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,7)*sin(6*lM)...
           + (prpgtr.const.GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,7)*sin(6*lS) );
    
    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    C(5,4) = C(5,4) + dCnm43;
    C(5,5) = C(5,5) + dCnm44;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;
    S(5,4) = S(5,4) + dSnm43;
    S(5,5) = S(5,5) + dSnm44;
    
    C(6,1) = C(6,1) + dCnm50;
    C(6,2) = C(6,2) + dCnm51;
    C(6,3) = C(6,3) + dCnm52;
    C(6,4) = C(6,4) + dCnm53;
    C(6,5) = C(6,5) + dCnm54;
    C(6,6) = C(6,6) + dCnm55;
    S(6,2) = S(6,2) + dSnm51;
    S(6,3) = S(6,3) + dSnm52;
    S(6,4) = S(6,4) + dSnm53;
    S(6,5) = S(6,5) + dSnm54;
    S(6,6) = S(6,6) + dSnm55;
    
    C(7,1) = C(7,1) + dCnm60;
    C(7,2) = C(7,2) + dCnm61;
    C(7,3) = C(7,3) + dCnm62;
    C(7,4) = C(7,4) + dCnm63;
    C(7,5) = C(7,5) + dCnm64;
    C(7,6) = C(7,6) + dCnm65;
    C(7,7) = C(7,7) + dCnm66;
    S(7,2) = S(7,2) + dSnm61;
    S(7,3) = S(7,3) + dSnm62;
    S(7,4) = S(7,4) + dSnm63;
    S(7,5) = S(7,5) + dSnm64;
    S(7,6) = S(7,6) + dSnm65;
    S(7,7) = S(7,7) + dSnm66;    
end

% Body-fixed position 
p_ecef = E * p;
% p_ecef = r;

% Auxiliary quantities
r = norm(p_ecef);                     % distance
latgc = asin(p_ecef(3)/r);
lon = atan2(p_ecef(2),p_ecef(1));

[pnm, dpnm] = Legendre(prpgtr.AuxParam.n+1,prpgtr.AuxParam.m+1,latgc);

gravity_potential = 0;
g_gravity_body = [0;0;0];    
a_gravity_inertial = [0;0;0];

%%/ Gravity Potential and ASPHERICAL ACCELERATION CALCULATION
if (prpgtr.AuxParam.a_grav)
    U = 0;
    dUdr = 0;
    dUdlatgc = 0;
    dUdlon = 0;
    q3 = 0; q2 = q3; q1 = q2;
    for n=0:prpgtr.AuxParam.n
        b1 = (-gm/r^2)*(r_ref/r)^n*(n+1);
        b2 =  (gm/r)*(r_ref/r)^n;
        b3 =  (gm/r)*(r_ref/r)^n;
        for m=0:prpgtr.AuxParam.m
            q1 = q1 + pnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
            q2 = q2 + dpnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
            q3 = q3 + m*pnm(n+1,m+1)*(S(n+1,m+1)*cos(m*lon)-C(n+1,m+1)*sin(m*lon));
        end
        U = U + q1*b2;
        dUdr     = dUdr     + q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;
        q3 = 0; q2 = q3; q1 = q2;
    end
    U = -U*mass;

    % Body-fixed acceleration
    r2xy = p_ecef(1)^2+p_ecef(2)^2;

    ax = (1/r*dUdr-p_ecef(3)/(r^2*sqrt(r2xy))*dUdlatgc)*p_ecef(1)-(1/r2xy*dUdlon)*p_ecef(2);
    ay = (1/r*dUdr-p_ecef(3)/(r^2*sqrt(r2xy))*dUdlatgc)*p_ecef(2)+(1/r2xy*dUdlon)*p_ecef(1);
    az =  1/r*dUdr*p_ecef(3)+sqrt(r2xy)/r^2*dUdlatgc;

    a_bf = [ax, ay, az]';
    a_gravity_inertial = E'*a_bf;
    gravity_potential = U;
end


%%% TORQUE CALCULATION
if (prpgtr.AuxParam.g_grav)
    % Matrices

    rr = r*r;
    rrr = rr*r;

    i2j2 = p(1)*p(1)+p(2)*p(2);

    d2rdr2(1:3,1:3) = [...
        p(1)*p(1)-rr, p(1)*p(2),      p(1)*p(3);...
        p(1)*p(2),     p(2)*p(2)-rr,  p(2)*p(3);...
        p(1)*p(3),     p(2)*p(3),      p(3)*p(3)-rr];

    d2rdr2 = -d2rdr2./(rrr);

    d20dr2(1:3,1:3) = [...
        p(3)*(2*p(1)^4+p(1)*p(1)*p(2)*p(2)-p(2)^4-p(2)*p(2)*p(3)*p(3)),...
            p(1)*p(2)*p(3)*(3*i2j2+p(3)*p(3)),...
            -p(1)*i2j2*(i2j2-p(3)*p(3));...
        p(1)*p(2)*p(3)*(3*i2j2+p(3)*p(3)),...
            p(3)*(-p(1)^4+p(1)*p(1)*p(2)*p(2)-p(1)*p(1)*p(3)*p(3)+2*p(2)^4),...
            -p(2)*i2j2*(i2j2-p(3)*p(3));...
        -p(1)*i2j2*(i2j2-p(3)*p(3)),...
            -p(2)*i2j2*(i2j2-p(3)*p(3)),...
            -2*p(3)*i2j2*i2j2];

    d20dr2 = d20dr2./(r^4*i2j2^(3/2.0));

    d2ldr2(1:3,1:3) = [...
        2*p(1)*p(2),           p(2)*p(2)-p(1)*p(1),    0;...
        p(2)*p(2)-p(1)*p(1),   -2*p(1)*p(2),           0;...
        0,                     0,                      0;];

    d2ldr2 = d2ldr2./(i2j2*i2j2);

    drdrT(1:3,1:3) = [ ...
        p(1)*p(1), p(1)*p(2),  p(1)*p(3); ...
        p(1)*p(2), p(2)*p(2),  p(2)*p(3); ...
        p(1)*p(3), p(3)*p(2),  p(3)*p(3)];


    drdrT = drdrT./(rr);


    d0d0T(1:3,1:3) = [ ...
        p(1)*p(1)*p(3)*p(3),   p(1)*p(2)*p(3)*p(3),    -p(1)*p(3)*i2j2;
        p(1)*p(2)*p(3)*p(3),   p(2)*p(2)*p(3)*p(3),    -p(2)*p(3)*i2j2;
        -p(1)*p(3)*i2j2,       -p(2)*p(3)*i2j2,        i2j2*i2j2];


    d0d0T = d0d0T./(rrr*r*i2j2*i2j2);


    dldlT(1:3,1:3) = [ ... 
        p(2)*p(2),     -p(1)*p(2), 0; ...
        -p(1)*p(2),    p(1)*p(1),  0; ...
        0,             0,          0];

    dldlT = dldlT./(i2j2*i2j2);

    drd0T2(1:3,1:3) = [ ...
        -2*p(1)*p(1)*p(3),         -2*p(1)*p(2)*p(3),          p(1)*i2j2-p(1)*p(3)*p(3); ...
        -2*p(1)*p(2)*p(3),         -2*p(2)*p(2)*p(3),          p(2)*i2j2-p(2)*p(3)*p(3); ...
        p(1)*i2j2-p(1)*p(3)*p(3),  p(2)*i2j2-p(2)*p(3)*p(3),   2*p(3)*i2j2 ];


    drd0T2 = drd0T2./(rrr*sqrt(i2j2));

    drdlT2(1:3,1:3) = [ ...
        -2*p(1)*p(2),          p(1)*p(1)-p(2)*p(2),    -p(2)*p(3); ...
        p(1)*p(1)-p(2)*p(2),   2*p(1)*p(2),            p(1)*p(3); ...
        -p(2)*p(3),            p(1)*p(3),              0];


    drdlT2 = drdlT2./(r*i2j2);


    dld0T2(1:3,1:3) = [ ...
        2*p(1)*p(2)*p(3),              p(3)*(p(2)*p(2)-p(1)*p(1)), -p(2)*i2j2; ...
        p(3)*(p(2)*p(2)-p(1)*p(1)),    -2*p(1)*p(2)*p(3),          p(1)*i2j2; ...
        -p(2)*i2j2,                    p(1)*i2j2,                  0];


    dld0T2 = dld0T2./(rr*i2j2^(3/2.0));

    % First and Second partial derivatives
    dUdr = 1;   % Takes into account spherical term
    dUd0 = 0;
    dUdl = 0;
    d2Udr2 = 2; % Takes into account spherical term
    d2Ud02 = 0;
    d2Udl2 = 0;
    d2Ud0dr = 0;
    d2Udldr = 0;
    d2Udld0 = 0;

    for l = 2:prpgtr.AuxParam.n
        c01 = (r_ref/r)^l;
        for m = 0:l
            dUdr = dUdr + c01*(l+2)*pnm(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
            dUd0 = dUd0 + c01*(pnm(l+1,m+2)-m*tan(latgc)*pnm(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
            dUdl = dUdl + c01*m*pnm(l+1,m+1)*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));

            d2Udr2 = d2Udr2 + c01*(l+1)*(l+2)*pnm(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
            d2Ud02 = d2Ud02 + c01*(pnm(l+1,m+2)-(2*m+1)*tan(latgc)*pnm(l+1,m+2)+m*(m*tan(latgc)-(1/(cos(latgc)*cos(latgc))))*pnm(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
            d2Udl2 = d2Udl2 + c01*m*m*pnm(l+1,m+1)*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));

            d2Ud0dr = d2Ud0dr + c01*(l+2)*(pnm(l+1,m+2)-m*tan(latgc)*pnm(l+1,m+1))*(C(l+1,m+1)*cos(m*lon)+S(l+1,m+1)*sin(m*lon));
            d2Udldr = d2Udldr + c01*m*(l+2)*pnm(l+1,m+1)*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));
            d2Udld0 = d2Udld0 + c01*m*(pnm(l+1,m+2)-m*tan(latgc)*pnm(l+1,m+1))*(S(l+1,m+1)*cos(m*lon)-C(l+1,m+1)*sin(m*lon));
        end
    end

    dUdr = -dUdr*gm/(rr);
    dUd0 = dUd0*gm/r;
    dUdl = dUdl*gm/r;

    d2Udr2 = d2Udr2*gm/(rrr);
    d2Ud02 = d2Ud02*gm/r;
    d2Udl2 = -d2Udl2*gm/r;

    d2Ud0dr = -d2Ud0dr*gm/(rr);
    d2Udldr = -d2Udldr*gm/(rr);
    d2Udld0 = d2Udld0*gm/r;

    % Acceleration Derivative in TEME Frame


    dadr = d2Udr2*drdrT + dUdr*d2rdr2 + d2Ud02*d0d0T + d2Udl2*dldlT + d2Ud0dr*drd0T2 + d2Udldr*drdlT2 + d2Udld0*dld0T2 + dUd0*d20dr2 + dUdl*d2ldr2;

    % Rotate into body-fixed frame
%         matrixmult(C_i2b,dadr,dadr_b);
    dadr_b = C_i2b*dadr;

    C_b2i = C_i2b';

%         matrixmult(dadr_b,C_b2i,G);
    G = dadr_b*C_b2i;

    % Gravity-Gradient Torque in body-fixed frame
    g_gravity_body(1) = G(2,3)*(Inertia(3,3)-Inertia(2,2)) - G(1,3)*Inertia(1,2) + G(1,2)*Inertia(1,3) + Inertia(2,3)*(G(2,2)-G(3,3));
    g_gravity_body(2) = G(1,3)*(Inertia(1,1)-Inertia(3,3)) + G(2,3)*Inertia(1,2) - G(1,2)*Inertia(2,3) + Inertia(1,3)*(G(3,3)-G(1,1));
    g_gravity_body(3) = G(1,2)*(Inertia(2,2)-Inertia(1,1)) - G(2,3)*Inertia(1,3) + G(1,3)*Inertia(2,3) + Inertia(1,2)*(G(1,1)-G(2,2));
end

end