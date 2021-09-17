function vec_out = j2000toteme(vec_in,MJD_UTC)

p = clPropagator.instance();
[x_pole,y_pole,UT1_UTC,LOD,ddpsi,ddeps,dx_pole,dy_pole,TAI_UTC] = IERS(p.eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_TT  = MJD_UTC + TT_UTC/86400;
JD_TT = MJD_TT + 2400000.5;
ttt = (JD_TT-(12*60*60.0))/(24*60*60.0*36525.0); 

[year,mon,day,hr,min,sec] = invjday(MJD_UTC+2400000.5);
[ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] ...
                    = convtime ( year, mon, day, hr, min, sec, 0, UT1_UTC, TAI_UTC );

% rGCRFtoTEME(JD_TT, 0, 0)
% 
% 
% arcsec2rad = pi/648000;
% 
% r_MOD_GCRF = rGCRFtoMOD_fk5( JD_TT);
% r_TEME_MOD = rMODtoTEME(T, JD_TT,deps*arcsec2rad,dpsi*arcsec2rad);
% 
% 
% 
% vec_out = r_MOD_GCRF*r_TEME_MOD*vec_in;


[prec,psia, wa, ea, xa] = precess ( ttt, '80' );

[deltapsi, trueeps, meaneps, omega, nut] = nutation  (ttt, ddpsi, ddeps );

% ------------------------ find eqeg ----------------------
% rotate teme through just geometric terms 
eqeg = deltapsi* cos(meaneps);

eqeg = rem (eqeg, 2.0*pi);

eqe(1,1) =  cos(eqeg);
eqe(1,2) =  sin(eqeg);
eqe(1,3) =  0.0;
eqe(2,1) = -sin(eqeg);
eqe(2,2) =  cos(eqeg);
eqe(2,3) =  0.0;
eqe(3,1) =  0.0;
eqe(3,2) =  0.0;
eqe(3,3) =  1.0;

tm = eqe * nut' * prec';

vec_out = tm * vec_in;
