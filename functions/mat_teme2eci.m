function mat = mat_teme2eci(year, mon, day, hr, mm, sec)

    p = clPropagator.instance();
    
    timezone = 0;
    mjd_utc = Mjday(year, mon, day, hr, mm, sec);
    [~,~,UT1_UTC,~,ddpsi,ddeps,~,~,TAI_UTC] = IERS(p.eopdata,mjd_utc,'l');
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
            = convtime ( year, mon, day, hr, mm, sec, timezone, UT1_UTC, TAI_UTC );
    
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

    mat = prec * nut * eqe';    