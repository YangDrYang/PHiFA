function [koe,perigee] = rv2coe_ephemerides(eph)
koe = zeros(8, length(eph));
perigee = zeros(2,length(eph));
for i=1:length(eph)
      [~,a,ecc,incl,omega,argp,nu,m,~,~,~ ] = rv2coe(eph(2:4,i)./1000,eph(5:7,i)./1000);
%       mu = 398600.4415;
%       oev = eci2orb_gooding (mu, eph(2:4,i)./1000,eph(5:7,i)./1000);
%       perigeediff(i) = a*(1-ecc)*1000 - oev(1)*(1-oev(2))*1000;
      koe(1,i) = eph(1,i);
      koe(2,i) = ecc;
      koe(3,i) = a*1000;
      koe(4,i) = incl;
      koe(5,i) = omega;
      koe(6,i) = argp;
      koe(7,i) = nu;
      koe(8,i) = m;
      perigee(1,i) = eph(1,i);
      perigee(2,i) = a*(1-ecc)*1000;
end
end
% oev = eci2orb_gooding (mu, eph(2:4,end)./1000,eph(5:7,end)./1000);
%  outputs       :
%    p           - semilatus rectum               km
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
function [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe(r,v)

    small = 1.0e-10;

    infinite  = 999999.9;
    undefined = 999999.1;

    % -------------------------  mathematical  --------------------
    rad    = 180.0 / pi;
    twopi  = 2.0 * pi;
    halfpi = pi * 0.5;

    % -------------------------  conversions  ---------------------
    ft2m    =    0.3048;
    mile2m  = 1609.344;
    nm2m    = 1852;
    mile2ft = 5280;
    mileph2kmph = 0.44704;
    nmph2kmph   = 0.5144444;
    % -----------------------  physical constants  ----------------
    % WGS-84/EGM-96 constants used here
    re         = 6378.137;         % km
    flat       = 1.0/298.257223563;
    omegaearth = 7.292115e-5;     % rad/s
    mu         = 398600.4418;      % km3/s2
    mum        = 3.986004418e14;   % m3/s2

    % derived constants from the base values
    eccearth = sqrt(2.0*flat - flat^2);
    eccearthsqrd = eccearth^2;

    renm = re / nm2m;
    reft = re * 1000.0 / ft2m;

    tusec = sqrt(re^3/mu);
    tumin = tusec / 60.0;
    tuday = tusec / 86400.0;

    omegaearthradptu  = omegaearth * tusec;
    omegaearthradpmin = omegaearth * 60.0;

    velkmps = sqrt(mu / re);
    velftps = velkmps * 1000.0/ft2m;
    velradpmin = velkmps * 60.0/re;
    %for afspc
    %velkmps1 = velradpmin*6378.135/60.0   7.90537051051763
    %mu1 = velkmps*velkmps*6378.135        3.986003602567418e+005        
    degpsec = (180.0 / pi) / tusec;
    radpday = 2.0 * pi * 1.002737909350795;

    speedoflight = 2.99792458e8; % m/s
    au = 149597870.0;      % km
    earth2moon = 384400.0; % km
    moonradius =   1738.0; % km
    sunradius  = 696000.0; % km

    masssun   = 1.9891e30;
    massearth = 5.9742e24;
    massmoon  = 7.3483e22;
    muin = mu; % this is the km version
    
    % -------------------------  implementation   -----------------
    magr = mag( r );
    magv = mag( v );
    % ------------------  find h n and e vectors   ----------------
    [hbar] = cross( r,v );
    magh = mag( hbar );
    if ( magh > small )
        nbar(1)= -hbar(2);
        nbar(2)=  hbar(1);
        nbar(3)=   0.0;
        magn = mag( nbar );
        c1 = magv*magv - muin /magr;
        rdotv= dot( r,v );
        for i= 1 : 3
            ebar(i)= (c1*r(i) - rdotv*v(i))/muin;
        end
        ecc = mag( ebar );

        % ------------  find a e and semi-latus rectum   ----------
        sme= ( magv*magv*0.5  ) - ( muin /magr );
        if ( abs( sme ) > small )
            a= -muin  / (2.0 *sme);
          else
            a= infinite;
        end
        p = magh*magh/muin;

        % -----------------  find inclination   -------------------
        hk= hbar(3)/magh;
        incl= acos( hk );

        % --------  determine type of orbit for later use  --------
        % ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit= 'ei';
        if ( ecc < small )
            % ----------------  circular equatorial ---------------
            if  (incl<small) || (abs(incl-pi)<small)
                typeorbit= 'ce';
              else
                % --------------  circular inclined ---------------
                typeorbit= 'ci';
            end
          else
            % - elliptical, parabolic, hyperbolic equatorial --
            if  (incl<small) || (abs(incl-pi)<small)
                typeorbit= 'ee';
            end
        end

        % ----------  find longitude of ascending node ------------
        if ( magn > small )
            temp= nbar(1) / magn;
            if ( abs(temp) > 1.0  )
                temp= sign(temp);
              end
            omega= acos( temp );
            if ( nbar(2) < 0.0  )
                omega= twopi - omega;
            end
          else
            omega= undefined;
        end

        % ---------------- find argument of perigee ---------------
        if (strcmp(typeorbit, 'ei') == 1)
            argp = angl( nbar,ebar);
            if ( ebar(3) < 0.0  )
                argp= twopi - argp;
            end
          else
            argp= undefined;
        end

        % ------------  find true anomaly at epoch    -------------
        if ( typeorbit(1:1) == 'e' )
            nu =  angl( ebar,r);
            if ( rdotv < 0.0  )
                nu= twopi - nu;
            end
          else
            nu= undefined;
        end

        % ----  find argument of latitude - circular inclined -----
        if (strcmp(typeorbit, 'ci') == 1)
            arglat = angl( nbar,r );
            if ( r(3) < 0.0  )
                arglat= twopi - arglat;
            end
            m = arglat;
          else
            arglat= undefined;
        end

        % -- find longitude of perigee - elliptical equatorial ----
        if  ( ecc>small ) && (strcmp(typeorbit, 'ee') == 1)
            temp= ebar(1)/ecc;
            if ( abs(temp) > 1.0  )
                temp= sign(temp);
            end
            lonper= acos( temp );
            if ( ebar(2) < 0.0  )
                lonper= twopi - lonper;
            end
            if ( incl > halfpi )
                lonper= twopi - lonper;
            end
          else
            lonper= undefined;
        end

        % -------- find true longitude - circular equatorial ------
        if  ( magr>small ) && (strcmp(typeorbit, 'ce') == 1)
            temp= r(1)/magr;
            if ( abs(temp) > 1.0  )
                temp= sign(temp);
            end
            truelon= acos( temp );
            if ( r(2) < 0.0  )
                truelon= twopi - truelon;
            end
            if ( incl > halfpi )
                truelon= twopi - truelon;
            end
            m = truelon;
          else
            truelon= undefined;
        end

        % ------------ find mean anomaly for all orbits -----------
        if ( typeorbit(1:1) == 'e' )
            [e,m] = newtonnu(ecc,nu );
        end

     else
       p    = undefined;
       a    = undefined;
       ecc  = undefined;
       incl = undefined;
       omega= undefined;
       argp = undefined;
       nu   = undefined;
       m    = undefined;
       arglat = undefined;
       truelon= undefined;
       lonper = undefined;
    end
end

function mag = mag ( vec )

    temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);

    if abs( temp ) >= 1.0e-16
        mag= sqrt( temp );
    else
        mag= 0.0;
    end
end

function [theta] = angl ( vec1,vec2 )

    small     = 0.00000001;
    undefined = 999999.1;

    magv1 = mag(vec1);
    magv2 = mag(vec2);

    if magv1*magv2 > small^2
        temp= dot(vec1,vec2) / (magv1*magv2);
        if abs( temp ) > 1.0
            temp= sign(temp) * 1.0;
        end
        theta= acos( temp );
    else
        theta= undefined;
    end
end

function [e0,m] = newtonnu ( ecc,nu );

    % ---------------------  implementation   ---------------------
    e0= 999999.9;
    m = 999999.9;
    small = 0.00000001;

    % --------------------------- circular ------------------------
    if ( abs( ecc ) < small  )
        m = nu;
        e0= nu;
    else
        % ---------------------- elliptical -----------------------
        if ( ecc < 1.0-small  )
            sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
            cose= ( ecc + cos(nu) ) / ( 1.0  + ecc*cos(nu) );
            e0  = atan2( sine,cose );
            m   = e0 - ecc*sin(e0);
        else
            % -------------------- hyperbolic  --------------------
            if ( ecc > 1.0 + small  )
                if (ecc > 1.0 ) && (abs(nu)+0.00001 < pi-acos(1.0 /ecc))
                    sine= ( sqrt( ecc*ecc-1.0  ) * sin(nu) ) / ( 1.0  + ecc*cos(nu) );
                    e0  = asinh( sine );
                    m   = ecc*sinh(e0) - e0;
                end
            else
                % ----------------- parabolic ---------------------
                if ( abs(nu) < 168.0*pi/180.0  )
                    e0= tan( nu*0.5  );
                    m = e0 + (e0*e0*e0)/3.0;
                end
            end
        end
    end

    if ( ecc < 1.0  )
        m = rem( m,2.0 *pi );
        if ( m < 0.0  )
            m= m + 2.0 *pi;
        end
        e0 = rem( e0,2.0 *pi );
    end
end