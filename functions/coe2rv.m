%
% ------------------------------------------------------------------------------
%
%                           function coe2rv
%
%  this function finds the position and velocity vectors in geocentric
%    equatorial (ijk) system given the classical orbit elements.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%    vallado     - add constant file use                         29 jun 2003
%
%  inputs          description                    range / units
%    p           - semilatus rectum               km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omega       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
%    truelon     - true longitude            (ce) 0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
%
%  outputs       :
%    r           - ijk position vector            km
%    v           - ijk velocity vector            km / s
%
%  locals        :
%    temp        - temporary real*8 value
%    rpqw        - pqw position vector            km
%    vpqw        - pqw velocity vector            km / s
%    sinnu       - sine of nu
%    cosnu       - cosine of nu
%    tempvec     - pqw velocity vector
%
%  coupling      :
%    mag         - magnitude of a vector
%    rot3        - rotation about the 3rd axis
%    rot1        - rotation about the 1st axis
%
%  references    :
%    vallado       2007, 126, alg 10, ex 2-5
%
% [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper );
% ------------------------------------------------------------------------------

function [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper )

        %% -------------------------  implementation   -----------------
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

        %% -------------------------------------------------------------
        %       determine what type of orbit is involved and set up the
        %       set up angles for the special cases.
        % -------------------------------------------------------------
        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) || ( abs(incl-pi)< small )
                argp = 0.0;
                omega= 0.0;
                nu   = truelon;
            else
                % --------------  circular inclined  ------------------
                argp= 0.0;
                nu  = arglat;
            end
        else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl<small) || (abs(incl-pi)<small) )
                argp = lonper;
                omega= 0.0;
            end
        end

        % ----------  form pqw position and velocity vectors ----------
        cosnu= cos(nu);
        sinnu= sin(nu);
        temp = p / (1.0  + ecc*cosnu);
        rpqw(1)= temp*cosnu;
        rpqw(2)= temp*sinnu;
        rpqw(3)=     0.0;
        if ( abs(p) < 0.0001)
            p= 0.0001;
        end
        vpqw(1)=    -sinnu*sqrt(mu)  / sqrt(p);
        vpqw(2)=  (ecc + cosnu)*sqrt(mu) / sqrt(p);
        vpqw(3)=      0.0;

        % ----------------  perform transformation to ijk  ------------
        [tempvec] = rot3( rpqw   , -argp );
        [tempvec] = rot1( tempvec, -incl );
        [r] = rot3( tempvec, -omega );

        [tempvec] =rot3( vpqw   , -argp );
        [tempvec] =rot1( tempvec, -incl );
        [v] = rot3( tempvec, -omega );

        r=r';
        v=v';
end

% ------------------------------------------------------------------------------
%
%                            function mag
%
%  this function finds the magnitude of a vector.  the tolerance is set to
%    0.000001, thus the 1.0e-12 for the squared test of underflows.
%
%  author        : david vallado                  719-573-2600   30 may 2002
%
%  revisions
%    vallado     - fix tolerance to match coe, eq, etc            3 sep 2002
%
%  inputs          description                    range / units
%    vec         - vector
%
%  outputs       :
%    mag         - magnitude
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
% mag = ( vec );
% ----------------------------------------------------------------------------- }

function mag = mag ( vec );

        temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);

        if abs( temp ) >= 1.0e-16
            mag= sqrt( temp );
          else
            mag= 0.0;
          end

end

% ------------------------------------------------------------------------------
%
%                                  rot1
%
%  this function performs a rotation about the 1st axis.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    vec         - input vector
%    xval        - angle of rotation              rad
%
%  outputs       :
%    outvec      - vector result
%
%  locals        :
%    c           - cosine of the angle xval
%    s           - sine of the angle xval
%    temp        - temporary extended value
%
%  coupling      :
%    none.
%
% [outvec] = rot1 ( vec, xval );
% ----------------------------------------------------------------------------- }

function [outvec] = rot1 ( vec, xval );

        temp= vec(3);
        c= cos( xval );
        s= sin( xval );

        outvec(3)= c*vec(3) - s*vec(2);
        outvec(2)= c*vec(2) + s*temp;
        outvec(1)= vec(1);


end

% ------------------------------------------------------------------------------
%
%                            function rot3
%
%  this function performs a rotation about the 3rd axis.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    vec         - input vector
%    xval        - angle of rotation              rad
%
%  outputs       :
%    outvec      - vector result
%
%  locals        :
%    c           - cosine of the angle xval
%    s           - sine of the angle xval
%    temp        - temporary extended value
%
%  coupling      :
%    none.
%
% [outvec] = rot3 ( vec, xval );
% ----------------------------------------------------------------------------- }

function [outvec] = rot3 ( vec, xval );

        temp= vec(2);
        c= cos( xval );
        s= sin( xval );

        outvec(2)= c*vec(2) - s*vec(1);
        outvec(1)= c*vec(1) + s*temp;
        outvec(3)= vec(3);

end
