function [U, dU] =  transMatEci2Ecef(MJD_UTC)
persistent last_mjd Umat dUmat
if isempty(last_mjd) || last_mjd~=MJD_UTC
    % global eopdata const
    p = clPropagator.instance();

    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(p.eopdata,MJD_UTC,'l');

%     x_pole = 0;
%     y_pole = 0;
%     UT1_UTC = 0;
%     LOD = 0;

    [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);

    MJD_UT1 = MJD_UTC + UT1_UTC/86400;
    MJD_TT  = MJD_UTC + TT_UTC/86400; 

    % ICRS to ITRS transformation matrix and derivative
    P      = PrecMatrix(p.const.MJD_J2000,MJD_TT);     % IAU 1976 Precession
    N      = NutMatrix(MJD_TT);                      % IAU 1980 Nutation
    Theta  = GHAMatrix(MJD_UT1);                     % Earth rotation
    Pi     = PoleMatrix(x_pole,y_pole);              % Polar motion

    S = zeros(3);
    S(1,2) = 1; S(2,1) = -1;                         % Derivative of Earth rotation 
%     Omega = 7292115.8553e-11+4.3e-15*( (MJD_UTC-p.const.MJD_J2000)/36525 ); % [rad/s]
    Omega = p.const.omega_Earth-0.843994809*1e-9*LOD;  % IERS
    dTheta = Omega*S*Theta; % matrix [1/s]

%     Umat      = Pi*Theta; 
%     dUmat     = Pi*dTheta;  
    Umat      = Pi*Theta*N*P;                           % ICRS to ITRS transformation
    dUmat     = Pi*dTheta*N*P;                          % Derivative [1/s]
    
    last_mjd=MJD_UTC;
end

U = Umat;
dU = dUmat;

end