%% test mex function gravity_field

numout = length(0:50:2000)*length(0:2*pi/20:2*pi)*length(-pi/2:pi/10:pi/2);
out = zeros(numout,9);
normp = zeros(numout,1);
norma = zeros(numout,1);
normg = zeros(numout,1);
phiout = zeros(numout,1);
latout = zeros(numout,1);
lonout = zeros(numout,1);
hout = zeros(numout,1);
nout = zeros(numout,1);
betaout = zeros(numout,1);

C_i2b = [0 0 1; 0 1 0; 1 0 0];
Inertia = rand(3);
[C, S] = loadDSPOSEGravCoef();
l_max_a = 100;
l_max_g = 100;
n_grav_a = 1;
n_grav_g = 1;
    
n = 1;
% dispstat('','init');
for height = 0:50:2000
    for phi = 0:2*pi/20:2*pi
        for beta = -pi/2:pi/10:pi/2
    %         dispstat(sprintf('%d%%',pot/maxpow*100));
            p = ((6738-1000) + height)*1000*[cos(phi) sin(phi) cos(beta)];
            p_ecef = ECI2ECEF(52388.9135185187, [p 0 0 0]);
        %     [latgc,latgd,lon,hellp] = ijk2ll ( p_ecef );
            wgs84 = wgs84Ellipsoid('meters');
            [lat,lon,h] = ecef2geodetic(wgs84,p_ecef(1),p_ecef(2),p_ecef(3), 'radians');
            LLA = [geocentricLatitude(lat, wgs84.Flattening, 'radians'),lat,lon,h/1000];

            [a_gravity_inertial, g_gravity_body] = ...
                mexGravityField(p_ecef,p,LLA,C_i2b,Inertia,C,S,l_max_a,l_max_g,n_grav_a,n_grav_g);

            out(n,1:3) = p;
            normp(n) = norm(p_ecef);
            hout(n) = h/1000;
            latout(n) = lat;
            lonout(n) = lon;
            out(n,4:6) = a_gravity_inertial;
            norma(n) = norm(a_gravity_inertial);
            out(n,7:9) = g_gravity_body;
            normg(n) = norm(g_gravity_body);
            phiout(n) = phi;
            nout(n) = n;
            betaout(n) = beta;
            n=n+1;
        end
    end
end
figure;
semilogy(nout, norma,'.', nout, normg,'x');
legend('norm acceleration','norm torque');
figure;
plot(nout, phiout,'.', nout, hout,'x', nout, betaout, '+');
legend('norm acceleration','norm torque');

% figure;
% semilogy(out(:,1), norma,'.', out(:,2), norma,'x', out(:,3), norma,'+');
% legend('norm acceleration','norm torque');
% figure;
% semilogy(out(:,4), norma,'.', out(:,5), norma,'x', out(:,6), norma,'+');
% legend('norm acceleration','norm torque');
% figure;
% semilogy(out(:,7), norma,'.', out(:,8), norma,'x', out(:,9), norma,'+');
% legend('norm acceleration','norm torque');