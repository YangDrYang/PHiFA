%% init
tic;
last = clLaserStation();
last.referencePower = 10e3;
last.Aperture = 1;
last.FWHM = 0.5;
last.beamType = eLaserBeam.Gaussian;
last.tumod = @tumodSLCDay;
last.lla(3)=0;
last.xv = [0 6371*sin(pi/4) 6371*cos(pi/4) 0 0 0];
last.init();

D2 = 2;

za_vec = [0 30 60];
% figure('Name', ['propBeam' num2str(i)]);
for i = 1:length(za_vec)
    %% prop beam
    targetDist = 700000;
    za = za_vec(i)*pi/180;
    h0 = 0;

    disp(num2str(za_vec(i)));
    [aborting,xn,yn,Uturb,propstruct] = ...
                        propBeam(last,targetDist,D2,last.tumod,0.05,0.01,za,h0);

    %% plotting distrib
%     close all;
    za_deg = za*180/pi;
    magnpt = (8.8542187817E-12*299792458/2).*(abs(Uturb)).^2./1000;
%     subplot(2,2,i);
    figure;
    contour3(xn,yn,magnpt,75);
    surf(xn,yn,magnpt);
    shading interp;
    h = colorbar;
%     title(['zenith angle ' num2str(za_deg) ' deg']);
    xlabel('x [m]');
    ylabel('y [m]');
    ylabel(h, 'Intensity [W/cm^2]');
    xlim([-2 2]);
    ylim([-2 2]);
    zlim([0 100]);
    view(2);
    plot2tikz(sprintf('bp_%dza', za_vec(i)),0.6);
end

%% uin
figure;
wvl=last.Wavelength;
k = 2*pi / wvl; % optical wavenumber [rad/m]
D1=1.5*last.Aperture;
epsilon=8.8542E-12;
c_light=299792458.000000;
delta1=0.05;
N=1024;
Dz=targetDist;
[x1, y1] = meshgrid((-N/2 : N/2-1) * delta1);
[theta1, r1] = cart2pol(x1, y1);
R = ( 1 + last.focusBias/100) * Dz * ( 1 - last.focusFluctuation*randn/100);
transmittedPower = last.referencePower * last.telescopeEff * ...
    (1 - last.powerFluctuation*randn(1)/100);
w0 = last.FWHM * ( 2 * sqrt( 2 * log( 2 ) ) );
% leaving beam distribution
Uin = sqrt( 4*transmittedPower/(c_light*epsilon*pi*w0^2) ) .* ...
    exp(-1i*k/(2*R) * r1.^2) ...
    .* exp(-(r1.^2/(w0^2)));
magnpt = (8.8542187817E-12*299792458/2).*(abs(Uin)).^2./1000;
surf(x1,y1,magnpt);
shading interp;
h = colorbar;
%     title(['zenith angle ' num2str(za_deg) ' deg']);
xlabel('x [m]');
ylabel('y [m]');
ylabel(h, 'Intensity [W/cm^2]');
xlim([-2 2]);
ylim([-2 2]);
zlim([0 100]);
view(2);
plot2tikz('bp_uin',0.6);