%% init
tic;
last = clLaserStation();
last.referencePower = 10e3;
last.Aperture = 1;
last.FWHM = 1;
last.beamType = eLaserBeam.Gaussian;
last.tumod = @tumodConst;
% last.tumod = @tumodSLCDay;
last.lla(3)=0;
last.xv = [0 6371*sin(pi/4) 6371*cos(pi/4) 0 0 0];
last.init();

deltan=0.005;
D2 = 2;
h0 = 0;
za = 60*pi/180;
targetDist = 160e3;
%%

[aborting,xn,yn,Uturb,propstruct] = ...
        propBeam(last,targetDist,D2,last.tumod,0.005,deltan,za,h0);

N = size(Uturb,1);
mask_circ = circ(xn/D2, yn/D2, 1);
mask_ones = ones(N);
MCF2 = zeros(N);

for i = 1:5
    %% prop beam
    [aborting,xn,yn,Uturb,propstruct] = ...
        propBeam(last,targetDist,D2,last.tumod,0.01,deltan,za,h0);
    
    MCF2 = MCF2 + corr2_ft(Uturb, Uturb, mask_circ, deltan);

end
MCDOC2 = abs(MCF2) / (MCF2(N/2+1,N/2+1));



%% eval
figure;
[thetan, rn] = cart2pol(xn, yn);
r0=propstruct.fried;
rn_lin = linspace(0, 2, 100);
% rn_lin = logspace(-5, 1, 100);
mu_k = exp(-0.5*6.88*(rn_lin).^(5/3));
plot(rn_lin, mu_k);
% r0=0.022610;
hold on
plot(rn(:)./r0, MCDOC2(:), '.','MarkerSize',5);
% plot(rn(:)./r0, MCDOC2(:));
% plot(rn_lin, mu_mvk);
xlim([0 1.5]);
pbaspect([2 1 1]);
grid;
ylabel('Coherence Factor [-]');
xlabel('r_n/r_0 [-]');
legend('theory', 'simulated');

cleanfigure();
plot2tikz('cohefac', 0.5);





