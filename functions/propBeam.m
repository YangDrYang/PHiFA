function [aborting,xn,yn,Uout,varargout] = propBeam(last,Dz,D2,tumod,delta1max,deltanmax,za,h0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
start=toc;
wvl=last.Wavelength;
D1=1.5*last.Aperture;
aborting=false;
epsilon=8.8542E-12;
c_light=299792458.000000;
%% INIT
k = 2*pi / wvl; % optical wavenumber [rad/m]
% bias and fluctuation of radius
R = ( 1 + last.focusBias/100) * Dz * ( 1 - last.focusFluctuation*randn/100);

Cn2_r0sw = @(x) tumod(x,za,h0)*(x/Dz)^(5/3);
intCn2_r0sw = integral(Cn2_r0sw,0,Dz,'ArrayValued',true);
% atmospheric properties
r0sw = (0.423 * k^2 * intCn2_r0sw)^(-3/5);
fprintf('FRIED parameter: %f cm\n', r0sw*100)
% log-amplitude variance
Cn2_rytov = @(x) tumod(x,za,h0)*x^(5/6)*(1-x/Dz)^(5/6);
intCn2_rytov = integral(Cn2_rytov,0,Dz,'ArrayValued',true);
rytov = 0.563 * k^(7/6) * intCn2_rytov;
fprintf('RYTOV Log-amplitude variance: %f\n', rytov);
% if rytov > 0.25 % theory of weak fluctuations doenst apply for that strong turbulence
%     fprintf('Log-amplitude variance high. Aborting pulse.\n');
%     aborting=true;
%     if nargout>4
%         varargout = cell(1,1);
%         varargout{1} = struct('fried',r0sw, 'rytov',rytov, 'maxint',0,...
%         'turberror',0, 'nscr',0);
%     end
%     xn=0; yn=0; Uout=0;
%     return;
% end

%% CONSTRAINTS
[d1,d2,N,nmin] = fulfilConstrains(D1,D2,wvl,Dz,r0sw,R,delta1max,deltanmax);
fprintf('Found Constrains: d1 = %f; d2 = %f; N = %i; nscr = %i\n',d1,d2,N,nmin);
nscr = 11; % number of screens

%% VACUUM PROPAGATION
delta1 = d1;    % source-plane grid spacing [m]
deltan = d2;    % observation-plane grid spacing [m]
n = nscr;         % number of planes

% coordinates
[x1, y1] = meshgrid((-N/2 : N/2-1) * delta1);
[theta1, r1] = cart2pol(x1, y1);

transmittedPower = last.referencePower * last.telescopeEff * ...
    (1 - last.powerFluctuation*randn(1)/100);
% w0 = last.FWHM * ( 2 * sqrt( 2 * log( 2 ) ) );
w0 = last.FWHM / ( 2 * sqrt( log( 2 ) ) );
% leaving beam distribution
Uin = sqrt( 4*transmittedPower/(c_light*epsilon*pi*w0^2) ) .* ...
    exp(-1i*k/(2*R) * r1.^2) ...
    .* exp(-(r1.^2/(w0^2)));
% partial prop planes
z = (1 : n-1) * Dz / (n-1);

sg = exp(-(x1/(0.47*N*d1)).^16) ...
        .* exp(-(y1/(0.47*N*d1)).^16);

% theory doesnt apply for this situation so we just do as if we propagate through vacuum
if r0sw>last.Aperture*100
    % simulate vacuum propagation 
    t = repmat(sg, [1 1 n]);    
    fprintf('Simulate vacuum propagation\n');
    [xn, yn, Uvac] = ang_spec_multi_prop(Uin, wvl, ...
        delta1, deltan, z, t);
    % collimate the beam
    
    Uout = Uvac .* exp(-1i*pi/(wvl*R)*(xn.^2+yn.^2));
else
    %% screen properties
    A = zeros(2, nscr); % matrix
    alpha = (0:nscr-1) / (nscr-1);
    A(1,:) = alpha.^(5/3);
    A(2,:) = (1 - alpha).^(5/6) .* alpha.^(5/6);
    b = [r0sw.^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
    % initial guess
    x0 = (nscr/3*r0sw * ones(nscr, 1)).^(-5/3);
    % objective function
    fun = @(X) sum((A*X(:) - b).^2);
    % constraints
    x1 = zeros(nscr, 1);
    rmax = 0.1; % maximum Rytov number per partial prop
    x2 = rmax/1.33*(k/Dz)^(5/6) ./ A(2,:);
    x2(A(2,:)==0) = 50^(-5/3);
    options = optimoptions('fmincon','Display','off');
    [X,fval,exitflag,output] ...
        = fmincon(fun,x0,[],[],[],[],x1,x2,[],options);
    % check screen r0s
    r0scrn = X.^(-3/5);
    r0scrn(isinf(r0scrn)) = 1e6;

    %% TURBULENCE PROPAGATION
    l0 = 0;     % inner scale [m]
    L0 = inf;   % outer scale [m]

    zt = [0 z];  % propagation plane locations
    Delta_z = zt(2:n) - zt(1:n-1);    % propagation distances
    % grid spacings
    alpha = zt / zt(n);
    delta = (1-alpha) * delta1 + alpha * deltan;

    % initialize array for phase screens
    phz = zeros(N, N, n);
%     Uturb = zeros(N);
    sg = repmat(sg, [1 1 n]);
    % loop over screens
    fprintf('Starting evaluating phase screens\n');
    parfor idxscr = 1 : 1 : n
        [phz_lo, phz_hi] ...
            = ft_sh_phase_screen ...
            (r0scrn(idxscr), N, delta(idxscr), L0, l0);
        phz(:,:,idxscr) = phz_lo + phz_hi;
    end
    % simulate turbulent propagation
    fprintf('Simulate turbulent propagation\n');
    [xn, yn, Uturb] = ang_spec_multi_prop(Uin, wvl, ....
        delta1, deltan, z, sg.*exp(1i*phz));
    % collimate the beam
    Uout = Uturb .* exp(-1i*pi/(wvl*R)*(xn.^2+yn.^2));
end
Uout = abs(Uout); % since its still complex after turbulence propagation
% fprintf('Maximum Intensity before propagation: %f W/m^2\n',max(abs(Uin(:))));
% fprintf('Maximum Intensity after propagation: %f W/m^2\n',max(abs(Uout(:))));
fprintf('Turbulence propagation done. Time taken: %f sec\n',toc-start);

if nargout > 4
    [intmax,poserror] = findMaxIntensity(Uout,xn,yn);
    varargout = cell(1,1);
    varargout{1} = struct('fried',r0sw, 'rytov',rytov, 'maxint',intmax,...
        'turberror',poserror/Dz, 'nscr',nscr, 'N', N);
end
end

