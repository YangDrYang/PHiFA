function [d1,d2,N,nmin] = fuFulfilPropagationConstrains(D1,D2,wvl,Dz,r0sw,R,d1max,d2max,Nmin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Npowmax = 12;
c = 3;
D1p = D1 + c*wvl*Dz/r0sw;
D2p = D2 + c*wvl*Dz/r0sw;

% constraint 1
deltan_max = @(delta1) -D2p/D1p*delta1 + wvl*Dz/D1p;
% constraint 3 + 1
d1 = wvl*Dz/D1p/(1+Dz/R+D2p/D1p); % built-in contraint 3
d2 = deltan_max(d1);
d1 = d1*0.9; d2 = d2*0.9; % win some distance to limit
% consider extern limits
d1 = min([d1 d1max]);
d2 = min([d2 d2max]);
% constraint 2
Nmin2 = (wvl * Dz + D1p*d2 + D2p*d1) ...
    ./ (2 * d1 .* d2);
% consider extern limits
Npow = nextpow2(max([Nmin2 Nmin])); % fft more efficient when N = 2^x
N = 2^Npow;
% constraint 4
mindsqr = min([d1 d2])^2;
zmax = mindsqr * N / wvl;
nmin = ceil(Dz / zmax) + 1;
% optimice N
while Npow>Npowmax
    nmin=nmin+1;
    dz=Dz/nmin;
    Npow = nextpow2(dz*wvl/mindsqr);
end
N = 2^Npow;
end

