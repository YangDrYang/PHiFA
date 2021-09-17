function [C, S] = loadMahootiGravCoef()
C = zeros(181,181);
S = zeros(181,181);
% fid = fopen('propagator/GGM03S.txt','r');
fid = fopen('propagator/GGM05S_ICGEM.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        C(n+1,m+1) = temp(3);
        S(n+1,m+1) = temp(4);
    end
end
fclose(fid);
end