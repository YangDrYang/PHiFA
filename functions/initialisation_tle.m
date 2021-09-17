function Y0 = initialisation_tle(TLEfile,epoch)

global opsmode
opsmode = 'm'; 
% Initialised by TLE
% TLEfile = 'tle.txt';
fid=fopen(TLEfile,'r');
TLE=fscanf(fid,'%69c%69c');
longstr1=TLE(1:69);
longstr2=TLE(72:140);
clear TLE
fclose(fid);
satrec = twoline2rvMOD(longstr1,longstr2);


if (satrec.epochyr < 57)
    Eyear= satrec.epochyr + 2000;
else
    Eyear= satrec.epochyr + 1900;
end
MjdTLE = Mjday(Eyear, 1, 0, 0, 0, 0);
MjdTLE = MjdTLE + satrec.epochdays;
% 
[~, xsat_eci, vsat_eci] = sgp4(satrec,(epoch-MjdTLE)*1440);
Y0 = [xsat_eci';vsat_eci']*1000;