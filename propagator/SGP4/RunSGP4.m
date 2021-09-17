%
%Set paths
clear
%set in sgp4init
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  

global opsmode

dtr=pi/180;
min_per_day=60*24;
path2TLE='.\DemoData';
TLEfile = 'TIROS10.TLE';
fid=fopen(TLEfile,'r');
TLE=fscanf(fid,'%69c%69c');
longstr1=TLE(1:69);
longstr2=TLE(72:140);
clear TLE
fclose(fid);

% %***************Get TLE Orbital Elements**********************
% if 1
%     TLEfiles=dir([path2TLE,'\*.TLE']);
%     if isempty(TLEfiles)
%         error('No TLE file')
%     else
%         fid=fopen(TLEfiles(1).name,'r');
%         TLE=fscanf(fid,'%69c%69c');
%         longstr1=TLE(1:69);
%         longstr2=TLE(72:140);
%         clear TLE
%         fclose(fid);
%     end
% else
%     load([path2TLE,'\SCION-TLE.mat']);
%     longstr1=line1;
%     longstr2=line2;
% end
fprintf('USING 2-Line Elements: \n')
fprintf('%s \n',longstr1)
fprintf('%s \n',longstr2)

%********************Initialize spg4****************************
satrec = twoline2rvMOD(longstr1,longstr2);

fprintf('\n')
fprintf('Satellite ID %5i \n',satrec.satnum)

if (satrec.epochyr < 57)
    Eyear= satrec.epochyr + 2000;
else
    Eyear= satrec.epochyr + 1900;
end
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays);

tsince = 56407.3656089243 - Mjday(Eyear,Emon,Eday,Ehr,Emin,Esec);
[satrec, xsat_eci, vsat_eci] = sgp4(satrec,tsince*86400);
save('IniOrb_TLE.mat','xsat_eci','vsat_eci');



    

