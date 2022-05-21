TLEfiles = 'H:\My_Research\uqlab_0.901.201508251817_ClosedBeta\SGP4\DemoData\test.txt';
fid=fopen(TLEfiles);
TLE=fscanf(fid,'%69c%69c');
longstr1=TLE(1:69);
longstr2=TLE(72:140);
satrec = twoline2rvMOD(longstr1,longstr2);

if (satrec.epochyr < 57)
    Eyear= satrec.epochyr + 2000;
else
    Eyear= satrec.epochyr + 1900;
end
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays);
UTsec=Ehr*3600+Emin*60+Esec;
%timehack=[Eyear,Emon,Eday,Ehr,Emin,Esec];
%GAST=siderealtime(timehack,0);
gst = gstime(satrec.jdsatepoch);
fprintf(' YEAR MO  DAY UTSEC \n')
fprintf('%5i %2i %4i %5.2f gst=%6.4f rad \n',Eyear,Emon,Eday,UTsec,gst);
fprintf('%5i %2i %4i %4i %4i %5.2f \n',Eyear,Emon,Eday,Ehr,Emin,Esec);