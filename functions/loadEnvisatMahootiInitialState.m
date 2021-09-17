function [Y0, mjd_utc, date] = loadEnvisatMahootiInitialState()
    fid = fopen('propagator/InitialState.txt','r');
    tline = fgetl(fid);
    date.year = str2num(tline(1:4));
    date.mon = str2num(tline(6:7));
    date.day = str2num(tline(9:10));
    date.hour = str2num(tline(12:13));
    date.min = str2num(tline(15:16));
    date.sec = str2num(tline(18:23));
    for i = 1:6
        tline = fgetl(fid);
        Y0(i) = str2num(tline);
    end
    fclose(fid);
    
    mjd_utc = Mjday(date.year, date.mon, date.day, date.hour, date.min, date.sec);
    Y0 = ECEF2ECI(mjd_utc, Y0);
    
end