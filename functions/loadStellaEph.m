function [eph, mjd0] = loadStellaEph(loc)

GPS_UTC = 14; % from first Jan 2006 to first Jan 2009
UTC_GPS = - GPS_UTC;

% Open file
fileID = fopen(loc,'r');

s=0;
% Read file
while (~feof(fileID))
    temp = fgetl(fileID);
    if ~isempty(temp)
        if strcmp(temp(1),'*')
            s = s+1;
        end
    end
end

% Close file
fclose(fileID);

eph = zeros(7,s);

% Open file
fileID = fopen(loc,'r');

% Read file
j=1;
while (~feof(fileID))
    temp = fgetl(fileID);
    if ~isempty(temp)
        if strcmp(temp(1),'*')
            date = sscanf(temp, '%s\t%d\t%d\t%d\t%d\t%d\t%f');
            mjd = Mjday(date(2), date(3), date(4), date(5), date(6), date(7)+UTC_GPS);
            temp = fgetl(fileID);
            pos = sscanf(temp, '%s\t%f\t%f\t%f\t%f');
            
            temp = fgetl(fileID);
            velo = sscanf(temp, '%s\t%f\t%f\t%f\t%f');
            
            eph(:,j) = [mjd; pos(5:7).*1000;velo(5:7)./10];
            j=j+1;
        end
    end
end
% Close file
fclose(fileID);

mjd0 = eph(1,1);
for i = 1:length(eph)
    eph(1,i) = (eph(1,i) - mjd0)*86400;
end

end