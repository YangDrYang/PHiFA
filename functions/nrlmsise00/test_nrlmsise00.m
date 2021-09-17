%% several nrlmsise00 test profiles
% 30/01/2016
% clc
% clear
close all;
format longe

path = fileparts(mfilename('fullpath'));
cd(path);

iyr = 2003;
iday = [172, 81, 172*ones(1,13)]';
ut = [29000, 29000, 75000, 29000*ones(1,12)]';
alt = [400,400,1000,100,400*ones(1,6), 0, 10,30,50,70]';
xlat = [60*ones(1,4), 0, 60*ones(1,10)]';
xlong = [-70*ones(1,5), 0,-70*ones(1,9)]';
% xlst = [16.0*ones(1,6),4,16.0*ones(1,8)]';
xlst = ut/3600 + xlong/15;
f107a = [150.0*ones(1,7),70,150.0*ones(1,7)]';
f107 = [150.0*ones(1,8),180,150.0*ones(1,6)]';
ap = [4.0*ones(1,9),40,4.0*ones(1,5)]'*ones(1,7);
aph = ones(15,7)*100;

flags = ones(1,23);
flags(9)=-1;
%% mex file
% i = 14;
% [d, t] = nrlmsise00(iday(i),ut(i),alt(i),xlat(i),xlong(i), f107a(i), f107(i), ap(i,:), flags);
% disp([d(6),t]);
%% stop
tic
if flags(9)==-1
  [d, t] = nrlmsise00(iday,ut,alt,xlat,xlong, ...
    f107a, f107, aph, flags,'NoOxygen');
else
  [d, t] = nrlmsise00(iday,ut,alt,xlat,xlong, ...
    f107a, f107, ap, flags,'NoOxygen');
end
toc
disp([d(1:15,:)]);
fprintf('\n');
fprintf('%8.3f  %8.3f\n',t(1,:)');
fprintf('\n');
%% built-in function
tic
if flags(9)==-1
  [T rho] = atmosnrlmsise00(alt*1000, xlat, xlong, 0, iday, ut, ...
    f107a, f107, aph, flags,'NoOxygen');
else
  [T rho] = atmosnrlmsise00(alt*1000, xlat, xlong, 0, iday, ut, ...
    f107a, f107, ap, flags,'NoOxygen');
end
toc
rho(1:15,6)
max(abs(d-rho),[],1)
return
%% test O+
ap=[6.375,7,7,12,7,4.875,3.111];
[d, t] = nrlmsise00(2009,247,43200,...
  500,39.03,117.68, ...
  69, 68.2, ap,flags(1:25),'NoOxygen');
d(:,6)
[T rho] = atmosnrlmsise00(500*1000, 39.03, 117.68, 2009, 247, 43200,...
  69, 68.2, ap, flags(1:23),'NoOxygen');
rho(:,6)
% with O+
[d, t] = nrlmsise00(2009,247,43200,...
  250,39.03,117.68, ...
  69, 68.2, ap,ones(1,25),'Oxygen');
d(:,6)

[T rho] = atmosnrlmsise00(250*1000, 39.03, 117.68, 2009, 247, 43200,...
  69, 68.2, ap, ones(1,23),'Oxygen');
rho(:,6)

