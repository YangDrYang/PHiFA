%{
This code performs pair-wise comparisons of Two-line Element sets (TLEs)
by propagating the TLEs using SGP4 within Satellite ToolKit.
Running this code requires MATLAB Aerospace Toolbox and STK 7.0. See
help files for "atb" functions used in this code.
This code is set up for evaluating 6 different satellites during 8 separate
time windows.
Three primary functions are performed by this code:
1) Position residuals resulting from the pair-wise comparisons are plotted.
2) Covariance matrices for the state vector that is computed for the most recent TLE of each time window.
3) An autocorrelation function is formulated by extending the time windows
of propagation and then combining all of them together.
Reference thesis document mentioned below for more information.
author : Capt Victor Osweiler 23 March 2006
In support of Masters Degree from Air Force Institute of Technology (AFIT).
Thesis Designator: AFIT/GSS/ENY/06-M09.
Thesis Title:
"COVARIANCE ESTIMATION AND AUTOCORRELATION OF NORAD TWO-LINE ELEMENT SETS"
inputs :
This script opens, reads in, and uses TLEs in their standard defined
format. TLEs are saved in one file for each satellite, in a specific
location on the available network (or hard drive).
The script can be modified to read TLEs from an alternate location.
outputs :
Plots -
Excel Files -
Covariance Matrix - 6x6 satellite-based coordinate system covariance
Covariance Matrix - 6x6 classical orbital elements covariance matrix
Autocorrelation Function -
Sub-functions called in the code :
Epoch2Date from Vallado
TwoDigit from Vallado
rv2vnc
rv2rtc
mag from Vallado
COVCT2CL from Vallado
%}

% --------------------------------------------------------------------
clear all; %#ok<CLSCR>
clc;
warning off MATLAB:xlswrite:AddSheet;
agiInit;
% -------------------------------------------------------
% Options in script: choose whether to plot, write excel
% files, do autocorrelation function and plot, as well as the
% satellite numbers to evaluate and the number of time windows to use
% NOTE: see below for time period windows
toplotornottoplot = true;
logfit = false;
writeexcel = true;
autocorrelate = false;
% ---------------------------------------------------------
% satellite # to work with (from 1 to 6)
% NOTE: see below for satellite catalog numbers
for satnr = 1:6
% time window # to work with (from 1 to 8)
% NOTE: see below for TLE epoch dates for the 8 windows
for timeloop = 1:8
% set STK internal clock to default value
atbSetEpoch(1972, 1, 1, 0, 0, 0, 'Z', 0);
t = timeloop;
dout = ['I:\My Documents\Thesis\Output\output_',...
num2str(satnr),num2str(t),'.txt'];
diary(dout)
satnr, t
temp1 = num2str(satnr);
temp2 = num2str(timeloop);
temp3 = [temp1,temp2];
temp4 = str2num(temp3);
fignr = temp4;
format long g;
TLEsize = 142;
buffer = 30;
% The following dates are the start and end dates in YYDDD.00000000 format
% that I will use for my time windows of propagation
startdates = ['03060.00000000';'03079.00000000';'03270.00000000';
'03290.00000000';'04045.00000000';'04065.00000000';
'04150.00000000';'04280.00000000']
enddates = ['03075.00000000';'03093.00000000';'03285.00000000';
'03305.00000000';'04060.00000000';'04080.00000000';
'04165.00000000';'04295.00000000']

%%%% This are the extended time windows which had to be implemented for
%%%% the autocorrelation function.
% startdates = ['03040.00000000';'03059.00000000';'03250.00000000';
% '03270.00000000';'04025.00000000';'04045.00000000';
% '04130.00000000';'04260.00000000']
% enddates = ['03095.00000000';'03113.00000000';'03305.00000000';
% '03325.00000000';'04080.00000000';'04100.00000000';
% '04185.00000000';'04315.00000000']
for i = 1:8
stimewindow(i) = Epoch2Date(startdates(i,:));
etimewindow(i) = Epoch2Date(enddates(i,:));
end;
timewindow = [stimewindow' etimewindow']
startdate = timewindow(t,1)
enddate = timewindow(t,2)
satids = ['08820';'22076';'25933';'27391';'27642';'24285'];
catnr = satids(satnr,:)
fname = ['sat',catnr,'.2le']
fi = fopen(['I:\My Documents\Thesis\TLEs\',fname],'rt')
% fi = fopen(fname,'rt');
% fi = fopen('I:\My Documents\Thesis\Main MATLAB Code\sat27642.txt','rt');
% fi = fopen('sat27642.2le','rt')
% if (fi = -1) % file does not exist, do what?
if (fi ~= -1) % file exists
fseek(fi,0,'eof');
totalcount = ftell(fi)/TLEsize;
fseek(fi,0,'bof');
line = textscan(fi,'%s%s','delimiter','\n');
anydata = false;
for (n = 1:totalcount)
line1 = char(line{1}(n,1));
edate = Epoch2Date(line1(19:32));
if (edate > startdate)
anydata = true;
break;
end; % if
end; % for n
n;
if (anydata) % if there is a TLE after the start date
count = 0;
n;
for i = n:totalcount
line1 = char(line{1}(i,1));
line1 = [line1(1:1),' ',line1(3:69)];
line2 = char(line{2}(i,1));
line2 = [line2(1:1),' ',line2(3:69)];
edate = Epoch2Date(line1(19:32));
if (edate > enddate)
break;
end; %if past enddate
count = count + 1;
name = ['TLE',TwoDigit(count)];
atbTLESetInfo(name,line1,line2)
end; % for i
m = 0;
count
%%% clear delta PosDiffVNC VelDiffVNC PosDiffRTC PosDiffRTC data;
points = zeros(15,1);
sumVNC = zeros(15,3);
sumVNC2 = zeros(15,3);
sumRTC = zeros(15,3);
sumRTC2 = zeros(15,3);
data = zeros(20,3,count,count);
%%%%%% Begin double loop of propagating each TLE in the time window to
%%%%%% every other TLE which follows it (is more recent in time).
%%%%%% Outer loop: Start with the last TLE in the time window, go backward.
for j = count:-1:1
jName = ['TLE',TwoDigit(j)];
jEpochSec = atbTLEEpoch(jName);
jDate = atbEpochSecToDate(jEpochSec);
atbSetEpoch(jDate);
[priTimes,priPosECF,priVelECF] = atbTLEProp(jName,0.0,0.0,1.0);
[priPosECI, priVelECI] = atbCbfToCbi('Earth',priTimes,priPosECF,priVelECF);
[priPosVNC,priVelVNC,transVNC] = rv2vnc(priPosECI(:,1),priVelECI(:,1));
[priPosRTC,priVelRTC,transRTC] = rv2rtc(priPosECI(:,1),priVelECI(:,1));
if j == count % when primary TLE is Nth TLE of time window
j
truthPosECI = priPosECI/1000
truthVelECI = priVelECI/1000
truthPosVNC = priPosVNC/1000
truthVelVNC = priVelVNC/1000
truthPosRTC = priPosRTC/1000
truthVelRTC = priVelRTC/1000
truthtransVNC = transVNC
truthtransRTC = transRTC
end;
%%%%%% Inner loop: Start with the last TLE in the time window, go backward.
for k = count:-1:1
if (j ~= k) % Don't compare the same TLE to itself
kName = ['TLE',TwoDigit(k)];
kEpochSec = atbTLEEpoch(kName);
if (kEpochSec < -0.001) % Only compare TLE if more
% than 0.001 seconds away
m = m + 1; % continuos counter for every j,k pair
% Calculate the number of days that kEpoch is away from jEpoch. This is
% the epoch time difference between the "K-th" TLE from the "J-th" TLE
delta(m) = -kEpochSec/86400;
% Notice below for the propagation of the "K-th" TLE, atbTLEProp has both
% a TimeStart and a TimeStop of 0.0 sec, meaning that the TLE will be
% propagated to the epoch tim, which was assigned above to the "J-th" TLE
[secTimes, secPosECF, secVelECF] = atbTLEProp(kName,0.0,0.0,1.0);
[secPosECI, secVelECI] = ...
atbCbfToCbi('Earth',secTimes,secPosECF,secVelECF);
% Calculate the relative position and velocity of the "K-th" TLE with
% respect to the "J-th" TLE, in the ECI coordinate frame
relPosECI = secPosECI(:,1) - priPosECI(:,1);
relVelECI = secVelECI(:,1) - priVelECI(:,1);
%old name: PosVNC(:,m) = transVNC*relPosECI/1000;
% Then convert these relative position & velocity vectors to the
% satellite-based coordinate system of the "J-th" TLE by using the
% transformation matrix from above (i.e. see the "rv2rtc" function).
% (Because vectors from atbTLEProp are in meters and
% meters/sec, divide the vectors below by 1000 to put in km and km/sec)
PosDiffVNC(:,m) = transVNC*relPosECI/1000;
VelDiffVNC(:,m) = transVNC*relVelECI/1000;
PosDiffRTC(:,m) = transRTC*relPosECI/1000;
VelDiffRTC(:,m) = transRTC*relVelECI/1000;
% Build data matrix of TLE numbers, TLE epochs, time differences,
% Positions and Velocities in RTC coordinate system for j and k TLES,
% Positions and Velocities in VNC coordinate system for j and k TLES,
% as well as their differences
data(1,1,j,k) = j;
data(1,2,j,k) = k;
data(1,3,j,k) = delta(m);
data(2,1,j,k) = jEpochSec;
data(2,2,j,k) = jEpochSec + kEpochSec;
data(2,3,j,k) = kEpochSec;
% RTC data
data((3:5),1,j,k) = priPosRTC(:,1)/1000;
data((3:5),2,j,k) = (priPosRTC(:,1)/1000) + PosDiffRTC(:,m);
data((3:5),3,j,k) = PosDiffRTC(:,m);
data((6:8),1,j,k) = priVelRTC(:,1)/1000;
data((6:8),2,j,k) = (priVelRTC(:,1)/1000) + VelDiffRTC(:,m);
data((6:8),3,j,k) = VelDiffRTC(:,m);
data((9:11),(1:3),j,k) = transRTC(:,:);
% VNC data
data((12:14),1,j,k) = priPosVNC(:,1)/1000;
data((12:14),2,j,k) = (priPosVNC(:,1)/1000) + PosDiffVNC(:,m);
data((12:14),3,j,k) = PosDiffVNC(:,m);
data((15:17),1,j,k) = priVelVNC(:,1)/1000;
data((15:17),2,j,k) = (priVelVNC(:,1)/1000) + VelDiffVNC(:,m);
data((15:17),3,j,k) = VelDiffVNC(:,m);
data((18:20),(1:3),j,k) = transVNC(:,:);
% Separate the Position Difference data into bins depending on the epoch
% time difference (in days). Separate the data into bins as follows:
% Bin 1: 0.0 - 0.49 days
% Bin 2: 0.5 - 1.49 days
%
% Bin 15: 13.5 - 14.49 days
bin = floor(delta(m) + 1.5);
if (bin < 16)
points(bin) = points(bin) + 1;
sumRTC(bin,1) = sumRTC(bin,1) + PosDiffRTC(1,m);
sumRTC(bin,2) = sumRTC(bin,2) + PosDiffRTC(2,m);
sumRTC(bin,3) = sumRTC(bin,3) + PosDiffRTC(3,m);
sumRTC2(bin,1) = sumRTC2(bin,1) + (PosDiffRTC(1,m))^2;
sumRTC2(bin,2) = sumRTC2(bin,2) + (PosDiffRTC(2,m))^2;
sumRTC2(bin,3) = sumRTC2(bin,3) + (PosDiffRTC(3,m))^2;
sumVNC(bin,1) = sumVNC(bin,1) + PosDiffVNC(1,m);
sumVNC(bin,2) = sumVNC(bin,2) + PosDiffVNC(2,m);
sumVNC(bin,3) = sumVNC(bin,3) + PosDiffVNC(3,m);
sumVNC2(bin,1) = sumVNC2(bin,1) + (PosDiffVNC(1,m))^2;
sumVNC2(bin,2) = sumVNC2(bin,2) + (PosDiffVNC(2,m))^2;
sumVNC2(bin,3) = sumVNC2(bin,3) + (PosDiffVNC(3,m))^2;
end; % if bin < 16
end; % if delta between J and K is negative
end; % if j <> k
end; %k
end; % for j
end; % if there is TLE data after the start date
end; % if file exists
fclose(fi);
% ---------------------------------------------------------------------
sat_alt_at_epoch = mag(truthPosECI)-6378;
% [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] =...
% RV2COE(truthPosECI(:,1),truthVelECI(:,1))
% [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] =...
% RV2COE(priPosECI(:,1)/1000,priVelECI(:,1)/1000)
sat_alt_at_epoch
min_delta_epoch = min(delta)
max_delta_epoch = max(delta)
% ---------------------------------------------------------------------
% This section calculates the variance of only the POSITION data
% differences in both the RTC and VNC coordinate frames FOR THE ENTIRE 2
% WEEK TIME WINDOW, USING THE BINS TO CALCULATE MEAN, VARIANCE, AND ST.DEV
RTCmean = zeros(15,3); % Initialize the 6 matrices
RTCvar = zeros(15,3);
RTCstdev = zeros(15,3);
VNCmean = zeros(15,3);
VNCvar = zeros(15,3);
VNCstdev = zeros(15,3);
for i = 1:15
if (points(i) > 1)
RTCmean(i,:) = sumRTC(i,:)/points(i);
RTCvar(i,1) = (sumRTC2(i,1) - points(i)*RTCmean(i,1)^2)/(points(i) - 1);
RTCvar(i,2) = (sumRTC2(i,2) - points(i)*RTCmean(i,2)^2)/(points(i) - 1);
RTCvar(i,3) = (sumRTC2(i,3) - points(i)*RTCmean(i,3)^2)/(points(i) - 1);
RTCstdev(i,1) = sqrt(RTCvar(i,1));
RTCstdev(i,2) = sqrt(RTCvar(i,2));
RTCstdev(i,3) = sqrt(RTCvar(i,3));
VNCmean(i,:) = sumVNC(i,:)/points(i);
VNCvar(i,1) = (sumVNC2(i,1) - points(i)*VNCmean(i,1)^2)/(points(i) - 1);
VNCvar(i,2) = (sumVNC2(i,2) - points(i)*VNCmean(i,2)^2)/(points(i) - 1);
VNCvar(i,3) = (sumVNC2(i,3) - points(i)*VNCmean(i,3)^2)/(points(i) - 1);
VNCstdev(i,1) = sqrt(VNCvar(i,1));
VNCstdev(i,2) = sqrt(VNCvar(i,2));
VNCstdev(i,3) = sqrt(VNCvar(i,3));
else
RTCmean(i,:) = NaN;
RTCvar(i,:) = NaN;
RTCstdev(i,:) = NaN;
VNCmean(i,:) = NaN;
VNCvar(i,:) = NaN;
VNCstdev(i,:) = NaN;
end; % if points(i) > 1
end; % for i
stats = [points, RTCmean, RTCvar,RTCstdev, VNCmean, VNCvar, VNCstdev];
augment = [delta;PosDiffRTC;VelDiffRTC;PosDiffVNC;VelDiffVNC];
data_at_epoch = data(:,:,count,:);
TLEsTimeInts = squeeze(data_at_epoch(1,3,1,:))
points, RTCmean, RTCvar, RTCstdev, VNCmean, VNCvar, VNCstdev,
% ---------------------------------------------------------------------
% This section calculates the covariance matrix from the position AND
% VELOCITY differences in the RTC and VNC coordinate frames.
% ** BUT ONLY AT THE EPOCH OF THE LAST TLE (most current).
% So only position differences from each TLE propagated to the last TLE,
% using that TLE as "truth" value.
RTCtemp = squeeze(data_at_epoch((3:8),3,1,(1:count-1)));
RTCmean_at_epoch = mean(RTCtemp')';
VNCtemp = squeeze(data_at_epoch((12:17),3,1,(1:count-1)));
VNCmean_at_epoch = mean(VNCtemp')';
RTCcov_at_epoch = zeros(6,6);
VNCcov_at_epoch = zeros(6,6);
for q = 1:(count-1)
RTCcov_at_epoch(:,:) = RTCcov_at_epoch(:,:) + ...
[(RTCtemp(:,q)-RTCmean_at_epoch(:))*...
(RTCtemp(:,q)-RTCmean_at_epoch(:))'];
VNCcov_at_epoch(:,:) = VNCcov_at_epoch(:,:) + ...
[(VNCtemp(:,q)-VNCmean_at_epoch(:))*...
(VNCtemp(:,q)-VNCmean_at_epoch(:))'];
end;
% (6 x 6) Covariance Matrices at epoch of more recent TLE using all the
% other TLEs propageted only to this epoch. In RTC and VNC coord. frames
RTCcov_at_epoch = (1/(count-1)) * RTCcov_at_epoch
VNCcov_at_epoch = (1/(count-1)) * VNCcov_at_epoch
% (6 x 1) column vector of the variances for the components of the
% RTC and VNC coordinate frames
RTCvariances_at_epoch = diag(RTCcov_at_epoch)
VNCvariances_at_epoch = diag(VNCcov_at_epoch)
% (6 x 1) column vector of standard deviations for the components of the
% RTC and VNC coordinate frames
RTCstdev_at_epoch = sqrt(RTCvariances_at_epoch)
VNCstdev_at_epoch = sqrt(VNCvariances_at_epoch)
% Calculate the correlation coefficients for this covariance matrix
for i = 1:6
for j = 1:6
if (i ~= j)
RTCcovar_corr_coeff(i,j) = RTCcov_at_epoch(i,j)...
/((sqrt(RTCcov_at_epoch(i,i)))*(sqrt(RTCcov_at_epoch(j,j))));
VNCcovar_corr_coeff(i,j) = VNCcov_at_epoch(i,j)...
/((sqrt(VNCcov_at_epoch(i,i)))*(sqrt(VNCcov_at_epoch(j,j))));
end;
if (i == j)
RTCcovar_corr_coeff(i,j) = 1;
VNCcovar_corr_coeff(i,j) = 1;
end;
end;
end;
RTCcovar_corr_coeff
VNCcovar_corr_coeff
% -----------------------------------------------------------------------
% Take the covariance_matrix_at_epoch in the satellite-based coordinate
% system, first rotate if back to cartesian coordinates (ECI), then
% transform it to a classical orbital elements covariance
% using Vallado's script, "COVCT2CL".
% ** This must use the "truthtrans" transformation matrix since it
% corresponds to the Nth TLE coordinate system.
cartstate = [truthPosECI;truthVelECI];
cartstate = cartstate*1000; % need to be passed in meters and meters/second
% from RTC
ECIcov1 = [truthtransRTC'*RTCcov_at_epoch((1:3),(1:3))*truthtransRTC,...
truthtransRTC'*RTCcov_at_epoch((1:3),(4:6))*truthtransRTC;...
truthtransRTC'*RTCcov_at_epoch((4:6),(1:3))*truthtransRTC,...
truthtransRTC'*RTCcov_at_epoch((4:6),(4:6))*truthtransRTC];
cartcov1 = ECIcov1*1000; % need to be passed in meters and meters/second
% from VNC
ECIcov2 = [truthtransVNC'*VNCcov_at_epoch((1:3),(1:3))*truthtransVNC,...
truthtransVNC'*VNCcov_at_epoch((1:3),(4:6))*truthtransVNC;...
truthtransVNC'*VNCcov_at_epoch((4:6),(1:3))*truthtransVNC,...
truthtransVNC'*VNCcov_at_epoch((4:6),(4:6))*truthtransVNC];
cartcov2 = ECIcov2*1000; % need to be passed in meters and meters/second
[COEcov1mean,tm1mean] = COVCT2CL(cartcov1,cartstate,'mean');
[COEcov1true,tm1true] = COVCT2CL(cartcov1,cartstate,'true');
[COEcov2mean,tm2mean] = COVCT2CL(cartcov2,cartstate,'mean');
[COEcov2true,tm2true] = COVCT2CL(cartcov2,cartstate,'true');
COE_variances_fromRTC_meananomaly = diag(COEcov1mean)
COE_variances_fromRTC_trueanomaly = diag(COEcov1true)
COE_variances_fromVNC_meananomaly = diag(COEcov2mean)
COE_variances_fromVNC_trueanomaly = diag(COEcov2true)
COEcheck1 = COEcov1mean - COEcov2mean
COEcheck2 = COEcov1true - COEcov2true
% ---------------------------------------
% Build giant matrix of covariance matrices and variances for all the
% respective coordinate systems
Epoch_Cov_All =[RTCcov_at_epoch RTCvariances_at_epoch RTCstdev_at_epoch;...
VNCcov_at_epoch VNCvariances_at_epoch VNCstdev_at_epoch;...
RTCcovar_corr_coeff zeros(6,2);...
VNCcovar_corr_coeff zeros(6,2);...
COEcov1mean COE_variances_fromRTC_meananomaly zeros(6,1);...
COEcov1true COE_variances_fromRTC_trueanomaly zeros(6,1);...
COEcov2mean COE_variances_fromVNC_meananomaly zeros(6,1);...
COEcov2true COE_variances_fromVNC_trueanomaly zeros(6,1)];
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Write data to Excel file
if (writeexcel)
AugCol = {'Delta-t TLE Epochs (days)';'Position Difference Radial (km)';...
'Position Difference Transverse (km)';
'Position Difference Cross-track (km)';...
'Velocity Difference Radial (km/sec)';...
'Velocity Difference Transverse (km/sec)';...
'Velocity Difference Cross-track (km/sec)';...
'Position Difference Velocity (km)';...
'Position Difference Normal (km)';...
'Position Difference Cross-track (km)';...
'Velocity Difference Velocity (km/sec)';...
'Velocity Difference Normal (km/sec)';...
'Velocity Difference Cross-track (km/sec)'};
StatRow1={'Statistics by Bin Number for only Position Differences (km)',...
'Mean values','in RTC frame','by component',...
'Variance values','in RTC frame','by component',...
'StDev values','in RTC frame','by component',...
'Mean values','in VNC frame','by component',...
'Variance values','in VNC frame','by component',...
'StDev values','in VNC frame','by component'};
StatRow2 = {'Sample Points in Bin',...
'Radial','Transverse','Cross Track',...
'Radial','Transverse','Cross Track',...
'Radial','Transverse','Cross Track',...
'Velocity','Normal','Cross Track',...
'Velocity','Normal','Cross Track',...
'Velocity','Normal','Cross Track'};
SheetDescrip = {'Jth TLE Number','Kth TLE Number',...
'Delta Time between epochs (days)';...
'Jth TLE epoch (sec)','Kth TLE epoch (sec)',...
'Delta Time between epochs (sec)';...
'Jth Radial Position (km)','Kth Radial Position (km)',...
'Difference Radial Position (km)';...
'Jth Transverse Position (km)','Kth Transverse Position (km)',...
'Difference Transverse Position (km)';...
'Jth Crosstrack Position (km)','Kth Crosstrack Position (km)',...
'Difference Crosstrack Position (km)';...
'Jth Radial Velocity (km/sec)','Kth Radial Velocity (km/sec)',...
'Difference Radial Velocity (km/sec)';...
'Jth Transverse Velocity (km/sec)','Kth Transverse Velocity (km/sec)',...
'Difference Transverse Velocity (km/sec)';...
'Jth Crosstrack Velocity (km/sec)','Kth Crosstrack Velocity (km/sec)',...
'Difference Crosstrack Velocity (km/sec)';...
'RTC transformation matrix - R1','R2','R3';'T1','T2','T3';...
'C1','C2','C3';...
'Jth Velocity Position (km)','Kth Velocity Position (km)',...
'Difference Velocity Position (km)';...
'Jth Normal Position (km)','Kth Normal Position (km)',...
'Difference Normal Position (km)';...
'Jth Crosstrack Position (km)','Kth Crosstrack Position (km)',...
'Difference Crosstrack Position (km)';...
'Jth Velocity Velocity (km/sec)','Kth Velocity Velocity (km/sec)',...
'Difference Velocity Velocity (km/sec)';...
'Jth Normal Velocity (km/sec)','Kth Normal Velocity (km/sec)',...
'Difference Normal Velocity (km/sec)';...
'Jth Crosstrack Velocity (km/sec)','Kth Crosstrack Velocity (km/sec)',...
'Difference Crosstrack Velocity (km/sec)';...
'VNC transformation matrix - V1','V2','V3';'N1','N2','N3';...
'C1','C2','C3'};
CovDes1 = {'POSITIONradial_POSITIONradial',...
'Pr_Pt','Pr_Pc','Pr_Vr','Pr_Vt','Pr_Vc';...
'Pt_Pr','Ptransverse_Ptransverse','Pt_Pc','Pt_Vr','Pt_Vt','Pt_Vc';...
'Pc_Pr','Pc_Pt','Pcrosstrack_Pcrosstrack','Pc_Vr','Pc_Vt','Pc_Vc';...
'Vr_Pr','Vr_Pt','Vr_Pc','Vradial_Vradial','Vr_Vt','Vr_Vc';...
'Vt_Pr','Vt_Pt','Vt_Pc','Vt_Vr','Vtransverse_Vtransverse','Vt_Vc';...
'Vc_Pr','Vc_Pt','Vc_Pc','Vc_Vr','Vc_Vt',...
'VELOCITYcrosstrack_VELOCITYcrosstrack';...
'POSITIONvelocity_POSITIONvelocity',...
'Pv_Pn','Pv_Pc','Pv_Vv','Pv_Vn','Pv_Vc';...
'Pn_Pv','Pnormal_Pnormal','Pn_Pc','Pn_Vv','Pn_Vn','Pn_Vc';...
'Pc_Pv','Pc_Pn','Pcrosstrack_Pcrosstrack','Pc_Vv','Pc_Vn','Pc_Vc';...
'Vv_Pv','Vv_Pn','Vv_Pc','Vvelocity_Vvelocity','Vv_Vn','Vv_Vc';...
'Vn_Pv','Vn_Pn','Vn_Pc','Vn_Vv','Vnormal_Vnormal','Vn_Vc';...
'Vc_Pv','Vc_Pn','Vc_Pc','Vc_Vv','Vc_Vn',...
'VELOCITYcrosstrack_VELOCITYcrosstrack'};
CovDes2 = {'6x6','RTC','Covariance','Correlation','Coefficients','Matrix'};
CovDes3 = {'6x6','VNC','Covariance','Correlation','Coefficients','Matrix'};
CovDes4 = {'R';'T';'C';'Variances';'at';'epoch';...
'V';'N';'C';'Variances';'at';'epoch'};
CovDes5 = {'R';'T';'C';'StDev';'at';'epoch';...
'V';'N';'C';'StDev';'at';'epoch'};
CovDes6 = {'6x6','COE','Covariance','Matrix','from RTC transform',...
'using Mean Anomaly','Variances Column Here','Zero fill column'};
CovDes7 = {'6x6','COE','Covariance','Matrix','from RTC transform',...
'using True Anomaly','Variances Column Here','Zero fill column'};
CovDes8 = {'6x6','COE','Covariance','Matrix','from VNC transform',...
'using Mean Anomaly','Variances Column Here','Zero fill column'};
CovDes9 = {'6x6','COE','Covariance','Matrix','from VNC transform',...
'using True Anomaly','Variances Column Here','Zero fill column'};
xfile = ['I:\My Documents\Thesis\Output\data_sat',...
catnr,'_t',num2str(t),'.xls'];
for ex = 1:count
xlswrite(xfile, data_at_epoch(:,:,1,ex),ex);
end;
xlswrite(xfile, SheetDescrip, 'SheetXX Description');
% xlswrite(xfile, TLEsTimeInts, 'times');
xlswrite(xfile, StatRow1, 'Statistics','A1');
xlswrite(xfile, StatRow2, 'Statistics','A2');
xlswrite(xfile, stats, 'Statistics','A3');
xlswrite(xfile, AugCol', 'AugmentMatrix','A1:M1');
xlswrite(xfile, augment', 'AugmentMatrix','A2');
xlswrite(xfile, Epoch_Cov_All, 'Epoch Covariance Data');
xlswrite(xfile, CovDes1, 'Epoch_Cov_All Sheet','A1');
xlswrite(xfile, CovDes2, 'Epoch_Cov_All Sheet','A13');
xlswrite(xfile, CovDes3, 'Epoch_Cov_All Sheet','A19');
xlswrite(xfile, CovDes4, 'Epoch_Cov_All Sheet','G1');
xlswrite(xfile, CovDes5, 'Epoch_Cov_All Sheet','H1');
xlswrite(xfile, CovDes6, 'Epoch_Cov_All Sheet','A25');
xlswrite(xfile, CovDes7, 'Epoch_Cov_All Sheet','A31');
xlswrite(xfile, CovDes8, 'Epoch_Cov_All Sheet','A37');
xlswrite(xfile, CovDes9, 'Epoch_Cov_All Sheet','A43');
end; % if writeexcel files
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% PLOTTING
if(toplotornottoplot) % PLOT RUN #1
x = 0.0:0.5:15.0;
b = 0:1:14;
if (logfit)
pV = polyfit(delta,log10(abs(PosDiffVNC(1,:))),2);
pN = polyfit(delta,log10(abs(PosDiffVNC(2,:))),2);
pC = polyfit(delta,log10(abs(PosDiffVNC(3,:))),2);
else
pV = polyfit(delta,(PosDiffVNC(1,:)),2);
pN = polyfit(delta,(PosDiffVNC(2,:)),2);
pC = polyfit(delta,(PosDiffVNC(3,:)),2);
end; % if logfit
if (logfit)
fitV = 10.^polyval(pV,x);
fitN = 10.^polyval(pN,x);
fitC = 10.^polyval(pC,x);
else
fitV = polyval(pV,x);
fitN = polyval(pN,x);
fitC = polyval(pC,x);
end; % if logfit
figure(fignr*10+1)
clf; figurenr = gcf;
subplot(3,1,1)
hold on;
title({'Positional Differences in NTW (VNC) Coordinate Frame';...
['NORAD Catalog Number ',catnr];...
['Time Window: ',datestr(startdate),' to ',datestr(enddate)]});
plot(delta,PosDiffVNC(1,:),'r.');
plot(x,fitV,'r:');
errorbar(b,VNCmean(:,1),VNCstdev(:,1),'ro','MarkerSize',6);
legend('Velocity (in-track)','Location','Best');
% axis tight;
xlim([0 15]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
subplot(3,1,2)
hold on;
plot(delta,PosDiffVNC(2,:),'b.');
plot(x,fitN,'b:');
errorbar(b,VNCmean(:,2),VNCstdev(:,2),'bo','MarkerSize',6);
legend('Normal (along-radial) ','Location','Best');
axis tight; xlim([0 15]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
subplot(3,1,3)
hold on;
plot(delta,PosDiffVNC(3,:),'k.');
plot(x,fitC,'k:');
errorbar(b,VNCmean(:,3),VNCstdev(:,3),'ko','MarkerSize',6);
legend('Cross-track','Location','Best');
% axis tight;
xlim([0 15]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\figure_sat',...
catnr,'_t',num2str(t),'_PosDiffVNC_subplots.jpg']);
% -----------------------------------------------------------
clear pC fitC;
if (logfit)
pR = polyfit(delta,log10(abs(PosDiffRTC(1,:))),2);
pT = polyfit(delta,log10(abs(PosDiffRTC(2,:))),2);
pC = polyfit(delta,log10(abs(PosDiffRTC(3,:))),2);
else
pR = polyfit(delta,(PosDiffRTC(1,:)),2);
pT = polyfit(delta,(PosDiffRTC(2,:)),2);
pC = polyfit(delta,(PosDiffRTC(3,:)),2);
end; % if logfit
if (logfit)
fitR = 10.^polyval(pR,x);
fitT = 10.^polyval(pT,x);
fitC = 10.^polyval(pC,x);
else
fitR = polyval(pR,x);
fitT = polyval(pT,x);
fitC = polyval(pC,x);
end; % if logfit
figure(fignr*10+2)
clf; figurenr = gcf;
subplot(3,1,1)
hold on;
title({'Positional Differences in RSW Coordinate Frame';...
['NORAD Catalog Number ',catnr];...
['Time Window: ',datestr(startdate),' to ',datestr(enddate)]},...
'Fontsize',12);
plot(delta,PosDiffRTC(2,:),'b.');
plot(x,fitT,'b:');
errorbar(b,RTCmean(:,2),RTCstdev(:,2),'bo','MarkerSize',6);
legend('Transverse (along-track)','Location','Best');
% axis tight;
xlim([0 15]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
subplot(3,1,2)
hold on;
plot(delta,PosDiffRTC(1,:),'r.');
plot(x,fitR,'r:');
errorbar(b,RTCmean(:,1),RTCstdev(:,1),'ro','MarkerSize',6);
legend('Radial','Location','Best');
% axis tight;
xlim([0 15]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
subplot(3,1,3)
hold on;
plot(delta,PosDiffRTC(3,:),'k.');
plot(x,fitC,'k:');
errorbar(b,RTCmean(:,3),RTCstdev(:,3),'ko','MarkerSize',6);
legend('Cross-track','Location','Best');
% axis tight;
xlim([0 15]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\figure_sat',...
catnr,'_t',num2str(t),'_PosDiffRTC_subplots.jpg']);
% --------------------------------------------------------
clear pC fitC;
if (logfit)
pV = polyfit(delta,log10(abs(PosDiffVNC(1,:))),2);
pN = polyfit(delta,log10(abs(PosDiffVNC(2,:))),2);
pC = polyfit(delta,log10(abs(PosDiffVNC(3,:))),2);
else
pV = polyfit(delta,(PosDiffVNC(1,:)),2);
pN = polyfit(delta,(PosDiffVNC(2,:)),2);
pC = polyfit(delta,(PosDiffVNC(3,:)),2);
end; % if logfit
figure(fignr*10+3)
clf; figurenr = gcf; hold on;
% plot(delta,(PosDiffVNC),'.');
plot(delta,PosDiffVNC(1,:),'rd');
plot(delta,PosDiffVNC(2,:),'bo');
plot(delta,PosDiffVNC(3,:),'kx');
if (logfit)
fitV = 10.^polyval(pV,x);
fitN = 10.^polyval(pN,x);
fitC = 10.^polyval(pC,x);
else
fitV = polyval(pV,x);
fitN = polyval(pN,x);
fitC = polyval(pC,x);
end; % if logfit
plot(x,fitV,'r:',x,fitN,'b:',x,fitC,'k:');
errorbar(b,VNCmean(:,1),VNCstdev(:,1),'r*','MarkerSize',6);
errorbar(b,VNCmean(:,2),VNCstdev(:,2),'b*','MarkerSize',6);
errorbar(b,VNCmean(:,3),VNCstdev(:,3),'k*','MarkerSize',6);
title({'Positional Differences in NTW (VNC) Coordinate Frame';...
['NORAD Catalog Number ',catnr];...
['Time Window: ',datestr(startdate),' to ',datestr(enddate)]},...
'Fontsize',12);
legend('Velocity (in-track)','Normal (along-radial)','Cross-track',...
'Location','Best');
% axis tight;
xlim([0 15]);
%ylim([-1 1]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\figure_sat',...
catnr,'_t',num2str(t),'_PosDiffVNC.jpg']);
% -------------------------------------------------------
clear pC fitC;
if (logfit)
pR = polyfit(delta,log10(abs(PosDiffRTC(1,:))),2);
pT = polyfit(delta,log10(abs(PosDiffRTC(2,:))),2);
pC = polyfit(delta,log10(abs(PosDiffRTC(3,:))),2);
else
pR = polyfit(delta,(PosDiffRTC(1,:)),2);
pT = polyfit(delta,(PosDiffRTC(2,:)),2);
pC = polyfit(delta,(PosDiffRTC(3,:)),2);
end; % if logfit
figure(fignr*10+4)
clf; figurenr = gcf; hold on;
% plot(delta,(PosDiffRTC),'.');
plot(delta,PosDiffRTC(1,:),'rd');
plot(delta,PosDiffRTC(2,:),'bo');
plot(delta,PosDiffRTC(3,:),'kx');
if (logfit)
fitR = 10.^polyval(pR,x);
fitT = 10.^polyval(pT,x);
fitC = 10.^polyval(pC,x);
else
fitR = polyval(pR,x);
fitT = polyval(pT,x);
fitC = polyval(pC,x);
end; % if logfit
plot(x,fitR,'b:',x,fitT,'r:',x,fitC,'k:');
errorbar(b,RTCmean(:,1),RTCstdev(:,1),'r*','MarkerSize',6);
errorbar(b,RTCmean(:,2),RTCstdev(:,2),'b*','MarkerSize',6);
errorbar(b,RTCmean(:,3),RTCstdev(:,3),'k*','MarkerSize',6);
title({'Positional Differences in RSW Coordinate Frame';...
['NORAD Catalog Number ',catnr];...
['Time Window: ',datestr(startdate),' to ',datestr(enddate)]},...
'Fontsize',12);
legend('Radial','Transverse (along-track)','Cross-track',...
'Location','Best');
% axis tight;
xlim([0 15]);
%ylim([-1 1]);
set(gca,'XTick',0:1:15); set(gca,'YGrid','on');
xlabel('Delta Epoch (days)'); ylabel('Delta Position (km)');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\figure_sat',...
catnr,'_t',num2str(t),'_PosDiffRTC.jpg']);
% ----------------------------------------------------------
figure(fignr*10+5)
clf;
figurenr = gcf;
subplot(2,1,1)
hold on;
plot(RTCvar(:,1),'ro-','LineWidth',1);
plot(RTCvar(:,3),'gd:','LineWidth',1);
plot(VNCvar(:,2),'bx--','LineWidth',1);
plot(VNCvar(:,3),'k*-.','LineWidth',1);
% text(0,0,['Number of sample points per bin = ',num2str(points(:)')],...
% 'HorizontalAlignment','left','VerticalAlignment','top',...
% 'BackgroundColor',[.7 .9 .7],'FontSize',16);
title({'Variance in Positional Differences';...
['for Entire Time Window: ',datestr(startdate),' to ',...
datestr(enddate)];...
['NORAD Catalog Number ',catnr]},'Fontsize',12);
legend('Radial','Cross-track','Normal (along-radial)',...
'Cross-track','Location','NorthWest');
% axis tight;
xlabel('Bin Number'); ylabel('Variance (km^2)');
set(gca,'XTick',0:1:15);
set(gca,'YGrid','on');
hold off;
subplot(2,1,2)
hold on;
title({'Number of sample points per bin = ';num2str(points(:)')},...
'FontSize',12);
plot(RTCvar(:,2),'bo-','LineWidth',1);
plot(VNCvar(:,1),'rx-','LineWidth',1);
legend('Transverse (along-track)','Velocity (in-track)','Location',...
'NorthWest');
% axis tight;
xlabel('Bin Number'); ylabel('Variance (km^2)');
set(gca,'XTick',0:1:15);
set(gca,'YGrid','on');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\figure_sat',...
catnr,'_t',num2str(t),'_Variances.jpg']);
end; % if toplotornottoplot #1 run
% END PLOTTING SECTION #1
% --------------------------------------------------------------------
% Save the Variance Values from the differenct covariance matrices
% computed above, all at the primary TLE epoch.
% Save for each timeloop.
RTC_Cov_Total(:,:,timeloop) = RTCcov_at_epoch;
VNC_Cov_Total(:,:,timeloop) = VNCcov_at_epoch;
RTC_Var_Total(:,timeloop) = RTCvariances_at_epoch;
VNC_Var_Total(:,timeloop) = VNCvariances_at_epoch;
COE_Var_total(:,timeloop) = COE_variances_fromVNC_meananomaly;
% ----------------------------------------------------
% Create Autocorrelation matrices for each timeloop
% They will be used after timeloop finishes
%%%% for RTC
% TLEperTimeWindow(timeloop) = count
% tempo1(timeloop) = size(PosDiffRTC,2)
% padamount = 7000 - tempo1(timeloop);
% tempAC1 = [delta' PosDiffRTC'];
% tempAC2 = padarray(tempAC1,[padamount 0],'post');
% AutoCorr(:,:,timeloop) = tempAC2;
%%%% for VNC
TLEperTimeWindow(timeloop) = count
tempo1(timeloop) = size(PosDiffVNC,2)
padamount = 7000 - tempo1(timeloop);
tempAC1 = [delta' PosDiffVNC'];
tempAC2 = padarray(tempAC1,[padamount 0],'post');
AutoCorr(:,:,timeloop) = tempAC2;
% ----------------------------------------------------
% close all; %close all figures --> they should be saving to hard drive
diary off;
end; % for timeloop 1 through 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
if(toplotornottoplot) % PLOT RUN #2
figure(fignr*10+6)
clf; figurenr = gcf;
subplot(3,1,1)
hold on;
title({'Variances at Primary TLE Epoch for all 8 Time Windows';...
['NORAD Catalog Number ',catnr]});
plot(RTC_Var_Total(2,:),'bo-','LineWidth',2);
legend('Transverse (along-track) Position','Location','BestOutside');
xlim([1 8]); set(gca,'XTick',0:1:8);
set(gca,'YGrid','on'); ylabel('Variance (km^2)');
hold off;
subplot(3,1,2)
hold on;
plot(RTC_Var_Total(1,:),'rd-','LineWidth',2);
plot(RTC_Var_Total(3,:),'kx-','LineWidth',2);
legend('Radial Postion','Cross-track Position','Location','BestOutside');
xlim([1 8]); set(gca,'XTick',0:1:8);
set(gca,'YGrid','on'); ylabel('Variance (km^2)');
hold off;
subplot(3,1,3)
hold on;
plot(RTC_Var_Total(4,:),'rd-','LineWidth',2);
plot(RTC_Var_Total(5,:),'bo-','LineWidth',2);
plot(RTC_Var_Total(6,:),'kx-','LineWidth',2);
legend('Radial Velocity','Transverse Velocity','Cross-track Velocity',...
'Location','BestOutside');
xlim([1 8]); set(gca,'XTick',0:1:8); set(gca,'YGrid','on');
xlabel('Time Window #'); ylabel('Variance (km^2/sec^2)');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\figure_sat',...
catnr,'_EpochVariances_subplots.jpg']);
end; % if toplotornottoplot #2 run
% ---------------------------------------------------------------------**
% ***********************************************************************
% AUTOCORRELATION FUNCTION
% ***********************************************************************
% ---------------------------------------------------------------------**
if(autocorrelate)
dout2 = ['I:\My Documents\Thesis\Output\AC_output_sat_',...
num2str(satnr),'.txt'];
diary(dout2)
ACin1 = AutoCorr((1:tempo1(1)),:,1);
ACin2 = AutoCorr((1:tempo1(2)),:,2);
ACin3 = AutoCorr((1:tempo1(3)),:,3);
ACin4 = AutoCorr((1:tempo1(4)),:,4);
ACin5 = AutoCorr((1:tempo1(5)),:,5);
ACin6 = AutoCorr((1:tempo1(6)),:,6);
ACin7 = AutoCorr((1:tempo1(7)),:,7);
ACin8 = AutoCorr((1:tempo1(8)),:,8);
% Combine all of the position residuals from every time window
% into one giant matrix
giantAC = [ACin1' ACin2' ACin3' ACin4' ACin5' ACin6' ACin7' ACin8'];
AC_overall_mean = mean(giantAC')';
% FOLLOWING 2 OPTIONS TO CHOOSE FOR AUTOCORRELATION FUNCTION
% Binwidth and # of days
binwidth = 0.5; % Width of each bin set to half day
ndays = 35; % number of days to pull residuals from
nob = ndays/binwidth; % total number of bins
ACsamples = zeros(nob,1);
ACsum = zeros(nob,3);
ACsum2 = zeros(nob,3);
for w = 1:length(giantAC)
ACbin = floor((giantAC(1,w)/binwidth) + 1.0);
if (ACbin < (nob+1))
ACsamples(ACbin) = ACsamples(ACbin) + 1;
ACsum(ACbin,1) = ACsum(ACbin,1) + giantAC(2,w);
ACsum(ACbin,2) = ACsum(ACbin,2) + giantAC(3,w);
ACsum(ACbin,3) = ACsum(ACbin,3) + giantAC(4,w);
ACsum2(ACbin,1) = ACsum2(ACbin,1) + ((giantAC(2,w))^2);
ACsum2(ACbin,2) = ACsum2(ACbin,2) + ((giantAC(3,w))^2);
ACsum2(ACbin,3) = ACsum2(ACbin,3) + ((giantAC(4,w))^2);
end;
end;
ACsamples
figure(10+7)
plot(ACsamples);
% stopper = 10*(nob/15);
stopper = nob;
% Initialize the matrices
AC_mean_by_bin = zeros(nob,3);
for i = 1:nob
if (ACsamples(i) > 0)
% AC_mean_by_bin(i,:) = ACsum(i,:)/ACsamples(i); % uses residual
AC_mean_by_bin(i,:) = ACsum2(i,:)/ACsamples(i); % uses (residual^2)
elseif (ACsamples(i) == 0)
AC_mean_by_bin(i,1) = AC_overall_mean(2);
AC_mean_by_bin(i,2) = AC_overall_mean(3);
AC_mean_by_bin(i,3) = AC_overall_mean(4);
% AC_mean_by_bin(i,:) = NaN;
% AC_mean_by_bin(i,:) = 0.0;
end; % if ACsamples(i) > 1
end; % for i
% Look for "NaN" values, and don't use them in calculating a mean value for
% each bin
% rowtemp = 0;
% tempNan = isnan(AC_mean_by_bin)
% for i = 1:nob
% if (tempNan(i) == 1)
% rowtemp = rowtemp + 1
% AC_mean_by_bin(i,:) = []
% deleted_timesteps(rowtemp) = i
% end;
% end;
% Reassign the mean value of each bin to be the observations used for
% Autocorreation
AC_obs = AC_mean_by_bin
AC_obs_mean = mean(AC_obs)
% ---------------------------------------------------------------------
% AUTOCORRELATION computations %% COEFF
% Using "xcov" MATLAB cross covariance function
xcov_c1_full = xcov(AC_obs(:,1),'coeff');
xcov_c2_full = xcov(AC_obs(:,2),'coeff');
xcov_c3_full = xcov(AC_obs(:,3),'coeff');
% xcov_c1 = xcov_c1_full([stopper:(stopper+(ceil(nob/5)))])
% xcov_c2 = xcov_c2_full([stopper:(stopper+(ceil(nob/5)))])
% xcov_c3 = xcov_c3_full([stopper:(stopper+(ceil(nob/5)))])
xcov_c1 = xcov_c1_full([stopper:length(xcov_c1_full)])
xcov_c2 = xcov_c2_full([stopper:length(xcov_c2_full)])
xcov_c3 = xcov_c3_full([stopper:length(xcov_c3_full)])
% AC_xfile = ['I:\My Documents\Thesis\Output\AutoCorr_sat',catnr,'.xls'];
% xlswrite(AC_xfile, AC_final, 'AutoCorrelation');
% ---------------------------------------------------------------------
% AUTOCORRELATION plot
% xsteps_old = [0:binwidth:(stopper-1)/(nob/ndays)]
leend = length(xcov_c1);
xsteps = [0:binwidth:((leend-1)*binwidth)]
figure(fignr*10+8)
clf; figurenr = gcf; hold on;
plot(xsteps,xcov_c1,'ro-','LineWidth',2);
plot(xsteps,xcov_c2,'bd--','LineWidth',2);
plot(xsteps,xcov_c3,'kx:','LineWidth',2);
ylim([-1 1]);
title({'Normalized Autocorrelation of Positional Differences';...
['NORAD Catalog Number ',catnr];...
'All Time Windows Grouped with Bin Averages as Observations'},...
'Fontsize',12);
legend('Velocity (in-track)','Normal (along radial)','Cross-track','Location','Best');
%legend('Radial','Transverse (along-track)','Cross-track','Location','Best');
set(gca,'XTick',0:1:ndays);
set(gca,'YGrid','on');
xlabel(['Time Shift (days) with Bin Width of: ',num2str(binwidth),' days)']);
ylabel('Correlation');
hold off;
saveas(figurenr,['I:\My Documents\Thesis\Output\AC_sat',...
catnr,'_Coeff_not_squared.jpg']);
%saveas(figurenr,['I:\My Documents\Thesis\Output\AC_sat',...
% catnr,'_Coeff_VNC_thin_extraTime.jpg']);
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% AUTOCORRELATION computations %% UNBIASED
% Using "xcov" MATLAB cross covariance function
% xcov_c4 = xcov(AC_obs((1:stopper),1),'unbiased');
% xcov_c5 = xcov(AC_obs((1:stopper),2),'unbiased');
% xcov_c6 = xcov(AC_obs((1:stopper),3),'unbiased');
%
% xcov_c4 = xcov_c4(stopper:length(xcov_c4))
% xcov_c5 = xcov_c5(stopper:length(xcov_c5))
% xcov_c6 = xcov_c6(stopper:length(xcov_c6))
%
%
% % AC_xfile = ['I:\My Documents\Thesis\Output\AutoCorr_sat',catnr,'.xls'];
% % xlswrite(AC_xfile, AC_final, 'AutoCorrelation');
%
% % ---------------------------------------------------------------------
% % AUTOCORRELATION plot
%
% xsteps = [0:binwidth:ndays-binwidth]
%
% figure(fignr*10+9)
% clf; figurenr = gcf; hold on;
% plot(xsteps,xcov_c4,'ro-','LineWidth',2);
% plot(xsteps,xcov_c5,'bd-','LineWidth',2);
% plot(xsteps,xcov_c6,'kx-','LineWidth',2);
% title({'Unbiased Autocorrelation of Positional Differences';...
% ['NORAD Catalog Number ',catnr];...
% 'All Time Windows Grouped with Bin Averages as Observations'},...
% 'Fontsize',12);
% legend('Velocity (in-track)','Normal (along radial)','Cross-track','Location','Best');
% %legend('Radial','Transverse (along-track)','Cross-track','Location','Best');
% set(gca,'XTick',0:1:ndays);
% set(gca,'YGrid','on');
% xlabel('Time Shift (days)'); ylabel('Correlation');
% hold off;
%
% saveas(figurenr,['I:\My Documents\Thesis\Output\AC_sat',...
% catnr,'_Unbiased_VNC_thin_extraTime.jpg']);
% -------------------------------------------------------------------
% close all;
diary off;
end; % if autocorrelate
% ---------------------------------------------------------------------**
% *********************************************************************
% AUTOCORRELATION FUNCTION
% *********************************************************************
% ---------------------------------------------------------------------**
%%%% Now CLEAR the variables before incrementing to next satellite #
clear TLEperTimeWindow count padamount tempAC1 tempAC2;
clear data_at_epoch AutoCorr tempNan N giantAC AC_overall_mean;
clear ACin1 ACin2 ACin3 ACin4 ACin5 ACin6 ACin7 ACin8;
clear ACsamples deleted_timesteps AC_obs_mean number_of_shifts AC_obs_sq;
clear AC_final
end; % for satnr 1:6
% close all;
% clc;
% end of file