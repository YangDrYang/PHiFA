% -------------------------------------------------------------
function e2p = Epoch2Date(epoch);
year2 = str2num(epoch(1:2));
doy = str2num(epoch(3:14));
if (year2 > 56)
e2p = datenum(1900 + year2,1,1) + doy - 1;
else
e2p = datenum(2000 + year2,1,1) + doy - 1;
end; % if
% -------------------------------------------------------------
function td = TwoDigit(n);
if (n < 10)
td = ['0',int2str(n)];
else
td = int2str(n);
end; % if
% -------------------------------------------------------------
function mag = mag ( vec );
temp= vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3);
if abs( temp ) >= 1.0e-16
mag= sqrt( temp );
else
mag= 0.0;
end

% ------------------------------------------------------------
%
% function rv2vnc
%
% This function converts position and velocity vectors into
% Velocity, Normal, and Cross-track coordinate frame (VNC)
%
% NOTE: that sometimes the second vector is called along-radial.
%
% author : david vallado 719-573-2600 5 jul 2002
%
% revisions
% - Capt Victor Osweiler 719-310-1801 4 Feb 2006
%
% inputs description range / units
% r - position vector km
% v - velocity vector km/s
%
% outputs :
% rvnc - position vector km
% vvnc - velocity vector km/s
%
% locals :
% temp - temporary position vector
%
% coupling :
% mag - magnitude of a vector
%
% references :
% vallado 2001, xx
%
% [rvnc,vvnc,transmat] = rv2vnc( r,v );
% -----------------------------------------------------------------
function [rvnc,vvnc,transmat] = rv2vnc( r,v );
% compute satellite position vector magnitude
rmag = mag(r);
% compute satellite velocity vector magnitude
vmag = mag(v);
% in order to work correctly each of the components must be
% unit vectors !
% in-velocity component
vvec = v / vmag;
% cross-track component
cvec = cross(r,v);
cvec = unit( cvec );
% along-radial component
nvec = cross(vvec,cvec);
nvec = unit( nvec );
% assemble transformation matrix from to vnc frame (individual
% components arranged in row vectors)
transmat(1,1) = vvec(1);
transmat(1,2) = vvec(2);
transmat(1,3) = vvec(3);
transmat(2,1) = nvec(1);
transmat(2,2) = nvec(2);
transmat(2,3) = nvec(3);
transmat(3,1) = cvec(1);
transmat(3,2) = cvec(2);
transmat(3,3) = cvec(3);
rvnc = transmat*r;
vvnc = transmat*v;

% -------------------------------------------------------------
%
% function rv2rtc
%
% This function converts position and velocity vectors into a
% Radial, Transverse, and Cross-track coordinate system frame.
%
% Radial positions are parallel to the position vector (along R axis)
% Transverse displacements are normal to position vector (along S axis).
% Cross-track positions are normal to the plane defined by the current
% position and velocity vectors (along the W axis)
%
% NOTE: sometimes the second vector (Transverse) is called the along-track
%
% Modified file from original file: rv2ivc by Vallado
% original author : david vallado 719-573-2600 5 Jul 2002
%
% revisions
% - Capt Victor Osweiler 719-310-1801 24 Jan 2006
%
% inputs description range / units
% r - position vector km
% v - velocity vector km/s
%
% outputs :
% rrtc - position vector km
% vrtc - velocity vector km/s
%
% locals :
%
%
% coupling :
% mag - magnitude of a vector
% unit - calculates the unit vector
%
% references :
% vallado 2001, xx
%
% [rrtc,vrtc,transmat] = rv2rtc( r,v );
% ------------------------------------------------------------------
function [rrtc,vrtc,transmat] = rv2rtc( r,v );
% compute satellite position vector magnitude
rmag = mag(r);
% compute satellite velocity vector magnitude
vmag = mag(v);
% in order to work correctly each of the components must be
% unit vectors !
% Radial component
rvec = r / rmag;
% Normal component
cvec = cross(r,v);
cvec = unit( cvec ); % calls function "unit.m" to calculate
% Transverse component
tvec = cross(cvec,rvec);
tvec = unit( tvec );
% assemble transformation matrix from to ivc frame (individual
% components arranged in row vectors)
transmat(1,1) = rvec(1);
transmat(1,2) = rvec(2);
transmat(1,3) = rvec(3);
transmat(2,1) = tvec(1);
transmat(2,2) = tvec(2);
transmat(2,3) = tvec(3);
transmat(3,1) = cvec(1);
transmat(3,2) = cvec(2);
transmat(3,3) = cvec(3);
rrtc = transmat*r;
vrtc = transmat*v;
% -------------------------------------------------------------