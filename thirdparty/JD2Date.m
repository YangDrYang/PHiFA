function datetime = JD2Date( jd, structOut )

%-------------------------------------------------------------------------------
%   Compute the calendar date from the Julian date. Uses the format
%   from clock. If no inputs are given it will output the current
%   date and time of the function call.
%-------------------------------------------------------------------------------
%   Form:
%   datetime = JD2Date( jd, structOut )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   jd           (1,1)  Julian date
%   structOut    (1,1)  If entered, output a structure
%
%   -------
%   Outputs
%   -------
%   datetime     (1,6)  [year month day hour minute seconds]
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   References: Montenbruck, O., T.Pfleger, Astronomy on the Personal
%               Computer, Springer-Verlag, Berlin, 1991, p. 13.
%-------------------------------------------------------------------------------
%   Copyright (c) 1993 Princeton Satellite Systems, Inc.
%   All rights reserved.
%-------------------------------------------------------------------------------

if( nargin < 1 )
  jd = [];
end

if( isempty(jd) )
  datetime = clock;
else

  seconds = (jd-R2P5(jd))*86400;

  jd0     = fix(jd+0.5);

  if( jd0 < 2299161 ) % Gregorian calendar
    c = jd0; 
  else
    b = fix(((jd0-1867216) - 0.25)/36524.25);
    c = jd0 + b - fix(b/4) + 1;
  end
  c = c + 1524;

  d = fix((c-122.1)/365.25); 
  e = 365*d + fix(d/4); 
  f = fix((c-e)/30.6001); 

  datetime(2) = f-1-12*fix(f/14) ;
  datetime(1) = d-4715-fix((7 + datetime(2))/10); 
  datetime(3) = fix(c-e+0.5)-fix(30.6001*f);

  datetime(4) = fix(seconds/3600);
  seconds     = seconds - 3600*datetime(4);
  datetime(5) = fix(seconds/60);
  datetime(6) = seconds - 60*datetime(5);

  if ( datetime(1) <= 0 ),
    datetime(1) = datetime(1)-1;
  end
end

if( nargin > 1 )
  datetime = DTAToDTS( datetime );
end

% PSS internal file version information
%--------------------------------------
% $Date: 2004-10-27 13:36:19 -0400 (Wed, 27 Oct 2004) $
% $Revision: 4251 $
