function flux = SolarFlx( au )
	
%-------------------------------------------------------------------------------
%   Solar flux from the sun as a function of distance in AU.
%   Inverse quadratic formula scales anywhere in the solar system. The nominal
%   value at 1 AU is 1367 W/m2.
%   Has a built-in demo.
%-------------------------------------------------------------------------------
%   Form:
%   flux = SolarFlx( au )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   au      (1,:) Astronomial units from the sun
%
%   -------
%   Outputs
%   -------
%   flux    (1,:) Solar flux (watts/m^2)
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   Copyright (c) 1993, 2006 Princeton Satellite Systems, Inc.
%   All rights reserved.
%-------------------------------------------------------------------------------

if( nargin == 0 )
   au = linspace(0.3,40);
end

if( nargout == 0 )
  Plot2D(au,1367./au.^2,'Distance (au)','Flux (W/m^2)','Solar Flux','ylog')
else
  flux = 1367./au.^2;
end


%--------------------------------------
% PSS internal file version information
%--------------------------------------
% $Date: 2007-07-26 09:08:53 -0400 (Thu, 26 Jul 2007) $
% $Revision: 10678 $
