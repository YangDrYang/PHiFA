function [tout, x] = RK4( Fun, tspan, x, varargin )

%-------------------------------------------------------------------------------
%   Fourth order Runge-Kutta. Called function is of the form:
%
%   Fun(x,t,varargin)
%
%   Accepts optional arguments that are passed through to Fun. Time is also
%   optional.
%   This version is streamlined and designed to take advantage of MATLAB's
%   function handles feature (MATLAB versions 6 and up only). Passing a 
%   string function name (version 5) will also work but is slower.
% 
%-------------------------------------------------------------------------------
%   Form:
%   x = RK4( Fun, x, h, t, varargin )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   Fun                Function    Fun(x,{t,...})
%   x                  State (column vector)
%   h                  Independent variable step
%   t                  Current time
%   p1...              Optional arguments
%
%   -------
%   Outputs
%   -------
%   x                  Updated state
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   Copyright 1993-1994, 2006 Princeton Satellite Systems, Inc.
%   All rights reserved.
%-------------------------------------------------------------------------------
if length(tspan)>2
    tspan = [tspan(1) tspan(end)];
end

h = tspan(2) - tspan(1);
t = tspan(1);
tout = tspan(2);

% if( nargin < 4 )
%   t = [];
% end

ho2 = 0.5*h;

if ~isempty(t)
  k1  = feval( Fun, x, t, varargin{:} );
  k2  = feval( Fun, x + ho2*k1, t+ho2, varargin{:} );
  k3  = feval( Fun, x + ho2*k2, t+ho2, varargin{:} );
  k4  = feval( Fun, x + h*k3, t+h, varargin{:} );
else
  k1  = feval( Fun, x, varargin{:} );
  k2  = feval( Fun, x + ho2*k1, varargin{:} );
  k3  = feval( Fun, x + ho2*k2, varargin{:} );
  k4  = feval( Fun, x + h*k3, varargin{:} );
end

x   = (x + h*(k1 + 2*(k2+k3) + k4)/6)';


%--------------------------------------
% PSS internal file version information
%--------------------------------------
% $Date: 2007-07-26 09:10:10 -0400 (Thu, 26 Jul 2007) $
% $Revision: 10679 $
