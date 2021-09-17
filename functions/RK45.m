% function [x, hLast] = RK45( Fun, x, h, hMax, hMin, tol, t, varargin )
function [tout, x] = RK45( Fun, tspan, x, tol, varargin )


%-------------------------------------------------------------------------------
%   Fourth/fifth order Runge-Kutta. 
%   Called function is of the form:
%
%   Fun(x,t,p1,p2,...)
%
%   Accepts optional arguments that are passed through to Fun.
%   Time is also optional.
%
%   This function will integrate Fun from the current t to t + hMax.
% 
%-------------------------------------------------------------------------------
%   Form:
%   [x, hLast] = RK45( Fun, x, h, hMax, hMin, tol, t, varargin )
%-------------------------------------------------------------------------------
%   ------
%   Inputs
%   ------
%   Fun                Function    Fun(x,t,p1,p2...)
%   x                  State (column vector)
%   h                  Independent variable step
%   hMax               Maximum step size
%   hMin               Minimum step size
%   tol                Tolerance on error
%   t                  Current time
%   varargin           Optional arguments
%
%   -------
%   Outputs
%   -------
%   x                  Updated state
%   hLast              Independent variable step
%
%-------------------------------------------------------------------------------
%	Reference: Cash, J.R., A. H. Karp, "A Variable Order Runge-Kutta Method for 
%              Initial Value Problems with Rapidly Varying Right-Hand Sides,"
%              ACM Trans. on Math. Soft., Vol. 16, No.3, Sept. 1990, pp 201-222.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   Copyright 1996 Princeton Satellite Systems, Inc.
%   All rights reserved.
%-------------------------------------------------------------------------------
if length(tspan)>2
    tspan = [tspan(1) tspan(end)];
end

hMax = tspan(2) - tspan(1);
h = hMax;
hMin = 0.000001;
t = tspan(1);

% if( nargin < 7 )
%   t = [];
% end

% Constants
%----------
scaleFactor = 0.9;
power       = 0.2;
maxErrCount = 10;

% Coefficients
%-------------
c   = [0.2 0.3 0.6 1.0 0.875];
a   = [       0.2       0         0            0        0;...
             3/40    9/40         0            0        0;...
	          0.3    -0.9       1.2            0        0;...
	       -11/54     2.5    -70/27        35/27        0;...
	   1631/55296 175/512 575/13824 44275/110592 253/4096];
	   
b   = [1.0 -1.5 2.5 19/54 0 -10/27 55/54 2825/27648 0 18575/48384 13525/55296 277/14336 0.25 37/378 0 250/621 125/594 0 512/1771];
	
del = [b(14)-b(8),0,b(16)-b(10),b(17)-b(11),-b(12),b(19)-b(13)];	
	
% Start of step with no accumulated step size
%--------------------------------------------
hAccum   = 0;
errCount = 0;

while( hAccum < hMax )
	
  % Adjust the step size
  %---------------------
  hLast = h;
  h     = min(hMax - hAccum,h);
	
  % Call the right hand the first time
  %-----------------------------------
  if ~isempty(t)
    k1  = feval(Fun, x, t, varargin{:});   
    k2  = feval(Fun, x+h*a(1,1)*k1, t+h*c(1), varargin{:});   
    k3  = feval(Fun, x+h*(a(2,1)*k1+a(2,2)*k2), t+h*c(2), varargin{:});   
    k4  = feval(Fun, x+h*(a(3,1)*k1+a(3,2)*k2+a(3,3)*k3), t+h*c(3), varargin{:});   
    k5  = feval(Fun, x+h*(a(4,1)*k1+a(4,2)*k2+a(4,3)*k3+a(4,4)*k4), t+h*c(4), varargin{:});   
    k6  = feval(Fun, x+h*(a(5,1)*k1+a(5,2)*k2+a(5,3)*k3+a(5,4)*k4 + a(5,5)*k5), t+h*c(5), varargin{:});   
  else
    k1  = feval(Fun, x, varargin{:});   
    k2  = feval(Fun, x+h*a(1,1)*k1, varargin{:});   
    k3  = feval(Fun, x+h*(a(2,1)*k1+a(2,2)*k2), varargin{:});   
    k4  = feval(Fun, x+h*(a(3,1)*k1+a(3,2)*k2+a(3,3)*k3), varargin{:});   
    k5  = feval(Fun, x+h*(a(4,1)*k1+a(4,2)*k2+a(4,3)*k3+a(4,4)*k4), varargin{:});   
    k6  = feval(Fun, x+h*(a(5,1)*k1+a(5,2)*k2+a(5,3)*k3+a(5,4)*k4 + a(5,5)*k5), varargin{:});   
  end
    
  % Find the difference between the 4th and 5th order solutions
  %------------------------------------------------------------
  errorN = norm(h*abs(del(1)*k1 + del(3)*k3 + del(4)*k4 + del(5)*k5 + del(6)*k6),inf);
	
  % Use relative error
  %-------------------
  tau = tol*max(norm(x,inf),1);
	
  % Update the equations if the error is acceptable
  %------------------------------------------------
  if( errorN < tau )
    hAccum   = hAccum + h;
    x        = (x + h*(b(14)*k1 + b(16)*k3 + b(17)*k4 + b(19)*k6));
    errCount = 0;
  else
	errCount = errCount + 1;
	if( errCount > maxErrCount)
	  error(sprintf('h = %12.4e, hMax = %12.4e, error = %12.4e, tau = %12.4e',h,hMax,errorN,tau));
	end
  end
	  
  % Update h and limit it to between hMin and hMax
  %-----------------------------------------------
  if( errorN ~= 0 )
	h = scaleFactor * h * ( tau / errorN )^power;
    if( h > hMax )
	  h = hMax;
	elseif ( h < hMin )
	  h = hMin;
    end	   
  end
  
  tout = t + hMax;
  
end

x = x';

%-------------------------------------------------------------------------------
% PSS internal file version information
%--------------------------------------
% $Date: 2007-08-01 23:58:00 -0400 (Wed, 01 Aug 2007) $
% $Revision: 10780 $
