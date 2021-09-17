function d = Dot ( w, y )

%-------------------------------------------------------------------------------
%   Dot product. The number of columns of w and y can be:
%   Both > 1 and equal
%   One can have one column and the other any number of columns
%-------------------------------------------------------------------------------
%   Form:
%   d = Dot ( w, y )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   w                 (:,:)  Vector
%   y                 (:,:)  Vector
%
%   -------
%   Outputs
%   -------
%   d                 (1,:)   Dot product of w and y
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%	 Copyright 1993-1997 Princeton Satellite Systems, Inc. All rights reserved.
%-------------------------------------------------------------------------------

if( nargin < 2 )
  y = w;
end

[rW,cW] = size(w);
[rY,cY] = size(y);

if( cW == cY )
  d = sum(w.*y); 
		 
elseif( cW == 1)
  d = w'*y; 
		 
elseif( cY == 1)
  d = y'*w;
   
else
  error('w and y cannot have different numbers of columns unless one has only one column')
end


%--------------------------------------
% PSS internal file version information
%--------------------------------------
% $Date: 2006-09-19 10:24:58 -0400 (Tue, 19 Sep 2006) $
% $Revision: 6945 $
