function [n, eclipseType] = Eclipse( rSC, rSource, rPlanet, radiusPlanet, radiusSource )

%-------------------------------------------------------------------------------
%   Computes eclipses.
%-------------------------------------------------------------------------------
%   Form:
%   [n, eclipseType] = Eclipse( rSC, rSource, rPlanet, radiusPlanet, radiusSource )
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   rSC             (3,:)      Position vector of the spacecraft
%   rSource         (3,: or 1) Position vector of the light source
%   rPlanet         (3,: or 1) Position vector of the planet causing the eclipse
%                              (default is zeros)
%   radiusPlanet    (1,1)      Radius of the planet causing the eclipse 
%                              (default is the earth)
%   radiusSource    (1,1)      Radius of the light source (default is the sun)
%
%   -------
%   Outputs
%   -------
%   n               (1,:)      Normalized Source intensity ( 0 < n < 1 )
%   eclipseType     (1,:)      0 = none 1 = partial, 2 = annular 3 = total
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   Copyright (c) 1995-1996, 2004 Princeton Satellite Systems, Inc.
%   All rights reserved.
%-------------------------------------------------------------------------------

% Input processing
%-----------------
if( nargin < 5 )
  radiusSource =  [];
end

if( nargin < 4 )
  radiusPlanet =  [];
end

if( nargin < 3 )
  rPlanet =  [];
end

% Defaults
%---------
if( isempty( radiusSource ) )
  radiusSource =  6.9499e5;
end

if( isempty( radiusPlanet )  )
  radiusPlanet =  6378.165;
end

mSC = size(rSC,2);
mS  = size(rSource,2);

if( mS == 1 )
  rSource = DupVect( rSource, mSC );
end

if( isempty( rPlanet) )
  dPlanet = -rSC;
  rPlanet = zeros(3,size(rSource,2));
else
  mP = size(rPlanet,2);
  if( mP == 1 )
    rPlanet = DupVect(rPlanet,mSC);
  end
  dPlanet = rPlanet - rSC;
end

dSource      = rSource - rSC;
nDS          = Mag(dSource);
nDP          = Mag(dPlanet);

S            = Mag(rSource - rPlanet);  % Distance between the source and planet
 
if( radiusSource == radiusPlanet )
  C = inf;
else
  C = (radiusPlanet/(radiusSource-radiusPlanet))*S; 
end

sinRhoPlanet = radiusPlanet./nDP;    
sinRhoSource = radiusSource./nDS; 

rhoPlanet    = asin(sinRhoPlanet); 
rhoSource    = asin(sinRhoSource);

cosRhoPlanet = cos(rhoPlanet); 
cosRhoSource = cos(rhoSource); 

cosTheta     = Dot(dPlanet,dSource)./(nDP.*nDS);   
theta        = acos(cosTheta);   
sinTheta     = sin(theta);  

deltaRho     = rhoPlanet - rhoSource;

% No Eclipse
%-----------  
m            = length(sinTheta);
n            = ones(1,m);
eclipseType  = zeros(1,m);

% Partial Eclipse
%----------------  
k  = find( (nDS > S) & (rhoSource+rhoPlanet > theta) & (theta > abs(deltaRho)) );
lK = length(k);
if( lK > 0 )
  f1             = (cosRhoPlanet(k) - cosRhoSource(k).*cosTheta(k)    )./(sinRhoSource(k).*sinTheta(k)); 
  f2             = (cosRhoSource(k) - cosRhoPlanet(k).*cosTheta(k)    )./(sinRhoPlanet(k).*sinTheta(k)); 
  f3             = (cosTheta(k)     - cosRhoSource(k).*cosRhoPlanet(k))./(sinRhoSource(k).*sinRhoPlanet(k)); 

  a              = pi - cosRhoSource(k).*acos(f1) - cosRhoPlanet(k).*acos(f2) - acos(f3);

  n(k)           = 1 - a./ ( pi*(1-cosRhoSource(k)) ); 
  eclipseType(k) = ones(1,lK);
end

% Annular Eclipse
%----------------  
k  = find( (abs(theta+rhoPlanet) < rhoSource) ); 
lK = length(k);
if( length(k) > 0 )
  n(k)           = 1 - (1-cosRhoPlanet(k))./(1-cosRhoSource(k));
  eclipseType(k) = DupVect(2,lK);
end
  
% Total Eclipse
%----------------  
k  = find( (deltaRho > theta) & (nDS > S) & ((S+C) > nDS) );
lK = length(k);
if( length(k) > 0 )
  n(k)           = zeros(1,lK);
  eclipseType(k) = DupVect(3,lK); 
end

% Plot if no output arguments
%----------------------------
if( nargout == 0 )
  Plot2D(1:m,n,'Sample','Eclipse Fraction');
  yLabels = [...
             'None   ';...
			 'Partial';...
			 'Annular';...
             'Total  ';...
			];
  NPlot( yLabels, eclipseType+1, 1:m, 'Sample', 'Eclipse Type', 'Eclipse Type', 'Eclipse Type' );
  clear n
end

% PSS internal file version information
%--------------------------------------
% $Date: 2007-03-30 09:48:44 -0400 (Fri, 30 Mar 2007) $
% $Revision: 9036 $
