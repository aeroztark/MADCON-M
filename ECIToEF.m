function c = ECIToEF( T )
	
%-------------------------------------------------------------------------------
%   Computes the matrix from mean of Aries 2000 to the earth fixed frame.
%   Accounts for Earth rotation, nutation, and precession.
%-------------------------------------------------------------------------------
%   Form:
%   c = ECIToEF( T )	
%-------------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   T                 Julian century
%
%   -------
%   Outputs
%   -------
%   c		          Matrix from Mean of 2000.0 to planet fixed
%
%-------------------------------------------------------------------------------
%   References:   Seidelmann, The Explanatory Supplement to the Astronomical
%                 Almanac, p. 705.
%                 Montenbruck, O., Pfleger, T., "Astronomy on the Personal
%                 Computer, Second Edition."
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%   Copyright 1999 Princeton Satellite Systems, Inc.
%   All rights reserved.
%-------------------------------------------------------------------------------

% Input processing
%-----------------
if( nargin < 1 )
  T = [];
end

% Defaults
%---------
if( isempty(T) )
  T = JD2T(Date2JD);
end

jD         = T2JD( T );
jDStandard = 2451545.0;
d          = jD - jDStandard;

c          = MECIToPlanet( -0.641*T, 90 - 0.557*T, 190.16 + 360.9856235*d ); % 190.16 from Ref 

%-------------------------------------------------------------------------------
%   References:   Montenbruck, O., Pfleger, T., "Astronomy on the Personal
%                 Computer, Second Edition."
%-------------------------------------------------------------------------------

function m = MECIToPlanet( alpha0, delta0, w )

degToRad = pi/180;
cW       = cos( w     *degToRad );
sW       = sin( w     *degToRad );
cA       = cos( alpha0*degToRad );
sA       = sin( alpha0*degToRad );
cD       = cos( delta0*degToRad );
sD       = sin( delta0*degToRad );

sWsA     = sW*sA;
sWcA     = sW*cA;
cWsA     = cW*sA;
cWcA     = cW*cA;

m        = [-cWsA - sWcA*sD   cWcA - sWsA*sD  sW*cD;...
             sWsA - cWcA*sD  -sWcA - cWsA*sD  cW*cD;...
		          cA*cD           sA*cD      sD];


