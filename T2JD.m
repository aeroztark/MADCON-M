function jd = T2JD( T )
	
%-------------------------------------------------------------------------------
%   Converts Julian centuries from J2000.0 to days
% INPUTS:
%   T      (1,:)  Julian centuries of 86400s dynamical time from j2000.0
%
% OUTPUTS:
%   jd     (1,:)  Julian date (days)


jd = T*36525 + 2451545;

