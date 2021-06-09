function T = JD2T( jd )
	
%-------------------------------------------------------------------------------
%   Converts Julian days to centuries from J2000.0

% INPUTS:
%   jd           Julian date (days)
%
% OUTPUTS:
%   T            Julian centuries of 86400s dynamical time from j2000.0


if( nargin < 1 )
  jd = Date2JD;
end

T = (jd - 2451545) / 36525;

