function x = R2P5( y )

%-------------------------------------------------------------------------------
%   Rounds towards zero to the nearest 1/2
% INPUTS:
%   y           Number
%
% OUTPUTS:
%   x           Rounded number
	
x = fix(y); 

d = y - x;

i   = find( d >= 0.5 );
j   = find( d <  0.5 );
x(i) = x(i) + 0.5;
x(j) = x(j) - 0.5;


