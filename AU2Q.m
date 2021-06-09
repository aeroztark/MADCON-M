function q = AU2Q( angle, u )

%   Converts an angle and a unit vector to a quaternion.
%
% INPUTS:
%   angle         (1,:)   Angle
%   u             (3,1)   Unit vector

halfAngle = angle/2;
c         = cos( halfAngle );
s         = sin( halfAngle );

q         = [ c; -[s*u(1);s*u(2);s*u(3)] ];


