function [q] = QUEST(v1b, v1eci, v2b, v2eci, w1, w2)

%-------------------------------------------------------------------------- 
%       This program simulates the QUEST algorithm
%
%       Created by: Sarthak Srivastava (6-JUN-2018)
%-------------------------------------------------------------------------
%   Inputs: v1b, v1eci, v2b, v2eci vectors in column form and scalar weights w1,w2
%   Output: q (attitude quaternion)
%--------------------------------------------------------------------------

%   convert all vectors to unit vectors
v1b = v1b/norm(v1b);
v1eci = v1eci/norm(v1eci);
v2b = v2b/norm(v2b);
v2eci = v2eci/norm(v2eci);


B = w1*v1b*v1eci' + w2*v2b*v2eci';
S = B + B';
Z = [(B(2,3)-B(3,2)); (B(3,1)-B(1,3)); (B(1,2)-B(2,1))];

sigma = trace(B);   %sum of diagonals
lambda = w1 + w2;    %sum of weights

p = ((lambda + sigma)*eye(3)-S)\Z;     % avoid computing matrix inverse, use LU decomposition instead

q(1) = 1/sqrt(1 + p'*p);
q(2) = q(1)*p(1);
q(3) = q(1)*p(2);
q(4) = q(1)*p(3);
q = [q(1); q(2); q(3); q(4)];

q = q./norm(q);   %added normalization
end