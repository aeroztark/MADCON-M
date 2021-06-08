function y = DupVect( x, n )

%---------------------------------------------------------------------------
%   Create a matrix with n rows or columns each of which equals the 
%   row or column vector x. When duplicating a scalar note that
%   For example,
%   DupVect(3,5) = [3 3 3 3 3]'
%---------------------------------------------------------------------------
%   Form:
%   y = DupVect( x, n )
%---------------------------------------------------------------------------
%
%   ------
%   Inputs
%   ------
%   x                  Vector to be duplicated
%
%   -------
%   Outputs
%   -------
%   y                  Matrix with n rows or n columns of x
%
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
%    Copyright 1995 Princeton Satellite Systems, Inc.
%     All rights reserved.
%---------------------------------------------------------------------------

if( n < 1 )
  error('n must be greater than 0')
end

[r,c] = size(x);

if( r > c )
  y = x(:,ones(1,n));
else
  y = x(ones(n,1),:);
end 

