
% 3x3 Cross Product Operator for a 3x1 Vector

function cross = Cross3x33(v)
cross = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]; 
end