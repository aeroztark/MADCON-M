
% 4x4 Cross Product Operator for a 3x1 Vector

function cross = Cross4x43(v)
cross = [0 -v(1) -v(2) -v(3); v(1) 0 v(3) -v(2); v(2) -v(3) 0 v(1); v(3) v(2) -v(1) 0]; 
end