
% 3x3 Cross Product Operator for a 4x1 Vector

function cross = Cross3x34(v)
cross = [v(1) -v(4) v(3); v(4) v(1) -v(2); -v(3) v(2) v(1)];
end