
% 4x4 Cross Product Operator for a 4x1 Vector

function cross = Cross4x44(p)
cross = [p(1) -p(2) -p(3) -p(4); p(2) p(1) -p(4) p(3); p(3) p(4) p(1) -p(2);    % From Yang eq 45a, same as in Q_Product
p(4) -p(3) p(2) p(1)]; 
end