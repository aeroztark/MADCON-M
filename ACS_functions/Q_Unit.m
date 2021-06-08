
% Quaternin Unit Function

function unit_q = Q_Unit(q)
unit_q = [q(1)/Q_Norm(q) q(2)/Q_Norm(q) q(3)/Q_Norm(q) q(4)/Q_Norm(q)].';
end