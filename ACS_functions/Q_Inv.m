
% Quaternion Inverse Function

function inv = Q_Inv(q)
inv = Q_Conj(q)/(Q_Norm(q)^2);  %
end