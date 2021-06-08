
% Quaternion Conjugation Function 

function conj = Q_Conj(q)
conj = [q(1) -q(2) -q(3) -q(4)].';  %Conjugates quaternion in column form
end