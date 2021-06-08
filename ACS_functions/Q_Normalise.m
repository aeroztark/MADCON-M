
% Quaternion Normalisation Function

% This function normalises quaternions in column form
function norm = Q_Normalise(q)
norm = q/Q_Norm(q);
end