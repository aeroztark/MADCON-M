    
% Rate of Change of Quaternion (Vector Part)

function q_dot = Q_Dot_Vector(q,omega)              % 4x1 quaternion input, 3x1 angular velocity input
q_dot = 0.5*Cross3x34(q)*omega;                     % Taken from Yang
end