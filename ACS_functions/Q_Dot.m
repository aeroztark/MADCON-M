
% Rate of Change of Quaternion

function q_dot = Q_Dot(q,omega)             % 4x1 quaternion input, 3x1 angular velocity input
q_dot = 0.5*Cross4x43(omega)*q;             % Taken from Yang
end