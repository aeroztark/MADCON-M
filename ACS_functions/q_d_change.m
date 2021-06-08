
% Change of Desired Quaternion

function q_d_new = q_d_change(q_d,q_d_dot,dt)             % 4x1 quaternion input, 4x1 quaternion output
q_d_new = Q_Normalise(q_d+q_d_dot*dt);             %
end