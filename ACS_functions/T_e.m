
% The T_e function from Oland PD+

function t_e = T_e(q)                        % 4x1 current quaternion input
eeta = q(1,1);
eps = q(2:end,1);
t_e = 0.5*[eps.'; eeta*eye(3)+Cross3x33(eps)];   % Taken from Oland
end