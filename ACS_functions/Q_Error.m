    
% Quaternion Error Computation

% Calculated from Kronecker Product
function q_err = Q_Error(q,q_d)             % 4x1 current quaternion input, 4x1 desired quaternion input
q_matrix = [q(1) -q(2) -q(3) -q(4);         % If both inputs are normalised, output is also a normalised q
            q(2)  q(1)  q(4) -q(3);
            q(3) -q(4)  q(1)  q(2);
            q(4)  q(3) -q(2)  q(1)] ;       % From Yang
q_err = q_matrix\q_d;                       % This works, stable outputs, for a while at least              
end