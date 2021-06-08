%This function computes the Xi matrix for a quaternion

%input: q as column vector with scalar component as the last row
%output: Xi matrix (source: Optimal Control of Dynamical Systems)


function[Xi] = Xi_matrix(q)
    q_scalar = q(4);
    q_vector = q(1:3);
    
    Xi = [q_scalar*eye(3) + CrossMatrix(q_vector);
                -q_vector'];
end