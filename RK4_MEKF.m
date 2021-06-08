% Runge-Kutta 4th order integration of covariance matrix and quaternion
% derivatives to obtain propagated (a priori) values of P and q for next
% step

% INPUTS:
% omega_est:
% gyro_ARW:
% P_k:
% quat_k:
% dt:

% OUTPUTS:
%P_k1_pred: propagated P matrix (6x6)
%q_k1_pred: propagated quaternion (scalar component is first element) (4x1)


function[P_k1_pred, q_k1_pred] = RK4_MEKF(omega_est,gyro_ARW,P_k,quat_k,dt)

F = [-CrossMatrix(omega_est)    -eye(3);
         zeros(3)             zeros(3)];
     
G = [-eye(3)  zeros(3);
     zeros(3)  eye(3)];
 
Q = [(gyro_ARW^2).*eye(3)     zeros(3);
         zeros(3)          (gyro_ARW^2).*eye(3)];

k1 = F*P_k + P_k*F' + G*Q*G';
k2 = F*(P_k+0.5*dt*k1) + (P_k+0.5*dt*k1)*F' + G*Q*G';
k3 = F*(P_k+0.5*dt*k2) + (P_k+0.5*dt*k2)*F' + G*Q*G';
k4 = F*(P_k+0.5*dt*k3) + (P_k+0.5*dt*k3)*F' + G*Q*G';

P_k1_pred = P_k + dt*((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4);


k1 = 0.5*Xi_matrix(quat_k)*omega_est';
k2 = 0.5*Xi_matrix(quat_k + 0.5*dt*k1)*omega_est';
k3 = 0.5*Xi_matrix(quat_k + 0.5*dt*k2)*omega_est';
k4 = 0.5*Xi_matrix(quat_k + 0.5*dt*k3)*omega_est';


quat_k1_pred = quat_k + dt*((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4); %scalar component is last
quat_k1_pred = quat_k1_pred/norm(quat_k1_pred);

q_k1_pred = [quat_k1_pred(4); quat_k1_pred(1:3)];   %scalar component is first

end