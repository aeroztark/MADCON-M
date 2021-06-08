function[P_k1_pred] = RK4_covariance(omega,gyro_ARW,P_k0,dt)

k1 = -CrossMatrix(omega)*P_k0 + P_k0*(-CrossMatrix(omega))' + (gyro_ARW^2).*eye(3,3);
k2 = -CrossMatrix(omega)*(P_k0+0.5*dt*k1) + (P_k0+0.5*dt*k1)*(-CrossMatrix(omega))' + (gyro_ARW^2).*eye(3,3);
k3 = -CrossMatrix(omega)*(P_k0+0.5*dt*k2) + (P_k0+0.5*dt*k2)*(-CrossMatrix(omega))' + (gyro_ARW^2).*eye(3,3);
k4 = -CrossMatrix(omega)*(P_k0+0.5*dt*k3) + (P_k0+0.5*dt*k3)*(-CrossMatrix(omega))' + (gyro_ARW^2).*eye(3,3);

P_k1_pred = P_k0 + dt*((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4);