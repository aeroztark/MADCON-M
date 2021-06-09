
function[new_w,new_q] = RK4_v2(I,w,torque,q,dt)

% Function to propagate attitude by integrating dynamics & kinematics equations through RK4 

%   Implementing RK4 for updating angular velocity
    k1 = w_dot(I,w,torque);
    k2 = w_dot(I,w+0.5*k1(1:3)*dt,torque);
    k3 = w_dot(I,w+0.5*k2(1:3)*dt,torque);
    k4 = w_dot(I,w+k3(1:3)*dt,torque);

    new_w = w + (1/6)*dt*(k1+2*k2+2*k3+k4);

%   Implementing RK4 for updating quaternion. The current updated value of
%   angular velocity (new_w) should be used and not old timestep omega (w).
%   The latter can be used but the timsetep has to be 0.1s or less.
    k1 = q_dot(new_w,q);
    k2 = q_dot(new_w,q+0.5*k1*dt);
    k3 = q_dot(new_w,q+0.5*k2*dt);
    k4 = q_dot(new_w,q+k3*dt);

    new_quat = q + (1/6)*dt*(k1+2*k2+2*k3+k4);
    new_q = new_quat/norm(new_quat);
    q = new_q;     

        %the ODEs are written as separate functions for neatness
        function [w_dot] = w_dot(I,w,torque) 
            w_dot = inv(I)*(torque - cross(w,I*w));   %Dynamics equation
        end
        
        function [q_dot] = q_dot(w,q)         %quaternion kinematic equation
            Q_matrix = [q(1) -q(2) -q(3) -q(4);
                        q(2)  q(1) -q(4)  q(3);
                        q(3)  q(4)  q(1) -q(2);
                        q(4) -q(3)  q(2)  q(1)];

            wk = [0;w(1);w(2);w(3)];
            q_dot = 0.5*Q_matrix*wk;
        
%           Alternatively, the other form of quaternion differential equation
%           could be used.
        
%           W_matrix = [   0      -w(1) -w(2) -w(3);         
%                         w(1)      0    w(3) -w(2);
%                         w(2)    -w(3)    0   w(1);
%                         w(3)     w(2) -w(1)   0];
%           q_dot = 0.5*W_matrix*q;        
        end

end
   