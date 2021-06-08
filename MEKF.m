
%Multiplicative quaternion based state error extended Kalman filter 

%inputs:
%         gyro_output:
%            gyro_ARW:
%                b_k0:
%          FSS_output:
%           FSS_stdev:
%          MAG_output:
%           MAG_stdev:
%           prev_quat:
%                P_k0:
%             sun_eci:
%             mag_eci:
%                  dt:

%outputs:       q_k:    filtered (a posteriori)quaternion for current step (4 x 1) with scalar as the first component
%               b_k:    gyro bias 3x1 row vector (a posteriori) for current step
%         omega_est:    estimated omega for current step - row vector in deg/s
%         P_k1_pred:    propagated covariance matrix (a priori) for next step
%         q_k1_pred:    propagated quaternion (a priori) for next step


function[q_k,b_k,omega_est,P_k1_pred,q_k1_pred,error_vector,P_k]  =  MEKF(gyro_output,gyro_ARW,b_k0,FSS_output,FSS_stdev,MAG_output, MAG_stdev,prev_quat,P_k0, sun_eci, mag_eci,SS_avail,MAG_avail,dt)

q_k0        = prev_quat';  %making col vector
gyro_output =  (gyro_output*(pi/180))';    %converting gyro reading to rad/s and making column vector
mag_eci     = mag_eci./norm(mag_eci);       %normalize magnetic vectors
MAG_output  = MAG_output./norm(MAG_output);

y_k  = [ SS_avail.*FSS_output;      %measurement vector 
          MAG_avail.*MAG_output];
      
%% Gain computation      
      H_k  = [ CrossMatrix(quat2rotm(q_k0')*(sun_eci)),    zeros(3);    %6x6 matrix
               CrossMatrix(quat2rotm(q_k0')*(mag_eci)),    zeros(3);];
           
     %Measurement error covariance
      R    = [(FSS_stdev^2).*eye(3,3)        zeros(3,3)  ;              %6x6 matrix
                   zeros(3,3)          (MAG_stdev^2).*eye(3,3)]; 
          
     %Compute Kalman Gain
      K_k  = P_k0*(H_k')*inv(H_k*P_k0*(H_k') + R);  % 6x6 matrix

      
%% State update
      %Update covariance matrix
      P_k  = (eye(6) - K_k*H_k)*P_k0;  
      
      %estimate output vector
      estimate_vector = [SS_avail.*quatrotate(q_k0',sun_eci')';
                         MAG_avail.*quatrotate(q_k0',mag_eci')'];
              
      %Error state vector update
      delta_x   = K_k*(y_k - estimate_vector);
      
      error_vector = delta_x';  %only for reporting purposes
      
      %convert quaternion notation [1 0 0 0] to [0 0 0 1] (q_k0 is a column vector)
      quat_k0 = [q_k0(2:4); q_k0(1)];   %also a column vector
      %This conversion is needed to follow reference text's notation for quaternion update   
      
      % A posteriori updates:
        %Updated quaternion
      quat_k = quat_k0 + 0.5*Xi_matrix(quat_k0)*delta_x(1:3);
      quat_k = quat_k./norm(quat_k);
      
      q_k = [quat_k(4);quat_k(1:3)]';  %back to simulator's format
      
        %Updated gyro bias vector (3x1)
      b_k = b_k0' + delta_x(4:6);
      
        %Updated angular rate
      omega_est = gyro_output - b_k';    %in rad/s
      
      
%% State Propagation

%A priori estimates for next timestep based on propagation by RK4
%integration
[P_k1_pred,q_k1_pred] = RK4_MEKF(omega_est,gyro_ARW,P_k,quat_k,dt);



omega_est = (rad2deg(omega_est));    %convert to deg/s and make row vector

end
