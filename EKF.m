

%inputs:


%outputs: q_k1: filtered quaternion for next step (1 x 4)
%         P_k1: covariance matrix for next step   (3 x 3)
function[q_k1,P_k1]  =  EKF(gyro_output,prev_quat,P_k0,gyro_ARW,FSS_stdev, MAG_stdev, sun_eci, mag_eci, FSS_output, MAG_output, dt)

q_k0        = prev_quat';  %making col vector
% converting degrees to rad
gyro_output =  gyro_output*(pi/180);
gyro_ARW = gyro_ARW*(pi/180);
FSS_stdev = FSS_stdev*(pi/180);
% normalizing magnetic field vector
mag_eci     = mag_eci./norm(mag_eci);
MAG_output  = MAG_output./norm(MAG_output);

%omega is now in rad/s
wx   = gyro_output(1);
wy   = gyro_output(2);
wz   = gyro_output(3);
w    = norm(gyro_output);

theta = 0.5*norm(gyro_output)*dt;
s     = sin(theta);
c     = cos(theta);

O     = CrossMatrix(gyro_output);

ExpMatrix = [    c     (-wx/w)*s  (-wy/w)*s  (-wz/w)*s;
              (wx/w)*s      c      (wz/w)*s   (-wy/w)*s;
              (wy/w)*s  (-wz/w)*s     c        (wx/w)*s;
              (wz/w)*s   (wy/w)*s  (-wx/w)*s      c    ];

% Projection for next step          
q_k1_pred = ExpMatrix*q_k0; 
q_k1_pred = (q_k1_pred./norm(q_k1_pred))';

%Pdot_k1_pred = -O*P_k0 + P_k0*((-O)') + (gyro_ARW^2)*eye(3,3);
P_k1_pred    = RK4_covariance(gyro_output,gyro_ARW,P_k0,dt);    %projecting covariance matrix

if isnan(norm(FSS_output)) % if sun vector reading is unavailable, reduce the formulation 
   H_k1 = CrossMatrix(quatrotate(q_k1_pred,mag_eci')); 
   R = (MAG_stdev^2).*eye(3,3);
   K_k1 = P_k1_pred*(H_k1')*inv(H_k1*P_k1_pred*(H_k1') + R);
   estimate_vector = quatrotate(q_k1_pred,mag_eci')';
   difference = (MAG_output - estimate_vector);
   correction   = K_k1*difference;

else
    % Both H_k1 formulation below work. Which is correct?
    H_k1 = [ CrossMatrix(quatrotate(q_k1_pred,sun_eci'));
             CrossMatrix(quatrotate(q_k1_pred,mag_eci'))];

    % H_k1 = [ CrossMatrix(quat2rotm(q_k1_pred)*sun_eci);
    %          CrossMatrix(quat2rotm(q_k1_pred)*mag_eci)];

    % H_k1 = [ quat2rotm(q_k1_pred)*CrossMatrix(sun_eci);
    %          quat2rotm(q_k1_pred)*CrossMatrix(mag_eci)];

    %Measurement noise covariance
    R = [(FSS_stdev^2).*eye(3,3)       zeros(3,3)  ;
             zeros(3,3)          (MAG_stdev^2).*eye(3,3)];
    %Compute Kalman Gain
    K_k1         = P_k1_pred*(H_k1')*inv(H_k1*P_k1_pred*(H_k1') + R);
    % estimate_vector = [quat2rotm(q_k1_pred)*sun_eci; 
    %                   quat2rotm(q_k1_pred)*mag_eci];
     %This gave a lot of error! Getting true body vector using quatrotate
     %instead of DCM worked!
    estimate_vector = [quatrotate(q_k1_pred,sun_eci')';
                      quatrotate(q_k1_pred,mag_eci')'];
    y_k1  = [ FSS_output;      %measurement vector 
              MAG_output];
    difference = (y_k1 - estimate_vector);
    %Error quaternion
    correction   = K_k1*difference;
end

% Update quaternion and covariance matrix
q_k1         = quatmultiply(q_k1_pred,[1; (correction./2)]');
q_k1         = q_k1./norm(q_k1);
P_k1         = (eye(3,3) - K_k1*H_k1)*P_k1_pred;

end
