
function[RW_torque, MTR_torque] = ACS(ACSmode,q,q_d,omega, omega_d,omega_dot_d,gains)
% function to simulate Attitude Control algorithm

% INPUTS: ACSmode (string) - specify the ACS operational mode
           % None, Idle, SunPointing, TargetPointing, NadirPointing
           % q -> (4x1) - current attitude quaternion
           % q_d -> (4x1)
% OUTPUTS: RW_torque (3x1 vector) - sum of all reaction wheel torques
%          MTR_torque (3x1 vector) - sum of torques from all torque rods

addpath('D:\Documents\GitHub\MADCON-M\ACS_functions');  % library of AC specific functions for better organization

%% Initialize ACS



% if ~(strcmp(ACSmode,'None')) % run detumbling & desaturation whenever ACSmode is not None
%% Detumbling check

%     %B-dot control for detumbling
%  %if sat_mode(count) ~= 0
%      
%     if sc.is_detumble == 0 || norm(omega(count,:)) < (pi/360)    %do not detumble if rate < 0.5 deg/s or if flag is set to 0
%         torque_B = [0;0;0];   %do not detumble
%     else
%         if is_MAG_on == 1   %keep torquer off when reading the MAG
%             torque_B = [0;0;0];
%             is_MAG_on = 0;  %flip the flag to enable MTR actuation in next instance.
%         else                %i.e. if is_MAG_on == 0, then actuate 
%             control_moment(count,:) = Bdot_moment(vectors.mag_b_noisy(count,:)',(omega(count,:))',sc.BDot_gain,sc.max_dipole);
%             torque_B = (cross((control_moment(count,:))',vectors.mag_b_true(count,:)'));
%             is_MAG_on = 1;  %flip the flag to enable MAG reading in next instance.
%         end
%     end
%  %end
   

%% RW desaturation check

% end

%% Mission modes

% if strcmp(ACSmode,'Idle') || strcmp(ACSmode,'None')  % no control torques are applied
    RW_torque = [0;0;0];
    MTR_torque = [0;0;0];
% end

% if strcmp(ACSmode,'SunPointing')
% end

% if strcmp(ACSmode,'TargetPointing')
    % Error computations
    q_error = Q_Error(q,q_d); % Attitude error quaternion (calculations based on Oland PD+ paper)
    omega_error = omega_d - omega;              
   % omega_dot_error = omega_dot_d - omega_dot;  
    
    if isequal(omega_d,[0; 0; 0])
    % The control law PD (based on Sarthak's pic, it works better, don't know why)
        RW_torque = gains.K_p*q_error(2:end)*q_error(1) + gains.K_d*omega_error;
    else
    % Alt PD control law for non-zero desired omega
        RW_torque = gains.K_p_w*omega_error + gains.K_d_w*omega_dot_error;
    end
    
    MTR_torque = [0;0;0];
%     % Decomposition into RWs
%     [H,omega_rw_new] = reaction_wheels(J_RW,rw_align,omega_rw,T,dt);
%     voltage_rw_new = omega_rw_new/rw_conversion;
% 
%     % Runge-Kutta Fourth Order numerical integration, output q normalisation
%     [omega_new,q_new,omega_dot_new] = RK4(J,J_inv,H,q,omega,dt,T);
%     q_new = Q_Normalise(q_new);
    
    
% end

% if strcmp(ACSmode,'NadirPointing')
% end
    