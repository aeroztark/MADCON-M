
%----------------------------------------------------------------------------- 
%% ADCS simulation - Main file
%Created by: Sarthak Srivastava

%%
%--------------------------------------------------------------------------
clear
clc

%% Importing stuff from param() file

sc              = param_satellite();      %call the appropriate param file

omega_init      = (pi/180).*sc.omega_init_deg;   %in rad/s
dynamics.omega  = omega_init;
step_in_sec     = sc.timestep;

% Convert initial Orbital elements to intial r_eci and v_eci
[r_init,v_init] = oe2rv(sc.oe);     %everything in m and m/s

%tle_datastruct = ConvertNORAD(sc.tle_file);    %if using TLE with SGP4

%% Initialization and Time properties

progressbar     %initializes progress bar
is_MAG_on             = 1;  %default flag value is true

% Convert to Julian Date
jd_start_time         = juliandate(sc.start_time);
jd_end_time           = juliandate(sc.end_time);
jd_timestep           = (sc.timestep/3600)/24;  %converting sec to days

% Define max index (no of rows) for the array to hold the results
count_max             = floor((jd_end_time - jd_start_time)/jd_timestep) + 1;

% initialize arrays
q                     = zeros(count_max,5);
omega                 = zeros(count_max,3);
omega(1,:)            = omega_init';
OP.r_eci              = zeros(count_max+1,3);
OP.v_eci              = zeros(count_max+1,3);
OP.lat                = zeros(count_max,1);
OP.lon                = zeros(count_max,1);
OP.alt                = zeros(count_max,1);
accel.totalgravity    = zeros(count_max,3);
accel.pointEarth      = zeros(count_max,3);
accel.nonsphEarth     = zeros(count_max,3);
accel.drag            = zeros(count_max,3);
rho                   = zeros(count_max,1);
dynamics.torque_total = zeros(count_max,3);
control_moment        = zeros(count_max,3);
torque_mag            = zeros(count_max,1);
eclipse_flag          = zeros(count_max, 1);
attitude.quest        = zeros(count_max,4);
attitude.q_EKF        = zeros(count_max+1,4);
attitude.quaternion   = zeros(count_max+1,4);
attitude.est_q   = zeros(count_max+1,4);
gyro_bias             = zeros(count_max,3);
omega_est             = zeros(count_max,3);
eclipse_flag          = zeros(count_max,1);
error_vec             = zeros(count_max,6);
P_matrix              = zeros(3,3,count_max);
is_eclipse            = 0;
RW_torque             = zeros(3,count_max);
MTR_torque             = zeros(3,count_max);
% Initialize the first rows
OP.r_eci(1,:)         = r_init';
OP.v_eci(1,:)         = v_init';
OP.alt(1)             = norm(r_init) - 6.378137e6;
rho(1)                = AtmDens2(OP.alt(1)/1000);
sunbvector            = [0,0,0];
P_prev                = eye(3);
gyro_bias(1,:)        = [-2 1 -0.5];
ADCS_init_status      = 0;%zeros(count_max,1);  %set ADCS initialization status to 0

% Global variables
global a_point a_nonsphEarth a_drag a_SRP a_thruster

%% Magic happens here

for count = 1:1:count_max
% Establish correspondence between index count and time (julian date format)    
    OP.jd_time          = jd_start_time + (count-1)*jd_timestep;
    OP.no_orbits(count) = count/sc.period;

% Start attitude calculation     
    q(1,2:5)            = sc.q_init';
    DCM                 = [1 0 0;0 1 0;0 0 1];  %initialize
    
%%  Orbit Propagation    
% Propagate orbit using TLE (provide TLE data struct, calendar time & write to array)
    %this is the format--> [r_eci,v_eci, lon, lat, alt] = OrbitPropagatorSGP4(tle,JD2Date(jd_time));
    %[OP.r_eci(count,:),OP.v_eci(count,:),OP.lon(count,:),OP.lat(count,:),OP.alt(count,:)] = OrbitPropagatorSGP4(tle_datastruct,JD2Date(OP.jd_time));
    
% Propagate Orbit using RK4
    
    [r_updated,v_updated] = RK4_orbit(OP.r_eci(count,:)',OP.v_eci(count,:)',step_in_sec,sc.area_deployed,sc.Cd,sc.m,rho(count),DCM,is_eclipse,sunbvector);
    OP.r_eci(count+1,:)   = r_updated';
    OP.v_eci(count+1,:)   = v_updated';
    
    DCM = LVLH2ECI(OP.r_eci(count,:)',OP.v_eci(count,:)');
    %q(count,2:5) = dcm2quat(DCM);  
    
% Convert ECI to LLA 
    [OP.lat(count),OP.lon(count),OP.alt(count)] = ECI2LLA(OP.r_eci(count,:)',OP.jd_time);
    
% Store the accelerations from different Force models
     accel.pointEarth(count,:)   = a_point';
     accel.nonsphEarth(count,:)  = a_nonsphEarth';
     accel.totalgravity(count,:) = accel.pointEarth(count,:) + accel.nonsphEarth(count,:);
     accel.drag(count,:)         = a_drag';    
     accel.SRP(count,:)          = a_SRP';
     
% Calculate Atmospheric density 
    rho(count+1) = AtmDens2(OP.alt(count)/1000);

%% Sun ECI vector using Sun Model    
% Calc Sun vector in ECI using sun.m and write to array
    vectors.sun_eci(count,:) = sun(OP.jd_time);                  %unit sun ECI vector using Sun Model
    
    %(Fun fact: one way to check if the Sun Model is properly working is to obtain sun vector...
    ...on vernal and autumnal equinox dates. The results should be [1,0,0] and [-1,0,0] respectively.)
    
%% Magnetic Field ECI vector using Matlab IGRF13 function 
    vectors.mag_eci(count,:) = igrfmagm(OP.alt(count),OP.lat(count),OP.lon(count),decyear(JD2Date(OP.jd_time)),13).*10^-9;  %mag eci vector in T
%% Eclipse model    
% Check for eclipse
%     eclipse_type = Eclipse((OP.r_eci(count,:)')/1000,((1.5e8)*vectors.sun_eci(count,:))');
%     eclipse_flag(count) = eclipse_type;
      eclipse_flag(count) = EclipseCheck(vectors.sun_eci(count,:)',OP.r_eci(count,:)/1000');
      is_eclipse = eclipse_flag(count);
%% True body frame vectors    
% Using initial attitude, calculate 'true' (i.e. without noise) initial body sun and mag vectors ...
    ...by transformation of eci vectors (1x3) to body frame
    ...matlab quatrotate function is used which requires inputs in row form
%     sunbvector = DCM*vectors.sun_eci(count,:)';   %unit 3x1 vector (deprecated)
    vectors.sun_b_true(count,:) = quatrotate(q(count,2:5),vectors.sun_eci(count,:));
    sunbvector = vectors.sun_b_true(count,:);
    vectors.mag_b_true(count,:) = (quatrotate(q(count,2:5),vectors.mag_eci(count,:))); %in T  3x1 vector

%% Sensor reading simulation
    if eclipse_flag(count) == 0
        SS_avail = 1;   % fine sun sensor is available
    else
        SS_avail = 0;   % fine sun sensor is not available
    end
    MAG_avail = 1;  % for now, MAG is always available
    
    % implement sensor models to get 'noisy' sensor readings
    vectors.sun_b_noisy(count,:) = FSSmodel(vectors.sun_b_true(count,:),sc.sensors.FSSmax_err,sc.sensors.FSS_SNR,eclipse_flag(count));  %FSS model
    vectors.mag_b_noisy(count,:) = MAGmodel(vectors.mag_b_true(count,:),sc.sensors.MAG_SNR); %MAG model
    vectors.gyro_noisy(count,:)  = Gyromodel(omega(count,:),sc.sensors.Gyro_drift,sc.sensors.Gyro_offset, sc.sensors.Gyro_SNR,step_in_sec);            %Gyro model
 

%% Attitude Estimation
%Initialization
    if ADCS_init_status == 0
        [attitude.est_q(count,:),init_mode] = ADS_init(SS_avail,MAG_avail,vectors.mag_b_noisy(count,:)',vectors.mag_eci(count,:)',vectors.sun_b_noisy(count,:)',vectors.sun_eci(count,:)');
        attitude.est_q(count+1,:) = attitude.est_q(count,:);    %in absence of EKF, the next state is the same as current
        ADCS_init_status = 1;   %initialization has happened
        disp('initialization\n')
    else
%         [attitude.est_q(count,:),gyro_bias(count,:),omega_est(count,:),P_pred,q_pred,error_vec(count,:),P_matrix(:,:,count)] = MEKF(vectors.gyro_noisy(count,:)',sc.sensors.Gyro_ARW,gyro_bias(count,:),vectors.sun_b_noisy(count,:)',sc.sensors.FSS_stdev,vectors.mag_b_noisy(count,:)', sc.sensors.MAG_stdev,attitude.est_q(count,:),P_prev, ...
%                                           vectors.sun_eci(count,:)',vectors.mag_eci(count,:)',SS_avail,MAG_avail,step_in_sec);
                                      
        [attitude.est_q(count,:),P_pred] = EKF(vectors.gyro_noisy(count,:)',attitude.q_EKF(count,:),P_prev,sc.sensors.Gyro_ARW, sc.sensors.FSS_stdev,...
                                sc.sensors.MAG_stdev,vectors.sun_eci(count,:)',vectors.mag_eci(count,:)',vectors.sun_b_noisy(count,:)',...
                                 vectors.mag_b_noisy(count,:)', step_in_sec);
                                      
        P_prev = P_pred;    %propagated P and q
        attitude.est_q(count+1,:) = attitude.est_q(count,:)';
        if (SS_avail && MAG_avail) ~= 1
            init_mode = 0;
        end
    end
    %ADCS_init_status(count+1) = 1;   %initialization has happened
    if (init_mode == 0) && (SS_avail && MAG_avail == 1)
        ADCS_init_status = 0;   %need to re-initialize
         disp('Re-initialization\n')
    end
 
 attitude.quest(count,:) = QUEST(vectors.mag_b_noisy(count,:)',vectors.mag_eci(count,:)',vectors.sun_b_noisy(count,:)',vectors.sun_eci(count,:)',0.4,0.6); 
 attitude.questtruevec(count,:) = QUEST(vectors.mag_b_true(count,:)',vectors.mag_eci(count,:)',vectors.sun_b_true(count,:)',vectors.sun_eci(count,:)',0.4,0.6);    
    if sc.use_EKF == 1
        if (SS_avail && MAG_avail == 1)
            attitude.q_EKF(1,:)     = attitude.quest(1,:);
        else
            attitude.q_EKF(1,:)     = attitude.est_q(1,:); % do not use quest attitude solution (will be NaN) if one of the vectors is unavailable
        end
        
        [attitude.q_EKF(count+1,:),P_new] = EKF(vectors.gyro_noisy(count,:)',attitude.q_EKF(count,:),P_prev,sc.sensors.Gyro_ARW, sc.sensors.FSS_stdev,...
                                           sc.sensors.MAG_stdev,vectors.sun_eci(count,:)',vectors.mag_eci(count,:)',vectors.sun_b_noisy(count,:)',...
                                           vectors.mag_b_noisy(count,:)', step_in_sec);
%         [attitude.q_EKF(count,:),gyro_bias(count,:),omega_est(count,:),P_new,q_pred] = MEKF(vectors.gyro_noisy(count,:)',sc.sensors.Gyro_ARW,gyro_bias(count,:),vectors.sun_b_noisy(count,:)',sc.sensors.FSS_stdev,vectors.mag_b_noisy(count,:)', sc.sensors.MAG_stdev,attitude.q_EKF(count,:),P_prev, ...
%                                           vectors.sun_eci(count,:)',vectors.mag_eci(count,:)',SS_avail,MAG_avail,step_in_sec);   

        P_prev = P_new;    %propagated P and q
        P_matrix(:,:,count) = P_new;    %saving covariance matrix for 3sig error plot
        attitude.est_q(count+1,:) = attitude.q_EKF(count+1,:);
    end

%% Accounting torques
% Computing Disturbance torques for attitude dynamics    
    [disturbance.torque_total(count,:),disturbance.SRP(count),disturbance.aero(count),disturbance.magnetic(count),disturbance.gg(count)] = disturbance_torques(sc.d,sc.m,sc.I,sc.D,sc.a,accel.drag(count,:),accel.SRP(count,:),vectors.mag_b_true(count,:),OP.r_eci(count,:)');

% Getting control torques from ACS algorithm
    [RW_torque(:,count), MTR_torque(:,count)] = ACS(sc.ACSmode,attitude.est_q(count,1:4)',sc.q_d,omega(count,1:3)',sc.omega_d,[0;0;0],sc.gains);
    torque = disturbance.torque_total(count,:)' + RW_torque(:,count) + MTR_torque(:,count);   %total torque
    
%% Attitude Propagation using Runge-Kutta 4th order    
% Quaternion and Omega propagation (all dynamics and kinematics is included in RK4 function)
     [new_omega,new_q] = RK4_attitude(sc.I,omega(count,:)',torque,q(count,2:5)',step_in_sec);
     
     dynamics.omega_magd(count) = norm(new_omega)*(180/pi);   %magnitude of omega in deg/s
     omega(count+1,:) = new_omega';    %update
     dynamics.omega(count,1:3) = new_omega';
     
     q(count+1,2:5) = new_q';          %update
     attitude.quaternion(count,:) = q(count,2:5);


    progressbar(count/count_max)    %update the progress bar - aesthetics!

end                                 %end the loop






%%      Some Post-Processing
for i = 1:count_max     %computing norms of disturbing accelerations
    mag.nonsphEarth(i) = norm(accel.nonsphEarth(i,:));
    mag.drag(i)        = norm(accel.drag(i,:));
    mag.SRP(i)         = norm(accel.SRP(i,:));
    mag.thruster(i)    = norm(accel.thruster(i,:));
    
end


%% Attitude Estimation error
attitude.quest_err_angles = QuatError(attitude.quaternion(2:end,:),attitude.quest); %Euler angle error in QUEST
attitude.err_angles   = QuatError(attitude.quaternion,attitude.est_q);

[X3sig,Y3sig,Z3sig] = ThreeSigmaCovariance(P_matrix);

% %angular error
% figure
% subplot(4,1,1)
% plot(OP.no_orbits,rad2deg(error_vec(:,1)))
% xlabel('orbits')
% ylabel('roll error (deg)')
% ylim([-1 1])
% subplot(4,1,2)
% plot(OP.no_orbits,rad2deg(error_vec(:,2)))
% xlabel('orbits')
% ylabel('pitch error (deg)')
% ylim([-1 1])
% subplot(4,1,3)
% plot(OP.no_orbits,rad2deg(error_vec(:,3)))
% xlabel('orbits')
% ylabel('yaw error (deg)')
% ylim([-1 1])
% subplot(4,1,4)
% plot(OP.no_orbits,eclipse_flag)
% ylabel('eclipse flag')
% 
% %gyro bias
% figure
% subplot(3,1,1)
% plot(OP.no_orbits,rad2deg(error_vec(:,4)))
% xlabel('orbits')
% ylabel('gyro bias X (deg/s)')
% subplot(3,1,2)
% plot(OP.no_orbits,rad2deg(error_vec(:,5)))
% xlabel('orbits')
% ylabel('gyro bias Y (deg/s)')
% subplot(3,1,3)
% plot(OP.no_orbits,rad2deg(error_vec(:,6)))
% xlabel('orbits')
% ylabel('gyro bias Z (deg/s)')
%Sun vector determination error

%B-field vector determination error

%%      Plotting

% ECI position
% figure
% plot(OP.no_orbits,OP.r_eci(2:count_max+1,:),'.')
% title('ECI position')
% xlabel('orbits')
% ylabel('meters')
% legend('r_x','r_y','r_z')

% ECI velocity
% figure
% plot(OP.no_orbits,OP.v_eci(2:count_max+1,:),'.')
% title('ECI velocity')
% xlabel('orbits')
% ylabel('m/s')
% legend('v_x','v_y','v_z')

% Groundtrack
% figure
% plot(OP.lon,OP.lat,'.')
% title('Orbit Groundtrack')
% xlabel('longitude')
% ylabel('latitude')
% grid on
% xlim([-180 180])
% ylim([-90 90])

% Orbit altitude
% figure
% plot(OP.no_orbits,OP.alt/1000,'.')
% title('Orbit Altitude')
% xlabel('orbits')
% ylabel('altitude (km)')
% ylim([100 500])

% Accelerations
% figure
% subplot(4,1,1)
% plot(OP.no_orbits,mag.nonsphEarth)
% title('Magnitude of disturbing accelerations')
% xlabel('orbits')
% ylabel('non-sph Earth m/(s^2)')
% subplot(4,1,2)
% plot(OP.no_orbits,mag.drag)
% xlabel('orbits')
% ylabel('aero drag m/(s^2)')
% subplot(4,1,3)
% plot(OP.no_orbits,mag.SRP)
% xlabel('orbits')
% ylabel('SRP m/(s^2)')
% subplot(4,1,4)
% plot(OP.no_orbits,mag.thruster)
% xlabel('orbits')
% ylabel('thruster m/(s^2)')

% % Environmental Disturbance torques
% figure
% subplot(4,1,1)
% plot(OP.no_orbits,disturbance.aero)
% title('Magnitude of disturbing torques')
% xlabel('orbits')
% ylabel('aero (Nm)')
% subplot(4,1,2)
% plot(OP.no_orbits,disturbance.SRP)
% xlabel('orbits')
% ylabel('SRP (Nm)')
% subplot(4,1,3)
% plot(OP.no_orbits,disturbance.magnetic)
% xlabel('orbits')
% ylabel('magnetic (Nm)')
% subplot(4,1,4)
% plot(OP.no_orbits,disturbance.gg)
% xlabel('orbits')
% ylabel('grav gradient (Nm)')
% 
% % Angular velocity of the SC
% figure
% subplot(2,1,1)
% plot(OP.no_orbits,dynamics.omega*(180/pi))
% legend('w_x','w_y','w_z')
% xlabel('orbits')
% ylabel('Angular velocity (deg/s)')
% subplot(2,1,2)
% plot(OP.no_orbits,dynamics.omega_magd,'.')
% xlabel('orbits')
% ylabel('magnitude (deg/s)')

% Visualize Earth as a sphere and plot orbit
% figure
% plot3(OP.r_eci(:,1),OP.r_eci(:,2),OP.r_eci(:,3),'.')
% hold on
% [x,y,z] = sphere();
% r = 6.378137e6;
% surf(r*x,r*y,r*z);  %plots sphere with radius r centered at (0,0,0)
% grid on
% title('Orbit visualization')
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';

% Attitude error
figure
subplot(2,1,1)
plot(OP.no_orbits,attitude.err_angles(1:end-1,:))
hold on
% plot(OP.no_orbits,ThreeSigmaCovariance(P_matrix),'LineWidth',1,'Color','k')
% plot(OP.no_orbits,-ThreeSigmaCovariance(P_matrix),'LineWidth',1,'Color','k')
xlabel('orbits')
ylabel('degrees')
ylim([-10 10])
legend('X','Y','Z')
title('Attitude error using EKF')
subplot(2,1,2)
plot(OP.no_orbits, attitude.quest_err_angles)
hold on
plot(OP.no_orbits,5.*eclipse_flag,'--','Color','k')
xlabel('orbits')
ylabel('degrees')
legend('X','Y','Z','EclipseFlag')
ylim([-10 10])
title('Attitude error using QUEST alone')

