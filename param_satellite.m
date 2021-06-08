%%Parameter file - returns a struct sc
%%Enter all parameters for the satellite here and call it in the main model
%%Create a unique param file for each satellite


%%  
function [sc] = param_satellite()

%----------Satellite Paramaters-----------------------------
sc.ID = 'SCOOB-2';

%mass
sc.m = 3;      %kg

%dimensions
sc.area_deployed = [0.3*0.1;0.1*0.1;0.3*0.3];  %[Ax;Ay;Az] area of yz,xz and xy faces respectively
sc.d = [0.001;0.005;-0.002];		% [dx,dy,dz] - (cp-cm) distance in x,y,z resepectively

%MoI (deployed)
sc.I = [0.029   0     0;   %kg-m^2
          0   0.007   0;
          0     0   0.03];
   
%Initial Keplerian Elements
sc.a              = 6928000;    %semi-major axis in m (~ 550 km)
sc.e              = 0.001;       %eccentricity (non-zero to avoid singularity)
sc.raan           = 30;  %RAAN in deg
sc.i              = 5;      %inclination in deg
sc.ap             = 5;      %argument of periapsis in deg
sc.ta             = 5;      %true anomaly in deg
sc.oe             = [sc.a sc.e sc.i sc.raan sc.ap sc.ta];

%---------Simulation Parameters-----------------------------
sc.start_time     = [2022 11 10 1 1 1];  %Start Epoch
sc.end_time       = [2022 11 10 5 1 1]; 
sc.period         = 2*pi*sqrt(((sc.a)^3)/3.986e14);
sc.timestep       = 0.1;	%in seconds
sc.q_init         = [0.5;0.5;0.5;0.5];	%initial quaternion
sc.omega_init_deg = [0.005;0.001;0.002];	%deg/s for deployed

% CLARIFY DCM implementation! It is reflecting LVLH conditions only!! ISolate it from quaternion%%
%--------Disturbance model parameters------------------------
sc.D              = [0.005;0.005;0.005];	%residual dipole moment (Am^2) of the spacecraft
sc.Cd             = 2.2;	%drag coefficient


%% Detumble parameters
sc.is_detumble    = 0;	%set flag to 1 if doing detumble, 0 if not
sc.BDot_gain      = 0.55;	
sc.max_dipole     = 0.2;	%max dipole moment of torquer (Am^2)

%% Sensor model parameters
sc.sensors.FSSmax_err  = 0.5;  %3-sig FSS err in degrees
sc.sensors.FSS_SNR     = 50;
sc.sensors.FSS_stdev   = 0.1; % deg
sc.sensors.MAG_SNR     = 50; 
sc.sensors.MAG_stdev   = 0.125e-2; % T
sc.sensors.Gyro_drift  = [0.01;0.02;0.015]./10; %deg/s^2
sc.sensors.Gyro_offset = [0.2;0.3;-0.5]./10;  %deg/s
sc.sensors.Gyro_ARW    = 0.1;      %deg/s^0.5
sc.sensors.Gyro_SNR    = 30;

%% ADCS algorithm paramaters
sc.use_EKF = 1;

%% ACS parameters
sc.ACSmode = 'TargetPointing';    

J_RW = [0.00000820  0   0;
        0   0.00000820  0;
        0   0   0.00001281];                        % Moment of inertia of reaction wheel kg-m^2
rw_align = eye(3);                                  % Reaction wheels alignment matrix
k_n = 543;                                          % BLDC speed constant in min^-1/V
rw_conversion = k_n*3.141592/30;                    % rad/s/V conversion factor

% Controller gains
sc.gains.K_p = 0.2;                                          % Proportional gain; 0.2 for quaternion control;
sc.gains.K_d = 0.1;                                          % Differential gain; 0.1 for quaternion control; (-) depending on error calc.
sc.gains.K_p_w = 0.7;                                        %  0.7 for omega control 
sc.gains.K_d_w = -0.005;                                     % -0.005 for omega control

% Initial Conditions
omega_rw_initial = [0 0 0]';                        % Initial angular velocity of RWs
sc.omega_d = [0;0;0];   % desired omega of the satellite (rad/s)
sc.q_d = [1; 0; 0; 0]; % desired quaternion
