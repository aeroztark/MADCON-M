%{

Convert a set of Orbital Elements (a,e,i,RAAN,w,v) to a State Vector (R,V)

Reference: Orbital Mechanics for Engineering Students by H. Curtis

%}
function [R_eci, V_eci] = oe2rv(orbital_elements) %returns 3x1 vectors r and v
%{
INPUTS: pay attention to the order!
    'a'     Semimajor Axis [m]
    'e'     Eccentricity Magnitude
    'i'     Inclination [deg]
    'RAAN'  Right Ascension of the Ascending Node [deg]
    'w'     Argument of Periapse [deg]
    'v'     True Anomaly [deg]

OUTPUTS:
R: Radius Vector in ECI [m]
V: Velocity Vector in ECI [m/s]
%}

%% Extracting input values
a = orbital_elements(1);
e = orbital_elements(2);
i = orbital_elements(3);
RAAN = orbital_elements(4);
w = orbital_elements(5);
v = orbital_elements(6);

%% Universal Variables
r_earth = 6378137; %m
u_earth = 3.986e14; %m^3/s^2

%Inertial Axes
I = [1; 0; 0];
J = [0; 1; 0];
K = [0; 0; 1];

%% Rotation Matrices
R1 = @(j) [1  0      0
           0  cos(j) sin(j)
           0 -sin(j) cos(j)];

R2 = @(j) [cos(j) 0 -sin(j)
           0      1  0
           sin(j) 0  cos(j)];

R3 = @(j) [ cos(j) sin(j) 0
           -sin(j) cos(j) 0
            0      0      1];
%% Convert All Values to Radians
i = (pi/180)*i;
RAAN = (pi/180)*RAAN;
w = (pi/180)*w;
v = (pi/180)*v;

%% Radius & Velocity in the Perifocal Frame [PQW]
R_pqw = [(a*(1-e^2))/(1+e*cos(v))*cos(v)
         (a*(1-e^2))/(1+e*cos(v))*sin(v)
          0];
     
V_pqw = [-sqrt(u_earth/(a*(1-e^2)))*sin(v)
          sqrt(u_earth/(a*(1-e^2)))*(cos(v)+e)
          0];

%% Rotation of Perifocal into ECI Frame
R_eci = (R3(-RAAN)*R1(-i)*R3(-w)*R_pqw);   
V_eci = (R3(-RAAN)*R1(-i)*R3(-w)*V_pqw);


