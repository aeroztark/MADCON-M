
function[total_torque,norm_SRP,norm_aero,norm_magnetic,norm_gg,norm_thruster] = disturbance_torques(sc_d,mass,sc_I,sc_D,sc_a,a_drag,a_SRP,a_thruster,magb_vector,r_eci)


T_aero = cross(sc_d,mass*a_drag);
norm_aero = norm(T_aero);
 
T_SRP = cross(sc_d,mass*a_SRP);
norm_SRP = norm(T_SRP);

T_magnetic = cross(sc_D,magb_vector);
norm_magnetic = norm(T_magnetic);

%Gravity gradient
mu = 3.986e14;
w = sqrt(mu/(sc_a^3)); %Orbital angular velocity
r_eci_unit = r_eci/norm(r_eci);
product = sc_I*r_eci_unit;
T_gg = (3*w^2)*cross(r_eci_unit,product);
T_gg = T_gg';
norm_gg = norm(T_gg);

% if thrusterflag == 1 && thruster is not exhausted
% T_thruster(count,:) = 
% else thrust = 0
% T_thruster = cross([0 0.005 0],mass*a_thruster);

T_thruster = [0;0;0]'; %no firing
norm_thruster = norm(T_thruster);
total_torque = T_aero + T_SRP + T_magnetic + T_gg + T_thruster;

end