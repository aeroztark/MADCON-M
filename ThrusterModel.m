function[a_thruster] = ThrusterModel(r,v,sat_mass,DCM)
    
    global a_thruster
    
    F_nominal = 0.0006;  %N
    total_impulse = 2500;   %Ns
    
    %v_eci_unit = v_eci/norm(v_eci); %get the unit vector for v_eci
    
%     a_thruster = F_nominal/sat_mass;
%     a_thruster = a_thruster*(-v_eci_unit);   %this is incorrect 

     a_lvlh = -[(F_nominal/sat_mass);0;0];   %thrust along -X in LVLH (to keep thrust tangential)
     %DCM = LVLH2ECI(r,v); %convert acceleration to ECI frame
     %a_thruster = DCM*a_lvlh;
     a_thruster = [0;0;0];  %no firing
end