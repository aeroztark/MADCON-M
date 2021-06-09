function [a_drag] = dragmodel(r_eci,v_eci,sc_area,sc_Cd,sc_mass,rho,DCM)

% Computes acceleration due to atmospheric drag

    % INPUTS: r_eci,v_eci -> vectors as 1x3 or 3x1 in SI units
    %         sc_area in m^2, sc_Cd, sc_mass -> from param() file
    %         rho -> atmospheric density from a suitable density model
    %         q -> current attitude quaternion as 1x4

    %   Relative veocity of the SC w.r.t the atmosphere
    w_earth = 0.000072921158553;    %rad/s
    
    v_relB_eci = [v_eci(1) + w_earth*r_eci(2);
                  v_eci(2) - w_earth*r_eci(1);
                  v_eci(3)];
    
    %   In body frame
    %v_relB = quatrotate(q,v_relB_eci'); %inputs to quatrotate must be 1x4 and 1x3
    v_relB = DCM*v_relB_eci;
    %   Components of the unit vector are direction cosines
    v_relB_unit = v_relB/norm(v_relB);
    
    k = (0.5*rho*sc_Cd*norm(v_relB));
    
    F_drag = -k*[v_relB(1)*sc_area(1)*v_relB_unit(1);
                 v_relB(2)*sc_area(2)*v_relB_unit(2);
                 v_relB(3)*sc_area(3)*v_relB_unit(3);];
             
    a_drag = F_drag/sc_mass;
    
end
    
             