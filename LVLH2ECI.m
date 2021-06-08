function[DCM] = LVLH2ECI(r_eci,v_eci)
% This function generates the DCM from LVLH frame to ECI frame

%Ref: Markley-Crassidis section 2.6.4

crossp = cross(r_eci,v_eci);

o3 = -r_eci/norm(r_eci);
o2 = -crossp/norm(crossp);
o1 = cross(o2,o3);

DCM = [o1 o2 o3];

% q = dcm2quat(DCM);
% 
% q = q/norm(q);

end
