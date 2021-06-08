function[a_SRP] = SRPmodel(eclipseflag,sun_body,sc_area,sc_mass)

%Solar Radiation Pressure basic model (elaborate further based on Markley-Crassidis)
phi = 1366;
c = 3e8;
q = 0.6;
mu = 3.986e14;  %mu of Earth
k = (phi/c)*(1+q);
    if eclipseflag == 1
        F_SRP = [0;0;0];
    else
        F_SRP = k.*[sc_area(1)*sun_body(1); sc_area(2)*sun_body(2);sc_area(3)*sun_body(3)];        
    end

a_SRP = F_SRP/sc_mass;

end
