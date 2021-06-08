function[eclipse_flag] = EclipseCheck(sun_eci_vector, SC_position_vector)

%Eclipse check function [ref: Markeley & Crassidis]
%This function determines shadowing of the cylindrical projection of
%Earth's diameter on the spacecraft - adequate for LEO

%inputs
% sun_eci_vector: 3x1 sun unit vector in ECI frame
% SC_position_vector: 3x1 spacecraft ECI position vector in km

%output:
%eclipse_flag: 0 if no eclipse
%               1 if eclipse

R_earth = 6378; %km

if dot(SC_position_vector,sun_eci_vector) < -sqrt((norm(SC_position_vector))^2 - R_earth^2)
    eclipse_flag = 1;
else
    eclipse_flag = 0;
end

end