% Another sun model based on algorithm in Markeley, Crassidis. 
%This sun model is simpler and yield results similar to the elaborate on being used currently
% Kept as backup


function[u_SunEarth] = sunmodel2(jdate)
%input: juliandate



T_UT1 = (jdate - 2451545)/36525;

phi_sun = mod(280.460 + 36000.771*T_UT1,360);    %mean longitude of the Sun in 0-360 deg range
M_sun   = mod(357.5277233 + 35999.05034*T_UT1,360);  %mean anomaly of the Sun in 0-360 deg range

phi_ecliptic = phi_sun + 1.914666471*sind(M_sun) + 0.019994643*sind(2*M_sun);   %longitude of ecliptic in degrees
e = 23.439291 - 0.0130042*T_UT1;    %obliquity of ecliptic in deg

u_SunEarth = [cosd(phi_ecliptic);   %ECI sun unit vector
              cosd(e)*sind(phi_ecliptic);
              sind(e)*sind(phi_ecliptic);]'; 
          
          






