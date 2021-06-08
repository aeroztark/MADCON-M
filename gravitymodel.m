
function[a_point,a_nonsphEarth] = gravitymodel(r_eci)
% This function computes ECI acceleration on the SC due to gravity of the...
  ...Earth. Spherical Harmonic coefficients upto J4 are modelled.

% INPUT: r_eci as 3x1 vector in metres
% OUTPUT: a_gravity as 3x1 vector

	mu = 3.986004418e14;
	Re = 6.378137e6;
    r = norm(r_eci);
    x = r_eci(1);
    y = r_eci(2);
    z = r_eci(3);
    
    %Zonal Spherical Harmonic coefficients
    J2 = +1.08262668355e-3;
    J3 = -2.53265648533e-6;
    J4 = -1.61962159137e-6;
    
	% Point mass for Earth as central body (Ideal Keplerian 2-Body case)
    a_point = (-mu/(r^3)).*r_eci;
    
	%Perturbations due to non-Spherical Earth
    a_J2 = (-1.5*J2*mu/(r^2))*((Re/r)^2).*[(1-5*((z/r)^2))*(x/r);
                                           (1-5*((z/r)^2))*(y/r);
                                           (3-5*((z/r)^2))*(z/r)];
                                       
    a_J3 = (-0.5*J3*mu/(r^2))*((Re/r)^3).*[(5*(7*((z/r)^3)-3*(z/r))*(x/r));
                                           (5*(7*((z/r)^3)-3*(z/r))*(y/r));
                                           (3*(10*((z/r)^2)-(35/3)*((z/r)^4)-1))];
                                       
    a_J4 = (-(5/8)*J4*mu/(r^2))*((Re/r)^4).*[(3-(42*(z/r)^2)+63*((z/r)^4))*(x/r);
                                             (3-(42*(z/r)^2)+63*((z/r)^4))*(y/r);
                                             -(15-(70*(z/r)^2)+63*((z/r)^4))*(z/r)];
                                       
    %Total acceleration due to Earth's gravity
    a_nonsphEarth = a_J2 + a_J3 + a_J4;
end