function [lat,lon,alt] = ECI2LLA(r_eci,jd_time)
    %This function computes LLA from ECI coordinates
    
    %INPUT: r_eci in 3x1
    %       time in juliandate
    
    %OUPUT: lat,lon in deg and alt in m
    
    % Author - Sarthak Srivastava (14-03-2019)
    
    Re = 6.378137e6;    %in metres
    % Get julian centuries
    jcenturies = JD2T(jd_time);
    
    % Get DCM from ECI to ECEF
    DCM_eci2ecef = ECIToEF(jcenturies);
    
    %ECEF coordinates
    r_ecef = DCM_eci2ecef*r_eci;    %ensure that r_eci is 3x1
    
    %ECEF to LLA
    l = r_ecef(1)/norm(r_ecef);
    m = r_ecef(2)/norm(r_ecef);
    n = r_ecef(3)/norm(r_ecef);
    
    lat = asind(n);
    lon = atan2d(m,l);
    alt = (norm(r_eci) - Re);
    
end
    