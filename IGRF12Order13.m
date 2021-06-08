%IGRF Magnetic Field Model 
%geocentric lat input will be converted to geodetic lat by the program

function B_ECI_model = IGRF12Order13(longitude_model, latitude_model, altitude_model, time)
    % Calculate Earth magnetic field IGRF12 order 13
    % Input:    longitude_model: true longitude (-180 +180 degree) 
    %           latitude_model: true latitude (-90 90 degree) geocentric
    %           altitude_model: altitude from Earth surface (km) height
    %           time: time in calendar form (will be converted to decyear
    %           format in the program later)
    %
    % Output:   B_ECI_model: 3 axis earth magnetic field (nTesla) (1 nTesla = 100mGauss)
    
     gd = geoc2geod( latitude_model, 6379136, 'WGS84');
     [B_RTP_model, H, DEC, DIP, F] = igrf12magm(altitude_model*1000, gd, longitude_model, decyear(time));   %geodetic latitude and height in metres
     B_RTP_model = B_RTP_model'/1000;
     B_EF_model = EF2RTP(gd, longitude_model) * B_RTP_model;    %conversion from NED to ECEF, input geodetic latitude (source: http://www.mathworks.com/help/aeroblks/directioncosinematrixeceftoned.html)
     %old R_ECI2EF = EarthRot(JD2T(Date2JD(time)))';
     julian = juliandate(time); %intermediate step to convert to Julian Century
     R_ECI2EF = ECIToEF(JD2T(julian));
     B_ECI_model = R_ECI2EF * B_EF_model;
end

function R_EF2RTP = EF2RTP(latitude, longitude)
    % Calculate rotation matrix from EF frame to RTP frame
    % RTP frame : +z toward Earth nadir
    %             +y toward east
    %             +x completes orthogonal frame
    % Input:    longitude: (0-360 degree) 
    %           geodetic latitude: (-90 90 degree) 
    % Output: Rotational matrix R
    z_RTP = - [cosd(latitude) * cosd(longitude); cosd(latitude) * sind(longitude); sind(latitude)];
    y_RTP = cross(z_RTP, [0;0;1]);
    y_RTP = y_RTP/norm(y_RTP);
    x_RTP = cross(y_RTP, z_RTP);
    x_RTP = x_RTP/norm(x_RTP);
    R_EF2RTP = [x_RTP y_RTP z_RTP];
end