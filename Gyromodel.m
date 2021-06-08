% Gyroscope model

% Written by Sarthak Srivastava 16-01-2020

% true rate is in rad/s

function [gyro_output] = Gyromodel(true_rate,gyro_drift,gyro_offset,gyro_SNR,timestep)

gyro_output = ((true_rate*(180/pi))' + gyro_drift*timestep+ gyro_offset);

gyro_output = awgn(gyro_output,gyro_SNR)';


end