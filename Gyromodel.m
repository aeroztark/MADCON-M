function [gyro_output] = Gyromodel(true_rate,gyro_drift,gyro_offset,gyro_SNR,timestep)

% Gyroscope model
% true rate is in rad/s

gyro_output = ((true_rate*(180/pi))' + gyro_drift*timestep+ gyro_offset);

gyro_output = awgn(gyro_output,gyro_SNR)';


end