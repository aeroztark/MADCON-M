% Magnetometer model 

%Written by Sarthak Srivastava, 16-01-2020


function [magnetometer_vector] = MAGmodel(mag_b_true,mag_SNR)

mag_noisy           = awgn(mag_b_true*(1e6),mag_SNR);

magnetometer_vector = mag_noisy*(1e-6);

end