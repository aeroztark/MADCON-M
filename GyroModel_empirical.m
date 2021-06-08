function [measurementvalue] = GyroModel_empirical(truevalue) %% Takes in true angular velocity and outputs a measurement value
measurementvalue= truevalue+ GyroOffset(truevalue)+ GyroNoise(truevalue) ; %% Measurement value is a function of the following 
end
