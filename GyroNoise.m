function [gaussiannoise] = GyroNoise(trueval) %% Takes in true value and outputs gaussian random variable for noise
variance = (trueval*  0.003688)  +0.002719 ; %% Variance for normal distribution as function of true velocity obtained from sample data
deviation = sqrt(variance) ; 
gaussiannoise = normrnd(0,deviation*1.2) ; %% Generate gaussian variable with safety factor of 1.2
end
