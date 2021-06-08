function [offsetval] = GyroOffset(trueval) %% Takes in true value and outputs offset at that value
offsetval = 0.02704*trueval+ 0.04086   ; %% Linear Fit Model for offset as a function of true velocity obtained from sample data 
end
