function[ThreeSig_x, ThreeSig_y, ThreeSig_z] = ThreeSigmaCovariance(P_matrix_array)
% This function computes the 3-sigma level from stored covariance matrices.
% After computation, + and - 3-sigma levels can be plotted

% INPUTS: 3D (3x3xtime dimension) P_matrix_array
% OUTPUTS: 3x1 vector of x, y and z 3-sigma levels (in degrees)

% compute square root of terms (to get tdev from variance)
P_sqrt = sqrt(P_matrix_array);

% compute 3 sig level and convert to deg
ThreeSig_x = (180/pi)*3*squeeze(P_sqrt(1,1,:));
ThreeSig_y = (180/pi)*3*squeeze(P_sqrt(2,2,:));
ThreeSig_z = (180/pi)*3*squeeze(P_sqrt(3,3,:));

end