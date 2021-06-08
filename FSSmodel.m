%FSS model

%created by Sarthak Srivastava 16-01-2020

function [FSS_sunvector] = FSSmodel(sunb_true,max_err_angle,FSS_SNR,eclipse_flag)

%model: FSS output = (true body sunvector + rand error rotation) + AWG

if eclipse_flag == 1
    FSS_sunvector = [NaN, NaN, NaN];
    
else

    %randomize error angle magnitude
    err_angle     = max_err_angle*(pi/180)*rand(1,1);

    %randomize error rotation vector
    axis          = rand(3,1);
    axis          = axis/norm(axis);

    %create error quaterion
    q_FSS_error   = AU2Q(err_angle,axis);

    %superimpose error to the true body angle
    FSS_sunvector = quatrotate(q_FSS_error',sunb_true);
    FSS_sunvector = awgn(FSS_sunvector,FSS_SNR);
%end
end



