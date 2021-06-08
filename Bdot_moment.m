
%reference: Magnetic satellite detumbling: the b-dot algorithm revisited,
% by Marco Lovera

function [control_m] = Bdot_moment(magb_noisy,ang_velocity,K,max_dipole)

B_dot = cross(magb_noisy,ang_velocity); %derovative of body frame mag vector

control_m = -K*max_dipole*sign(B_dot);  %simple controller


%Need to ensure that the control_m magnitude does not exceed MTR physical
%specs!
