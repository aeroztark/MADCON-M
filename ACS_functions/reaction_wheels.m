
% Reaction Wheel model

function [H_rw_new,omega_rw_new_RPM] = reaction_wheels(J_RW,rw_align,omega_rw,T_control,dt)
  omega_rw_new = omega_rw + omega_rw_dot(J_RW,T_control)*dt;     %
  H_rw_new = rw_align*J_RW*omega_rw_new;                        
  omega_rw_new_RPM = omega_rw_new/(2*pi);
end

function omega_rw_dot = omega_rw_dot(J_RW,T)
omega_rw_dot = J_RW\T;
end