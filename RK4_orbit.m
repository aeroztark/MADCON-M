
function [updated_r,updated_v] = RK4_orbit(r_old,v_old,dt,A,Cd,m,rho,DCM,eclipse,sunb)
% RK4 integration for orbit propagation. 
% Input: r_old and v_old (ECI frame 3x1 vectors) in SI units, dt in sec
% Output: updated_r and updated_v are 3x1 vectors

        %first propagate v and r alternately OR can be written in matrix
        %form
        k1v = vdot(r_old,v_old,A,Cd,m,rho,DCM,eclipse,sunb);
        k1r = rdot(v_old);
        
        k2v = vdot(r_old+0.5*k1r*dt,v_old+0.5*k1v*dt,A,Cd,m,rho,DCM,eclipse,sunb);
        k2r = rdot(v_old+0.5*k1v*dt);
        
        k3v = vdot(r_old+0.5*k2r*dt,v_old+0.5*k2v*dt,A,Cd,m,rho,DCM,eclipse,sunb);
        k3r = rdot(v_old+0.5*k2v*dt);
        
        k4v = vdot(r_old+k3r*dt,v_old+k3v*dt,A,Cd,m,rho,DCM,eclipse,sunb);
        k4r = rdot(v_old+k3v*dt);
      
                
        updated_v = v_old + (1/6)*dt*(k1v + 2*k2v + 2*k3v + k4v);
        updated_r = r_old + (1/6)*dt*(k1r + 2*k2r + 2*k3r + k4r);
end 
