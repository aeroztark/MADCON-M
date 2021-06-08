function[vdot] = vdot(r,v,A,Cd,m,rho,DCM,eclipse,sunb)
    global a_point a_nonsphEarth a_drag a_SRP a_thruster
    [a_point,a_nonsphEarth] = gravitymodel(r);
    a_gravity = a_point + a_nonsphEarth;
    a_drag = dragmodel(r,v,A,Cd,m,rho,DCM);
    a_SRP = SRPmodel(eclipse,sunb,A,m);
    a_thruster = ThrusterModel(r,v,m,DCM);
    
    vdot = a_gravity + a_drag + a_SRP + a_thruster;
end