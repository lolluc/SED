function [Rot_ge_pf,Rot_pf_ge] = rot (OMEGA,omega,i)


    Rot_OMEGA = [cos(OMEGA) sin(OMEGA) 0; -sin(OMEGA) cos(OMEGA) 0; 0 0 1];
    Rot_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    Rot_omega = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    Rot_ge_pf = Rot_omega*Rot_i*Rot_OMEGA;
    Rot_pf_ge = Rot_ge_pf';

end