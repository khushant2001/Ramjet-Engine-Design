% The following is the code for the adiabatic pressure loss of the
% flameholder section. 

function [P3_prime,T3_prime] = flameholder(P_3, M_3,T_3)
    disp("Calculating Properties Across Flameholder ...")
    R = 287;
    gaama = 1.4;    
    P3_prime = P_3 - P_3*.81*gaama*M_3^2;
    T3_prime = P3_prime*sqrt(T_3)/P_3;
end