% The following is the code for the adiabatic pressure loss of the
% flameholder section. 

function [P_3_prime] = flameholder(P_3, M_3)
    disp("Calculating Properties Across Flameholder ...")
    gaama = 1.4;    
    P_3_prime = P_3 - P_3*.81*gaama*M_3^2;
end