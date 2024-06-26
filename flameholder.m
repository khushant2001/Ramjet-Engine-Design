% The following is the code for the adiabatic pressure loss of the
% flameholder section. 

function [P3_prime,T3_prime,M3_prime] = flameholder(P3, M3,T3)
    %disp("Calculating Properties Across Flameholder ...")
    
    % Declaring some constants!
    R = 287;
    gaama = 1.4; 
    cp = gaama*R/(gaama-1);

    P3_prime = P3 - P3*.81*gaama*M3^2;
    tolerance = 0.00001;
    exp1 = sqrt(2*cp*T3 + M3^2*gaama*R*T3)*P3*M3/sqrt(T3);
    M3_prime = 0;
    while 1
        exp2 = P3_prime*M3_prime*sqrt(2*cp + M3_prime^2*gaama*R);
        if exp1 - exp2 < tolerance
            break
        else
            M3_prime = M3_prime + 0.001;
        end
    end
    T3_prime = (P3_prime^2*M3_prime^2*T3)/(P3^2*M3^2);
    
end