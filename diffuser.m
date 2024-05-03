% This is for subsonic flow and we are expanding the flow! 

% Things to take care of!
% 1. What is unstart!

function [x,A2,M_req,P3,T3] = diffuser(M_in, T_in, P_in,A_in,T_inf,M_inf)
    
    disp("Calculating Properties Across Diffuser ...")
    % Defininng initial properties!
    gaama = 1.4;
    M2 = M_in;

    % Finding stagnation values at the inlet. 
    T_0 = T_in*(1 + .5*(gaama-1)*M_in^2);
    P_0 = P_in*(1 + .5*(gaama-1)*M_in^2)^(gaama/(gaama-1));
    
    % Finding M3_final, at exit of diffuser, to maintain isentropic flow!
    M_req = .2;
    T3_final = T_inf/((M_req^2)/(2/(gaama-1))+1)/(1+.5*(gaama-1)*M_inf^2);
    % Defining arrays to store future state properties. 
    A2 = [A_in];
    T3 = [];
    P3 = [];

    if imag(M_req) ~= 0
        disp("...ERROR! DIFFUSER ERROR. CAN'T ATTAIN THE REQUIRED TEMPERATURE!")
    end

    M_out = [];
    
    if M2 <= M_req
        disp("ERROR: DIFFUSER CANT INCREASE FLOW SPEED!!");
    else
        while M2 >= M_req % This NEEDS TO BE FIXED!
            M_out = [M_out,M2];
            i = i+1;
            M2 = M2 - .0005;
            area_ratio = (M2/M_in)*((1+.5*(gaama-1)*M_in^2)/(1+.5*(gaama-1)*M2^2))^((gaama+1)/(2*(gaama-1)));
            A2 = [A2,A_in/area_ratio];
            T3 = [T3,T_0/(1+.5*(gaama-1)*M2^2)];
            P3 = [P3,P_0/((1+.5*(gaama-1)*M2^2)^(gaama/(gaama-1)))];
        end
    end
    x = linspace(0,2*(A2(end)-A2(1)),length(A2));
    x = [x,x];
    A2 = [A2,-A2];
end