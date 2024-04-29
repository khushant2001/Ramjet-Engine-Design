% This is for subsonic flow. So we can consider the fluid to be
% incompressible. And use those estimations for our stuff. 

% Things to take care of!
% 1. What is unstart!
% 2. What is the limit on the temperature for.
% 3. How do we find the x coordinates. 
% 4. How do we include the temperature limits. 
function [x,A2,M_req,P3,T3] = diffuser(M_in, T_in, P_in,A_in,T_inf,M_inf)
    
    disp("Calculating Properties Across Diffuser ...")
    % Defininng initial properties!
    gaama = 1.4;
    M2 = M_in;

    % Finding stagnation values at the inlet. 
    T_0 = T_in*(1 + .5*(gaama-1)*M_in^2);
    P_0 = P_in*(1 + .5*(gaama-1)*M_in^2)^(gaama/(gaama-1));

    % Variables the need to be deleted later!
    i = 0; % Need to figure out how to get the x points. So this can be deleted. 
    
    % Finding T3, at exit of diffuser, to maintain isentropic flow!
    T3_final = 670; %T_inf*(1 + .5*(gaama-1)*M_inf^2); % This is to completely stop the flow!

    % Defining arrays to store future state properties. 
    A2 = [A_in];
    T3 = [];
    P3 = [];
    M_req = sqrt((2/(gaama-1))*((T_inf/T3_final)*(1+.5*(gaama-1)*M_inf^2)-1));
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
            M2 = M2 - .005;
            area_ratio = (M2/M_in)*((1+.5*(gaama-1)*M_in^2)/(1+.5*(gaama-1)*M2^2))^((gaama+1)/(2*(gaama-1)));
            A2 = [A2,A_in/area_ratio];
            T3 = [T3,T_0/(1+.5*(gaama-1)*M2^2)];
            P3 = [P3,P_0/((1+.5*(gaama-1)*M2^2)^(gaama/(gaama-1)))];
        end
    end
    x = linspace(0,1*(A2(end)-A2(1)),length(A2));
    x = [x,x];
    A2 = [A2,-A2];
end