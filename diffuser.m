% This is for subsonic flow. So we can consider the fluid to be
% incompressible. And use those estimations for our stuff. 

% Things to take care of!
% 1. What is unstart!
% 2. What is the limit on the temperature for.
% 3. How do we find the x coordinates. 
% 4. How do we include the temperature limits. 
function [x,A2] = diffuser(M_in, T_in, P_in,A_in)
    
    % Defininng initial properties!
    gaama = 1.4;
    M2 = M_in;
    T_0 = T_in*(1 + .5*(gaama-1)*M_in^2);
    P_0 = P_in*(1 + .5*(gaama-1)*M_in^2)^(gaama/(gaama-1));

    % Variables the need to be deleted later!!
    A1 = A_in; % A1 = A_in
    i = 0; % Need to figure out how to get the x points. So this can be deleted. 

    % Defining arrays to store future state properties. 
    A2 = [];
    T2 = [];
    P2 = [];
    M_out = [];

    while M2 >= 0.4 % This .4 is arbritary. TODO
        M_out = [M_out,M2];
        i = i+1;
        M2 = M2 - .005;
        area_ratio = (M2/M_in)*((1+.5*(gaama-1)*M_in^2)/(1+.5*(gaama-1)*M2^2))^((gaama+1)/(2*(gaama-1)));
        A2 = [A2,A1/area_ratio];
        T2 = [T2,T_0/(1+.5*(gaama-1)*M2^2)];
        P2 = [P2,P_0/((1+.5*(gaama-1)*M2^2)^(gaama/(gaama-1)))];
    end
    x = linspace(0,i,i);
end