%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ramjet Engine Design %
    % Angel, Christos, Khushant %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

% Thing to think about!
% 1. Should have an x and y array passing around the functions to get the
% location of every single component. 
clc
clear 
close all

%% Inlet condtions
elevation = 16764; %m
mach_in = 3.25; 
[T_1, P_1] = atmospheric(elevation);

%% Algorithm implementation!

[T_2, P_2, M_2,x_cowl,y_cowl,x,y,inlet_area] = inlet_design(3.25,P_1,T_1,1,2);
[x_diffuser,A_diffuser] = diffuser(M_2,T_2,P_2,inlet_area);

figure
hold on
scatter(x,y)
scatter(x_cowl,y_cowl)
scatter(x_diffuser,A_diffuser)
