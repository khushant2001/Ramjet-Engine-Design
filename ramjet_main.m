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

%% Inlet condtions!

elevation = 16764; %m
mach_in = 3.25; 
[T1, P1] = atmospheric(elevation);  

%% Imposed conditions!

%T3 = 1500; % K. To prevent high temperatures! Need to talk about this
T_tungsten = 1900; %K. Melting point of tungsten!
%% Algorithm implementation!

% Inlet
[T2, P2, M2,x_cowl,y_cowl,x,y,inlet_area] = inlet_design(3.25,P1,T1,5,1);

% Diffuser
[x_diffuser,A_diffuser,M3,P3,T3] = diffuser(M2,T2,P2,inlet_area/2,T1,mach_in);
x_diffuser = x_diffuser + x(end);
A_diffuser = A_diffuser + (y_cowl(2) - A_diffuser(1));
% Flameholder
P_3_prime = flameholder(P3(end),M3);
% Combustor

% Converging section

% Nozzle

%% Plotting the engine!
figure
hold on
plot(x, y, 'ro-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_cowl,y_cowl,'bo-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_diffuser(1:4), A_diffuser(1:4), 'go-', 'MarkerSize', 10, 'LineWidth', 2)
plot(x_diffuser(5:end), A_diffuser(5:end), 'go-', 'MarkerSize', 10, 'LineWidth', 2)





